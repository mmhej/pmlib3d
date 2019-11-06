!---------------------------------------------------------------------------------!
! pmlib_topology_cuboid.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_topology_cuboid(patch,topo,ierr)

USE pmlib_mod_patch

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch),INTENT(INOUT)                         :: patch
  TYPE(class_topology),DIMENSION(:),POINTER,INTENT(OUT)   :: topo
  INTEGER, INTENT(OUT)                                    :: ierr

!---------------------------------------------------------------------------------!
!  Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=21)                     :: caller = 'pmlib_topology_cuboid'
  INTEGER                               :: iproc
  INTEGER                               :: isub,jsub,ksub
  INTEGER                               :: i,j,k

  INTEGER                               :: nsub_max
  INTEGER, DIMENSION(pmlib_ndim)        :: nsub

  INTEGER, DIMENSION(pmlib_ndim)        :: ncell, icell
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin, xmax, dx 
  INTEGER, DIMENSION(pmlib_ndim)        :: bc_ext
  INTEGER, DIMENSION(2*pmlib_ndim)      :: bc, nghost
  INTEGER, DIMENSION(2*pmlib_ndim)      :: ineigh
  INTEGER                               :: nghost_ext, nghost_int

  INTEGER, DIMENSION(pmlib_ndim)        :: rem_ncell

  REAL(MK), DIMENSION(:,:,:), POINTER   :: sub2proc

  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_sub
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin_sub, xmax_sub

!---------------------------------------------------------------------------------!
! Allocate 
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(topo) ) THEN
    DEALLOCATE(topo,stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to deallocate topology.')
      GOTO 9999
    ENDIF
  END IF
  ALLOCATE( topo(0:mpi_nproc-1) )


!---------------------------------------------------------------------------------!
! Exterior boundary conditions
!---------------------------------------------------------------------------------!
  bc_ext = patch%bound_cond

!---------------------------------------------------------------------------------!
! Determine topology parameters
!---------------------------------------------------------------------------------!
  dx   = patch%dx
  xmin = patch%xmin
  xmax = patch%xmax

!---------------------------------------------------------------------------------!
! Number of mesh and ghost cells
!---------------------------------------------------------------------------------!
  nghost_ext = patch%nghost
  IF(pmlib_interpolation_order .EQ. 3)THEN
    nghost_int = 2
  ELSEIF(pmlib_interpolation_order .EQ. 4)THEN
    nghost_int = 3
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'Can not define number of ghost cells.')
    GOTO 9999    
  END IF
  ncell = patch%ncell

!---------------------------------------------------------------------------------!
! Pass patch related data to topology
!---------------------------------------------------------------------------------!
  DO iproc = 0,mpi_nproc-1
    topo(iproc)%ptype  = patch%ptype
    topo(iproc)%level  = patch%level
    topo(iproc)%parent = patch%parent
  END DO

!---------------------------------------------------------------------------------!
! The topology is calculated on processor 0 only after wich it is 
! broadcasted to the remaining processors. This insures an identical 
! truncation and defination of subdomain boundaries.
!---------------------------------------------------------------------------------!
  IF(mpi_rank .EQ. 0)THEN 
!---------------------------------------------------------------------------------!
! Find number of subdomains in each direction
!---------------------------------------------------------------------------------!
    nsub_max = mpi_nproc

    DO i = mpi_nproc,1,-1
      DO j = 1,mpi_nproc
        DO k = 1,mpi_nproc
          IF( i*j*k .EQ. mpi_nproc .AND. & 
            & MAXVAL( (/ i , j , k /) ) .LE. nsub_max)THEN
            nsub = (/ i , j , k /) 
            nsub_max = MAXVAL(nsub)
          END IF

        END DO
      END DO
    END DO

!---------------------------------------------------------------------------------!
! Check the total number of subdomains
!---------------------------------------------------------------------------------!
    patch%nsubs = nsub(1)*nsub(2)*nsub(3)
    IF(patch%nsubs .NE. mpi_nproc)THEN
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Wrong number of subdomains.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Assign subdomains and keep track of the local index of the subdomain on
! the processor
!---------------------------------------------------------------------------------!
    ALLOCATE(sub2proc( nsub(1),nsub(2),nsub(3)) )
    sub2proc = 0
    iproc = 0
    DO ksub = 1,nsub(3)
      DO jsub = 1,nsub(2)
        DO isub = 1,nsub(1)
          sub2proc(isub,jsub,ksub) = iproc
          iproc = iproc + 1
        END DO
      END DO
    END DO

!---------------------------------------------------------------------------------!
! Find extend of each subdomain and number of grid points
!---------------------------------------------------------------------------------!
    icell(3) = 1
    xmin_sub(3) = xmin(3)
    DO ksub = 1,nsub(3)

      icell(2) = 1
      xmin_sub(2) = xmin(2)
      DO jsub = 1,nsub(2)

        icell(1) = 1
        xmin_sub(1) = xmin(1)
        DO isub = 1,nsub(1)
!---------------------------------------------------------------------------------!
! Divide the grid points to subdomain (ghost not included)
! If the grid points does not equally divide with the number of processors
! the remaining grid points are distributed equally over the last processors
!---------------------------------------------------------------------------------!
          rem_ncell = ncell - nsub*FLOOR(REAL(ncell,MK)/REAL(nsub,MK))

          IF(isub .GT. nsub(1)-rem_ncell(1))THEN
            ncell_sub(1) = FLOOR(REAL(ncell(1),MK)/REAL(nsub(1),MK)) + 1
          ELSE
            ncell_sub(1) = FLOOR(REAL(ncell(1),MK)/REAL(nsub(1),MK))
          ENDIF

          IF(jsub .GT. nsub(2)-rem_ncell(2))THEN
            ncell_sub(2) = FLOOR(REAL(ncell(2),MK)/REAL(nsub(2),MK)) + 1
          ELSE
            ncell_sub(2) = FLOOR(REAL(ncell(2),MK)/REAL(nsub(2),MK))
          ENDIF

          IF(ksub .GT. nsub(3)-rem_ncell(3))THEN
            ncell_sub(3) = FLOOR(REAL(ncell(3),MK)/REAL(nsub(3),MK)) + 1
          ELSE
            ncell_sub(3) = FLOOR(REAL(ncell(3),MK)/REAL(nsub(3),MK))
          ENDIF

!---------------------------------------------------------------------------------!
! Extend of subdomain (ghost not included)
!---------------------------------------------------------------------------------!  
          xmax_sub = xmin_sub + REAL(ncell_sub,MK)*dx

          IF(isub .EQ. nsub(1))THEN
            xmax_sub(1) = xmax(1)
          ENDIF
          IF(jsub .EQ. nsub(2))THEN
            xmax_sub(2) = xmax(2)
          ENDIF
          IF(ksub .EQ. nsub(3))THEN
            xmax_sub(3) = xmax(3)
          ENDIF

!---------------------------------------------------------------------------------!
! Find boundary conditions, neighbours and number of ghost cells
!---------------------------------------------------------------------------------! 
          IF(isub .EQ. 1)THEN ! exterior west
            IF(bc_ext(1) .EQ. 0)THEN
              bc(1)     = 0
              ineigh(1) = -1
              nghost(1) = nghost_ext
            ELSE
              bc(1)     = 1
              ineigh(1) = sub2proc(nsub(1),jsub,ksub)
              nghost(1) = nghost_int
            END IF
          ELSE ! interior west
            bc(1)     = 1
            ineigh(1) = sub2proc(isub-1,jsub,ksub)
            nghost(1) = nghost_int
          END IF

          IF(isub .EQ. nsub(1))THEN ! exterior east
            IF(bc_ext(1) .EQ. 0)THEN
              bc(2)     = 0
              ineigh(2) = -1
              nghost(2) = nghost_ext
            ELSE
              bc(2)     = 1
              ineigh(2) = sub2proc(1,jsub,ksub)
              nghost(2) = nghost_int
            END IF
          ELSE
            bc(2)     = 1
            ineigh(2) = sub2proc(isub+1,jsub,ksub)
            nghost(2) = nghost_int
          END IF

          IF(jsub .EQ. 1)THEN ! exterior south
            IF(bc_ext(2) .EQ. 0)THEN
              bc(3)     = 0
              ineigh(3) = -1          
              nghost(3) = nghost_ext
            ELSE
              bc(3)     = 1 
              ineigh(3) = sub2proc(isub,nsub(2),ksub)
              nghost(3) = nghost_int
            END IF
          ELSE ! interior south
            bc(3)     = 1
            ineigh(3) = sub2proc(isub,jsub-1,ksub)
            nghost(3) = nghost_int
          END IF

          IF(jsub .EQ. nsub(2))THEN ! exterior north
            IF(bc_ext(2) .EQ. 0)THEN
              bc(4)     = 0 
              ineigh(4) = -1         
              nghost(4) = nghost_ext
            ELSE
              bc(4)     = 1
              ineigh(4) = sub2proc(isub,1,ksub)
              nghost(4) = nghost_int
            END IF
          ELSE ! interior north
            bc(4)     = 1
            ineigh(4) = sub2proc(isub,jsub+1,ksub)
            nghost(4) = nghost_int
          END IF

          IF(ksub .EQ. 1)THEN ! exterior bottom
            IF(bc_ext(3) .EQ. 0)THEN
              bc(5)     = 0
              ineigh(5) = -1      
              nghost(5) = nghost_ext    
            ELSE
              bc(5)     = 1
              ineigh(5) = sub2proc(isub,jsub,nsub(3))
              nghost(5) = nghost_int
            END IF
          ELSE ! interior bottom
            bc(5)     = 1
            ineigh(5) = sub2proc(isub,jsub,ksub-1)
            nghost(5) = nghost_int
          END IF

          IF(ksub .EQ. nsub(3))THEN ! exterior top
            IF(bc_ext(3) .EQ. 0)THEN
              bc(6)     = 0
              ineigh(6) = -1
              nghost(6) = nghost_ext
            ELSE
              bc(6)     = 1
              ineigh(6) = sub2proc(isub,jsub,1)
              nghost(6) = nghost_int
            END IF
          ELSE ! interior top
            bc(6)     = 1
            ineigh(6) = sub2proc(isub,jsub,ksub+1)
            nghost(6) = nghost_int
          END IF
        
!---------------------------------------------------------------------------------!
! Store the topology
!---------------------------------------------------------------------------------!
          iproc = sub2proc(isub,jsub,ksub)

          topo(iproc)%xmin  = xmin_sub
          topo(iproc)%xmax  = xmax_sub
          topo(iproc)%ncell = ncell_sub
          topo(iproc)%icell = icell
          topo(iproc)%dx    = dx

          topo(iproc)%bound_cond = bc
          topo(iproc)%nghost     = nghost
          topo(iproc)%ineigh     = ineigh

          icell(1) = icell(1) + ncell_sub(1)
          xmin_sub(1) = xmax_sub(1)
        END DO
        icell(2) = icell(2) + ncell_sub(2)
        xmin_sub(2) = xmax_sub(2)
      END DO   
      icell(3) = icell(3) + ncell_sub(3)
      xmin_sub(3) = xmax_sub(3)
    END DO

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
    DEALLOCATE(sub2proc,stat=ierr)

  END IF ! mpi_rank0

!---------------------------------------------------------------------------------!
! Broadcast the topology
!---------------------------------------------------------------------------------!
  DO iproc = 0,mpi_nproc-1
    CALL MPI_BCAST(topo(iproc)%xmin, pmlib_ndim, & 
             & mpi_prec_real, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%xmax, pmlib_ndim, & 
             & mpi_prec_real, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%ncell, pmlib_ndim, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%icell, pmlib_ndim, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%dx, pmlib_ndim, & 
             & mpi_prec_real, 0, mpi_comm, ierr)

    CALL MPI_BCAST(topo(iproc)%bound_cond, 2*pmlib_ndim, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%nghost, 2*pmlib_ndim, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
    CALL MPI_BCAST(topo(iproc)%ineigh, 2*pmlib_ndim, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_topology_cuboid
