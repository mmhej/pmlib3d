!---------------------------------------------------------------------------------!
! pmlib_topology_kernel.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_topology_kernel(patch,topo,pencil_dir,ierr,fextend,nextend)

USE pmlib_mod_patch

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch),INTENT(IN)                              :: patch
  TYPE(class_topology),DIMENSION(:),POINTER,INTENT(OUT)     :: topo
  INTEGER,INTENT(IN)                                        :: pencil_dir
  INTEGER, INTENT(OUT)                                      :: ierr
  INTEGER,DIMENSION(pmlib_ndim),OPTIONAL,INTENT(IN)         :: fextend
  INTEGER,DIMENSION(pmlib_ndim),OPTIONAL,INTENT(IN)         :: nextend

!---------------------------------------------------------------------------------!
!  Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=21)                           :: caller = 'pmlib_topology_kernel'
  INTEGER                                     :: iproc
  INTEGER                                     :: isub,jsub,ksub, lsub
  INTEGER                                     :: i,j,k
  INTEGER                                     :: nsub_dif_min
  INTEGER, DIMENSION(pmlib_ndim)              :: nsub

  INTEGER, DIMENSION(pmlib_ndim)              :: ncell, icell
  REAL(MK), DIMENSION(pmlib_ndim)             :: xmin, xmax, dx 
  INTEGER, DIMENSION(pmlib_ndim)              :: bc_ext
  INTEGER, DIMENSION(2*pmlib_ndim)            :: bc, nghost
  INTEGER, DIMENSION(2*pmlib_ndim)            :: ineigh
  INTEGER                                     :: nghost_ext, nghost_int

  INTEGER, DIMENSION(pmlib_ndim)              :: rem_ncell

  REAL(MK), DIMENSION(:,:,:), POINTER         :: sub2proc
  REAL(MK), DIMENSION(:,:,:), POINTER         :: proc_sub

  INTEGER, DIMENSION(pmlib_ndim)              :: ncell_sub
  REAL(MK), DIMENSION(pmlib_ndim)             :: xmin_sub, xmax_sub

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

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
  bc_ext = 0

!---------------------------------------------------------------------------------!
! Determine topology parameters
!---------------------------------------------------------------------------------!
  dx   = patch%dx
  xmin = patch%xmin
  xmax = patch%xmax

!---------------------------------------------------------------------------------!
! Include ghost cells in the domain if free boundary conditions
! Correct extend of topology and number of mesh and ghost cells
!---------------------------------------------------------------------------------!
  nghost_ext = 0
  nghost_int = 0

  DO i = 1,pmlib_ndim
    IF( patch%bound_cond(i) .EQ. 0 )THEN
!---------------------------------------------------------------------------------!
! Include ghost cells
!---------------------------------------------------------------------------------!
      ncell(i) = patch%ncell(i) + 2*patch%nghost
!---------------------------------------------------------------------------------!
! Extend domain in the specified dimensions
!---------------------------------------------------------------------------------!
      IF(PRESENT(fextend) .AND. fextend(i) .GT. 1)THEN
        ncell(i) = ncell(i) * fextend(i)
      END IF
      IF(PRESENT(nextend) .AND. nextend(i) .GT. 0)THEN
        ncell(i) = ncell(i) + nextend(i)
      END IF

    ELSE

!---------------------------------------------------------------------------------!
! Exclude ghost cells in the domain if periodic boundary conditions
!---------------------------------------------------------------------------------!
      ncell(i)  = patch%ncell(i)
    END IF

!---------------------------------------------------------------------------------!
! Set domain maximum and minimum
!---------------------------------------------------------------------------------!
    xmin(i) = - ( 0.5_MK * REAL(ncell(i),MK) + 0.5_MK)*dx(i)
    xmax(i) =   ( 0.5_MK * REAL(ncell(i),MK) - 0.5_MK)*dx(i)
  END DO

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
    nsub_dif_min = mpi_nproc

    DO i = mpi_nproc,1,-1
      DO j = 1,mpi_nproc

        IF (pencil_dir .EQ. 1)THEN
          IF(i*j .EQ. mpi_nproc .AND. ABS(i-j) .LE. nsub_dif_min)THEN
            nsub = (/ 1 , i , j /) 
            nsub_dif_min = ABS(i-j)
          END IF
        ELSEIF (pencil_dir .EQ. 2)THEN
          IF(i*j .EQ. mpi_nproc .AND. ABS(i-j) .LE. nsub_dif_min)THEN
            nsub = (/ i , 1 , j /) 
            nsub_dif_min = ABS(i-j)
          END IF
        ELSEIF (pencil_dir .EQ. 3)THEN
          IF(i*j .EQ. mpi_nproc .AND. ABS(i-j) .LE. nsub_dif_min)THEN
            nsub = (/ i , j , 1 /) 
            nsub_dif_min = ABS(i-j)
          END IF
        END IF

      END DO
    END DO

!---------------------------------------------------------------------------------!
! Check the total number of subdomains
!---------------------------------------------------------------------------------!
    IF(nsub(1)*nsub(2)*nsub(3) .NE. mpi_nproc)THEN
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Wrong number of subdomains.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Assign subdomains and keep track of the local index of the subdomain on
! the processor
!---------------------------------------------------------------------------------!
    ALLOCATE(sub2proc( nsub(1),nsub(2),nsub(3)) )
    ALLOCATE(proc_sub( nsub(1),nsub(2),nsub(3)) )
    sub2proc = 0
    proc_sub = 0
    iproc = 0
    lsub  = 1
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
    DEALLOCATE(proc_sub,stat=ierr)

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

END SUBROUTINE pmlib_topology_kernel

