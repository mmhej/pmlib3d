!---------------------------------------------------------------------------------!
! pmlib_repatch.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_repatch(patch,topo_all,part,mesh,min_vort,buf,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_particles
USE pmlib_mod_mesh
USE pmlib_mod_communication
USE pmlib_mod_poisson

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                         :: patch
  TYPE(class_topology_all)                  :: topo_all
  TYPE(class_particles)                     :: part
  TYPE(class_mesh)                          :: mesh
  REAL(MK)                                  :: min_vort
  INTEGER, INTENT(IN)                       :: buf
  INTEGER, INTENT(INOUT)                    :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=13)                             :: caller = 'pmlib_repatch'
  CHARACTER(LEN=256)                            :: msg
  INTEGER                                       :: i,j
  REAL(MK),DIMENSION(pmlib_ndim)                :: xmin,xmax
  REAL(MK),DIMENSION(pmlib_ndim)                :: dx
  INTEGER                                       :: npart, npart_send
  INTEGER                                       :: npart_old 
  INTEGER                                       :: npart_glo1, npart_glo2 
  REAL(MK)                                      :: mag
  LOGICAL                                       :: map_done, map_done_gl

!---------------------------------------------------------------------------------!
! Initiate subroutine 
!---------------------------------------------------------------------------------!
  ierr = 0
  xmin =  1.0E5_MK
  xmax = -1.0E5_MK

! check number of particles
!  CALL MPI_ALLREDUCE( topo(mpi_rank)%npart,npart_glo1,1, &
!                 & mpi_prec_int,MPI_SUM,mpi_comm,ierr)

!---------------------------------------------------------------------------------!
! Search particles for the min and max
!---------------------------------------------------------------------------------!
  npart_old = topo_all%cuboid(mpi_rank)%npart
  IF (patch%ptype .EQ. 2) THEN

    dx   = patch%dx
    DO i = 1,npart_old

      mag = SQRT( part%vort(1,i)**2 + part%vort(2,i)**2 + &
                & part%vort(3,i)**2 )

      IF( mag .GT. min_vort ) THEN
        DO j = 1,pmlib_ndim
          IF( patch%bound_cond(j) .EQ. 0)THEN
            xmin(j) = MIN(xmin(j),part%pos(j,i))
            xmax(j) = MAX(xmax(j),part%pos(j,i))
          END IF
        END DO
      END IF

    END DO

! Extend with a buffer of n cells
    DO j = 1,pmlib_ndim
      IF( patch%bound_cond(j) .EQ. 0)THEN
        xmin(j) = xmin(j) - (REAL(buf,MK) + 0.5_MK) * dx(j)
        xmax(j) = xmax(j) + (REAL(buf,MK) + 0.5_MK) * dx(j)
      ELSE
        xmin(j) = patch%xmin(j)
        xmax(j) = patch%xmax(j)
      END IF
    END DO

!---------------------------------------------------------------------------------!
! Gather patch extend to find global min/max
!---------------------------------------------------------------------------------!
    CALL MPI_ALLREDUCE( xmin,patch%xmin,pmlib_ndim, &
                      & mpi_prec_real,MPI_MIN,mpi_comm,ierr)
    CALL MPI_ALLREDUCE( xmax,patch%xmax,pmlib_ndim, &
                      & mpi_prec_real,MPI_MAX,mpi_comm,ierr)

!---------------------------------------------------------------------------------!
! Adjust new patches
!---------------------------------------------------------------------------------!
    CALL pmlib_patch_adjust(patch,ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to adjust patches.')
      GOTO 9999
    ENDIF

  END IF

!---------------------------------------------------------------------------------!
! Setup new cuboid topology
!---------------------------------------------------------------------------------!
  CALL pmlib_topology_cuboid(patch,topo_all%cuboid,ierr)
  topo_all%cuboid(mpi_rank)%npart = npart_old
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to create new topology.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Re-allocate mesh
!---------------------------------------------------------------------------------!
  CALL pmlib_mesh_allocate( topo_all%cuboid,mesh,ierr, &
                          & vort = .TRUE., dvort = .TRUE., &
                          & vel = .TRUE., mask = .FALSE.)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to allocate mesh.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Set up Poisson solver
!---------------------------------------------------------------------------------!
  CALL pmlib_poisson_setup(patch,topo_all,mesh,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write( mpi_rank,caller, &
                    & 'Failed to set up Poisson solver.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Set up structure
!---------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------!
! Map particles
!---------------------------------------------------------------------------------!
  DO i = 1,pmlib_ndim
    map_done_gl = .FALSE.
    DO WHILE(.NOT. map_done_gl)
! setup particle communication
      CALL pmlib_particles_map( patch,topo_all%cuboid,part,i, &
                              & pmlib_nvort+pmlib_ndim,ierr)
! pack
      CALL pmlib_comm_pack(part%pos,ierr)
      CALL pmlib_comm_pack(part%vort,ierr)
! send buffers
      CALL pmlib_comm_send(ierr)
! unpack
      CALL pmlib_comm_unpack(topo_all%cuboid,part%vort,ierr)
      CALL pmlib_comm_unpack(topo_all%cuboid,part%pos,map_done,i,ierr)
! finalise
      CALL pmlib_comm_finalise(ierr)
      CALL MPI_ALLREDUCE( map_done,map_done_gl,1,MPI_LOGICAL, &
                        & MPI_LAND,mpi_comm,ierr)
    END DO
  END DO

!---------------------------------------------------------------------------------!
! Check the total number of particles
!---------------------------------------------------------------------------------!
!  CALL MPI_ALLREDUCE(topo_all%cuboid(mpi_rank)%npart, & 
!         & npart_glo2,1,mpi_prec_int,MPI_SUM,mpi_comm,ierr)
!  IF(mpi_rank .EQ. 0)THEN
!    WRITE(*,*)npart_glo1,npart_glo2
!  END IF

!---------------------------------------------------------------------------------!
! Reallocate particle time derivatives
!---------------------------------------------------------------------------------!
  DEALLOCATE(part%vel)
  ALLOCATE(part%vel(pmlib_ndim,topo_all%cuboid(mpi_rank)%npart))

  DEALLOCATE(part%dvort)
  ALLOCATE(part%dvort(pmlib_nvort,topo_all%cuboid(mpi_rank)%npart))

  DEALLOCATE(part%vel_rk)
  ALLOCATE(part%vel_rk(pmlib_ndim,topo_all%cuboid(mpi_rank)%npart))
  part%vel_rk = 0.0_MK

  DEALLOCATE(part%dvort_rk)
  ALLOCATE(part%dvort_rk(pmlib_nvort,topo_all%cuboid(mpi_rank)%npart))
  part%dvort_rk = 0.0_MK

!---------------------------------------------------------------------------------!
! Output new domain size
!---------------------------------------------------------------------------------!
  WRITE(msg,'(A,3I6)')'Re-patching - # of cells ',patch%ncell
  CALL pmlib_write(mpi_rank,caller,TRIM(msg),verb=.TRUE.)

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN


END SUBROUTINE pmlib_repatch

