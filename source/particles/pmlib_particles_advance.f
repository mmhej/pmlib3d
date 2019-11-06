!---------------------------------------------------------------------------------!
! pmlib_particles_advance.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_particles_advance( patch,topo,part,itimestage,dtime,ierr, &
                                  & euler_vort)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_communication

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch), INTENT(IN)                            :: patch
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(INOUT) :: topo
  TYPE(class_particles), INTENT(INOUT)                     :: part
  INTEGER, INTENT(IN)                                      :: itimestage
  REAL(MK), INTENT(IN)                                     :: dtime
  INTEGER, INTENT(OUT)                                     :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                            :: euler_vort

!---------------------------------------------------------------------------------!
!  Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=23)                   :: caller = 'pmlib_particles_advance'
  INTEGER                                   :: i,j,k

  REAL(MK),DIMENSION(pmlib_ndim)            :: xmin, xmax
  INTEGER,DIMENSION(pmlib_ndim)             :: bc

  REAL(MK),DIMENSION(3)                     :: a1, a2, b1, b2, c1, c2

  REAL(MK)                                  :: vort_max_loc

  LOGICAL                                   :: map_done, map_done_gl
  INTEGER                                   :: npart_glo1, npart_glo2 

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Set time integration scheme coefficients
!---------------------------------------------------------------------------------!
  IF(pmlib_time_integration_order .EQ. 1)THEN
! Position
    a1 = (/ 0.0_MK, 0.0_MK, 0.0_MK /)
    b1 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
    c1 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
! Vorticity
    a2 = (/ 0.0_MK, 0.0_MK, 0.0_MK /)
    b2 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
    c2 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
  ELSEIF(pmlib_time_integration_order .EQ. 2)THEN
! Position
    a1 = (/ 0.0_MK, -0.5_MK, 0.0_MK /)
    b1 = (/ 1.0_MK,  1.0_MK, 0.0_MK /)
    c1 = (/ 0.5_MK,  1.0_MK, 0.0_MK /)
! Vorticity
    IF( PRESENT(euler_vort) .AND. euler_vort )THEN
      a2 = (/ 0.0_MK, 1.0_MK, 0.0_MK /)
      b2 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
      c2 = (/ 0.5_MK, 0.5_MK, 0.0_MK /)
    ELSE 
      a2 = (/ 0.0_MK, -0.5_MK, 0.0_MK /)
      b2 = (/ 1.0_MK,  1.0_MK, 0.0_MK /)
      c2 = (/ 0.5_MK,  1.0_MK, 0.0_MK /)
    END IF
  ELSEIF(pmlib_time_integration_order .EQ. 3)THEN
! Position
    a1 = (/ 0.0_MK, -5.0_MK/9.0_MK, -153.0_MK/128.0_MK /)
    b1 = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
    c1 = (/ 1.0_MK/3.0_MK, 15.0_MK/16.0_MK, 8.0_MK/15.0_MK /)
! Vorticity
    IF( PRESENT(euler_vort) .AND. euler_vort )THEN
      a2 = (/ 0.0_MK, 1.0_MK, 1.0_MK /)
      b2 = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
      c2 = (/ 1.0_MK/3.0_MK, 1.0_MK/3.0_MK, 1.0_MK/3.0_MK /)
    ELSE 
      a2 = (/ 0.0_MK, -5.0_MK/9.0_MK, -153.0_MK/128.0_MK /)
      b2 = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
      c2 = (/ 1.0_MK/3.0_MK, 15.0_MK/16.0_MK, 8.0_MK/15.0_MK /)
    END IF
  END IF

!---------------------------------------------------------------------------------!
! Advance particle possition and vorticity
!---------------------------------------------------------------------------------!
  vort_max_loc = 0.0_MK

  DO i = 1,topo(mpi_rank)%npart
! the storage 
    part%vel_rk(1:pmlib_ndim,i) = &
                  &   a1(itimestage)*part%vel_rk(1:pmlib_ndim,i) &
                  & + b1(itimestage)*dtime*part%vel(1:pmlib_ndim,i)

    part%dvort_rk(1:pmlib_nvort,i) = &
                  &   a2(itimestage)*part%dvort_rk(1:pmlib_nvort,i) &
                  & + b2(itimestage)*dtime*part%dvort(1:pmlib_nvort,i)


! the time integration
    part%pos(1:pmlib_ndim,i) = part%pos(1:pmlib_ndim,i) &
                  & + c1(itimestage)*part%vel_rk(1:pmlib_ndim,i)

    part%vort(1:pmlib_nvort,i) = part%vort(1:pmlib_nvort,i) &
                  & + c2(itimestage)*part%dvort_rk(1:pmlib_nvort,i)

  END DO !ipart

!  CALL MPI_ALLREDUCE( topo(mpi_rank)%npart,npart_glo1,1,mpi_prec_int, &
!                    & MPI_SUM,comm,ierr)

!---------------------------------------------------------------------------------!
! Communicate particles
!---------------------------------------------------------------------------------!
  IF(pmlib_time_integration_order .EQ. 1)THEN
    DO i = 1,pmlib_ndim
      map_done_gl = .FALSE.
      DO WHILE(.NOT. map_done_gl)
! setup particle communication
        CALL pmlib_particles_map( patch,topo,part,i, &
                                & pmlib_nvort+pmlib_ndim,ierr)
! pack
        CALL pmlib_comm_pack(part%pos,ierr)
        CALL pmlib_comm_pack(part%vort,ierr)
! send buffers
        CALL pmlib_comm_send(ierr)
! unpack
        CALL pmlib_comm_unpack(topo,part%vort,ierr)
        CALL pmlib_comm_unpack(topo,part%pos,map_done,i,ierr)
! finalise
        CALL pmlib_comm_finalise(ierr)
        CALL MPI_ALLREDUCE( map_done,map_done_gl,1,MPI_LOGICAL, &
                          & MPI_LAND,mpi_comm,ierr)
      END DO
    END DO

  ELSE
    DO i = 1,pmlib_ndim
      map_done_gl = .FALSE.
      DO WHILE(.NOT. map_done_gl)
! setup particle communication
        CALL pmlib_particles_map( patch,topo,part,i, &
                                & 2*(pmlib_nvort+pmlib_ndim),ierr)
! pack
        CALL pmlib_comm_pack(part%pos,ierr)
        CALL pmlib_comm_pack(part%vort,ierr)
        CALL pmlib_comm_pack(part%vel_rk,ierr)
        CALL pmlib_comm_pack(part%dvort_rk,ierr)
! send buffers
        CALL pmlib_comm_send(ierr)
! unpack
        CALL pmlib_comm_unpack(topo,part%dvort_rk,ierr)
        CALL pmlib_comm_unpack(topo,part%vel_rk,ierr)
        CALL pmlib_comm_unpack(topo,part%vort,ierr)
        CALL pmlib_comm_unpack(topo,part%pos,map_done,i,ierr)
! finalise
        CALL pmlib_comm_finalise(ierr)
        CALL MPI_ALLREDUCE( map_done,map_done_gl,1,MPI_LOGICAL, &
                          & MPI_LAND,mpi_comm,ierr)
      END DO
    END DO

  END IF

!---------------------------------------------------------------------------------!
! Check the total number of particles
!---------------------------------------------------------------------------------!
!  IF(SIZE(part%pos,2) .NE. topo(mpi_rank)%npart)THEN
!    CALL pmlib_write( mpi_rank,caller, &
!                    & 'Particle array does not correspond.')
!  END IF

!  CALL MPI_ALLREDUCE( topo(mpi_rank)%npart,npart_glo2,1,mpi_prec_int, &
!                    & MPI_SUM,comm,ierr)

!  IF( npart_glo1 .NE. npart_glo2 )THEN
!    ierr = -1
!    CALL pmlib_write( mpi_rank,caller, &
!                    & 'Number of particles are not conserved.',verb=.TRUE.)
!    GOTO 9999
!  END IF

!---------------------------------------------------------------------------------!
! Reallocate particle time derivatives
!---------------------------------------------------------------------------------!
  DEALLOCATE(part%vel)
  ALLOCATE(part%vel(pmlib_ndim,topo(mpi_rank)%npart))

  DEALLOCATE(part%dvort)
  ALLOCATE(part%dvort(pmlib_ndim,topo(mpi_rank)%npart))

  IF(pmlib_time_integration_order .EQ. 1)THEN
    DEALLOCATE(part%vel_rk)
    ALLOCATE(part%vel_rk(pmlib_ndim,topo(mpi_rank)%npart))
    part%vel_rk = 0.0_MK

    DEALLOCATE(part%dvort_rk)
    ALLOCATE(part%dvort_rk(pmlib_ndim,topo(mpi_rank)%npart))
    part%dvort_rk = 0.0_MK
  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_particles_advance

