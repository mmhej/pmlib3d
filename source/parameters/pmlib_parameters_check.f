!---------------------------------------------------------------------------------!
! pmlib_parameters_check.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_parameters_check(ierr)

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  INTEGER, INTENT(OUT) :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=22)                 :: caller = 'pmlib_parameters_check'

!---------------------------------------------------------------------------------!
! Initialise subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! MPI parameters
!---------------------------------------------------------------------------------!
  IF( mpi_rank .EQ. -99 )THEN
    ierr = -1
    CALL pmlib_write(0,caller,'mpi_rank is incorrectly set.')
    GO TO 9999
  END IF
  IF( mpi_nproc .EQ. -99 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'mpi_nproc is incorrectly set.')
    GO TO 9999
  END IF
  IF( mpi_comm .EQ. -99 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'mpi_comm is incorrectly set.')
    GO TO 9999
  END IF

!---------------------------------------------------------------------------------!
! Numerical parameters
!---------------------------------------------------------------------------------!
  IF( pmlib_poisson_kernel .NE. 1 .AND. & 
    & pmlib_poisson_kernel .NE. 2 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, & 
        & 'pmlib_poisson_kernel is incorrectly set.')
    GO TO 9999
  END IF

  IF( pmlib_poisson_order .NE. 0   .AND. &
    & pmlib_poisson_order .NE. 2   .AND. &
    & pmlib_poisson_order .NE. 4   .AND. &
    & pmlib_poisson_order .NE. 6   .AND. &
    & pmlib_poisson_order .NE. 8   .AND. &
    & pmlib_poisson_order .NE. 10  .AND. &
    & pmlib_poisson_order .NE. 100 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
        & 'pmlib_poisson_order is incorrectly set.')
    GO TO 9999
  END IF

  IF( pmlib_poisson_order .GT. 0 .AND. &
    & pmlib_regularisation_radius .LT. 1.0_MK )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
        & 'pmlib_regularisation_radius is incorrectly set.')
    GO TO 9999
  END IF

  IF( pmlib_fd_order .NE. 2 .AND. &
    & pmlib_fd_order .NE. 4 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
        & 'pmlib_fd_order is incorrectly set.')
    GO TO 9999
  END IF

  IF( pmlib_interpolation_order .NE. 3 .AND. &
    & pmlib_interpolation_order .NE. 4 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
        & 'pmlib_interpolation_order is incorrectly set.')
    GO TO 9999
  END IF

  IF( pmlib_time_integration_order .NE. 1 .AND. &
    & pmlib_time_integration_order .NE. 2 .AND. &
    & pmlib_time_integration_order .NE. 3 )THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
        & 'pmlib_interpolation_order is incorrectly set.')
    GO TO 9999
  END IF


!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_parameters_check
