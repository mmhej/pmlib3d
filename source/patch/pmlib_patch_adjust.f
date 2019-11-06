!---------------------------------------------------------------------------------!
! pmlib_patch_adjust.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_patch_adjust(patch,ierr)

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                  :: patch
  INTEGER, INTENT(INOUT)             :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=18)                 :: caller = 'pmlib_patch_adjust' 
  INTEGER                           :: i

  INTEGER,DIMENSION(pmlib_ndim)     :: extend

  REAL(MK)                          :: centre
  INTEGER                           :: ncell 

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Adjust patch on mpi_rank 0 and communicate
!---------------------------------------------------------------------------------!
  IF(mpi_rank .EQ. 0)THEN 
    DO i = 1,pmlib_ndim

      IF( patch%bound_cond(i) .EQ. 0 )THEN

!---------------------------------------------------------------------------------!
! Calculate the patch centre
!---------------------------------------------------------------------------------!
        centre = 0.5_MK * (patch%xmax(i) + patch%xmin(i))

!---------------------------------------------------------------------------------!
! Calculate the minimum number of interiour mesh cells
!---------------------------------------------------------------------------------!
        ncell = CEILING( (patch%xmax(i) - patch%xmin(i))/patch%dx(i) )

!---------------------------------------------------------------------------------!
! Extend the number of cells to a number such that the total number of 
! mesh cells is a multiple of 16 (to make the fft more efficient) 
!---------------------------------------------------------------------------------!
        ncell = ncell + (16 - MOD(ncell + 2*patch%nghost,16) )

!---------------------------------------------------------------------------------!
! Adjust extend of patch to the number of cells
!---------------------------------------------------------------------------------!
        patch%xmin(i) = centre &
                & - REAL( NINT( 0.5_MK * REAL(ncell,MK)),MK ) * patch%dx(i)
        patch%xmax(i) = centre &
                & + REAL( NINT( 0.5_MK * REAL(ncell,MK)),MK ) * patch%dx(i)

!---------------------------------------------------------------------------------!
! Recalculate the number of mesh cells
!---------------------------------------------------------------------------------!
        patch%ncell(i) = NINT( (patch%xmax(i) - patch%xmin(i))/patch%dx(i) )
        IF( patch%ncell(i) .NE. ncell )THEN
          ierr = - 1
          CALL pmlib_write( mpi_rank,caller, &
                          & 'Failed to adjust the number of mesh cells')
          GO TO 9999
        END IF

      ELSE
!---------------------------------------------------------------------------------!
! Extend of patch (fit to integer number of cells)
!---------------------------------------------------------------------------------!
        patch%xmax(i) = patch%xmin(i) + patch%dx(i) * &
             & REAL(NINT( (patch%xmax(i) - patch%xmin(i))/patch%dx(i) ),MK )

!---------------------------------------------------------------------------------!
! Calculate the number of interiour mesh cells
!---------------------------------------------------------------------------------!
        patch%ncell(i) = NINT( (patch%xmax(i) - patch%xmin(i))/patch%dx(i) )

      END IF
    END DO

!---------------------------------------------------------------------------------!
! Find fft-sequence
!---------------------------------------------------------------------------------!
    DO i = 1,pmlib_ndim
      IF( patch%bound_cond(i) .EQ. 0 )THEN
        extend(i) = 2*patch%ncell(i)
      ELSE
        extend(i) = patch%ncell(i)
      END IF
    END DO

! For memory saving N1 < N2 < N3

! For speed N1 > N2 > N3 
    IF( extend(1) .GE. extend(2) .AND. &
      & extend(1) .GE. extend(3) )THEN

      IF( extend(2) .GE. extend(3) )THEN
        patch%fft_seq(1) = 1 
        patch%fft_seq(2) = 2
        patch%fft_seq(3) = 3
      ELSE
        patch%fft_seq(1) = 1 
        patch%fft_seq(2) = 3
        patch%fft_seq(3) = 2
      END IF

    ELSEIF( extend(2) .GE. extend(1) .AND. &
          & extend(2) .GE. extend(3) )THEN

      IF( extend(1) .GE. extend(3) )THEN
        patch%fft_seq(1) = 2 
        patch%fft_seq(2) = 1
        patch%fft_seq(3) = 3
      ELSE
        patch%fft_seq(1) = 2 
        patch%fft_seq(2) = 3
        patch%fft_seq(3) = 1
      END IF
    ELSEIF( extend(3) .GE. extend(1) .AND. &
          & extend(3) .GE. extend(2) )THEN

      IF( extend(1) .GE. extend(2) )THEN
        patch%fft_seq(1) = 3 
        patch%fft_seq(2) = 1
        patch%fft_seq(3) = 2
      ELSE
        patch%fft_seq(1) = 3 
        patch%fft_seq(2) = 2
        patch%fft_seq(3) = 1
      END IF
    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Cant figure out fft sequence.')
      GOTO 9999        
    END IF

  END IF ! mpi_rank0

!---------------------------------------------------------------------------------!
! Broadcast the patch layout
!---------------------------------------------------------------------------------!
  CALL MPI_BCAST( patch%xmin, pmlib_ndim, mpi_prec_real, 0, mpi_comm, ierr)
  CALL MPI_BCAST( patch%xmax, pmlib_ndim, mpi_prec_real, 0, mpi_comm, ierr)
  CALL MPI_BCAST(   patch%dx, pmlib_ndim, mpi_prec_real, 0, mpi_comm, ierr)

  CALL MPI_BCAST(      patch%ncell,pmlib_ndim,mpi_prec_int,0,mpi_comm,ierr)
  CALL MPI_BCAST(     patch%nghost,         1,mpi_prec_int,0,mpi_comm,ierr)
  CALL MPI_BCAST( patch%bound_cond,pmlib_ndim,mpi_prec_int,0,mpi_comm,ierr)
  CALL MPI_BCAST(    patch%fft_seq,pmlib_ndim,mpi_prec_int,0,mpi_comm,ierr)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_patch_adjust
