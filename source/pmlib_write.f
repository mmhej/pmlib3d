!---------------------------------------------------------------------------------!
! pmlib_write.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_write

IMPLICIT NONE

CONTAINS

SUBROUTINE pmlib_write(rank,caller,mesg,verb)
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  INTEGER, INTENT(IN)            :: rank
  CHARACTER(LEN=*), INTENT(IN)   :: caller
  CHARACTER(LEN=*), INTENT(IN)   :: mesg
  LOGICAL,OPTIONAL, INTENT(IN)   :: verb

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER            :: ios
  INTEGER            :: icaller,imsg

!---------------------------------------------------------------------------------!
! Get length of messages
!---------------------------------------------------------------------------------!
  icaller = LEN_TRIM(caller)
  imsg    = LEN_TRIM(mesg)

!---------------------------------------------------------------------------------!
! Write to terminal
!---------------------------------------------------------------------------------!
  IF(PRESENT(verb) .AND. verb)THEN
    IF (rank .EQ. 0) THEN
      WRITE(*,'(4A)',IOSTAT=ios)  '(', caller(1:icaller), ') : ', &
      &      mesg(1:imsg)
    END IF
  ELSE
    WRITE(*,'(A,I4.4,4A)',IOSTAT=ios) '[',rank,'](',caller(1:icaller),') : ', &
    &      mesg(1:imsg)
  END IF

!---------------------------------------------------------------------------------!
!  Return
!---------------------------------------------------------------------------------!
RETURN
END SUBROUTINE pmlib_write
END MODULE pmlib_mod_write
