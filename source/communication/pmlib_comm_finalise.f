!---------------------------------------------------------------------------------!
! pmlib_comm_finalise.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
! Sets up the inter topology communication
!   Defines the communication order/sequence
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_comm_finalise(ierr)

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  INTEGER, INTENT(OUT)                        :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                           :: caller = 'pmlib_comm_finalise'
  INTEGER                                     :: i

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Deallocate buffer
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(comm_buffer) ) THEN
    DO i = 1,communication%ncomm
      IF( ASSOCIATED(comm_buffer(i)%mesh_send) ) THEN
        DEALLOCATE(comm_buffer(i)%mesh_send)
      END IF  
      IF( ASSOCIATED(comm_buffer(i)%mesh_recieve) ) THEN
        DEALLOCATE(comm_buffer(i)%mesh_recieve)
      END IF  
      IF( ASSOCIATED(comm_buffer(i)%part_send) ) THEN
        DEALLOCATE(comm_buffer(i)%part_send)
      END IF  
      IF( ASSOCIATED(comm_buffer(i)%part_recieve) ) THEN
        DEALLOCATE(comm_buffer(i)%part_recieve)
      END IF  
    END DO

    DEALLOCATE(comm_buffer)
  END IF  

!---------------------------------------------------------------------------------!
! Deallocate communication parameters
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(communication%proc_send) ) THEN
    DEALLOCATE(communication%proc_send)
  END IF  
  IF( ASSOCIATED(communication%proc_recieve) ) THEN
    DEALLOCATE(communication%proc_recieve)
  END IF 

  IF( ASSOCIATED(communication%imin_send) ) THEN
    DEALLOCATE(communication%imin_send)
  END IF 

  IF( ASSOCIATED(communication%imax_send) ) THEN
    DEALLOCATE(communication%imax_send)
  END IF 

  IF( ASSOCIATED(communication%imin_recieve) ) THEN
    DEALLOCATE(communication%imin_recieve)
  END IF 

  IF( ASSOCIATED(communication%imax_recieve) ) THEN
    DEALLOCATE(communication%imax_recieve)
  END IF 

  IF( ASSOCIATED(communication%ncell_send) ) THEN
    DEALLOCATE(communication%ncell_send)
  END IF

  IF( ASSOCIATED(communication%ncell_recieve) ) THEN
    DEALLOCATE(communication%ncell_recieve)
  END IF

  IF( ASSOCIATED(communication%npart_send) ) THEN
    DEALLOCATE(communication%npart_send)
  END IF

  IF( ASSOCIATED(communication%npart_recieve) ) THEN
    DEALLOCATE(communication%npart_recieve)
  END IF

  IF( ASSOCIATED(communication%ipart0) ) THEN
    DEALLOCATE(communication%ipart0)
  END IF

  IF( ASSOCIATED(communication%ipart1) ) THEN
    DEALLOCATE(communication%ipart1)
  END IF

  IF( ASSOCIATED(communication%ipart2) ) THEN
    DEALLOCATE(communication%ipart2)
  END IF

  IF( ASSOCIATED(communication%comm_dir) ) THEN
    DEALLOCATE(communication%comm_dir)
  END IF

  IF( ASSOCIATED(communication%icommseq) ) THEN
    DEALLOCATE(communication%icommseq)
  END IF

  communication%npacked  = 0
  communication%ncomm    = 0
  communication%ncommseq = 0

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_finalise

