!---------------------------------------------------------------------------------!
! pmlib_comm_send.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
! Sends the communication buffer
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_comm_send(ierr)

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  INTEGER, INTENT(OUT)                     :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=15)                        :: caller = 'pmlib_comm_send'
  INTEGER                                  :: isend, icomm, communicator
  INTEGER                                  :: iproc, jproc
  INTEGER                                  :: stat(MPI_STATUS_SIZE)
  INTEGER                                  :: i,j,k,l

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Send/recieve to/from other processors
!---------------------------------------------------------------------------------!
  DO isend = 1,communication%ncommseq ! for number of communications

    communicator = communication%icommseq(1,mpi_rank,isend)
    icomm        = communication%icommseq(2,mpi_rank,isend)

!---------------------------------------------------------------------------------!
! Senders
!---------------------------------------------------------------------------------!
    IF(communicator .GT. 0)THEN

      iproc = mpi_rank
      jproc = communicator - 1

      IF( communication%ctype .EQ. 1)THEN

        CALL MPI_SEND( comm_buffer(icomm)%mesh_send, &
                     & comm_buffer(icomm)%nsend,mpi_prec_real, &
                     & jproc,jproc,mpi_comm,stat,ierr )

      ELSE IF(communication%ctype .EQ. 2)THEN

        CALL MPI_SEND( comm_buffer(icomm)%part_send, &
                     & comm_buffer(icomm)%nsend,mpi_prec_real, &
                     & jproc,jproc,mpi_comm,stat,ierr )

      ELSE

        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Communication type unknown.')
        GOTO 9999

      END IF
!---------------------------------------------------------------------------------!
! Recievers
!---------------------------------------------------------------------------------!
    ELSEIF(communicator .LT. 0)THEN

      iproc = -(communicator + 1)  
      jproc = mpi_rank

      IF( communication%ctype .EQ. 1)THEN

        CALL MPI_RECV( comm_buffer(icomm)%mesh_recieve, &
                     & comm_buffer(icomm)%nrecieve,mpi_prec_real, &
                     & iproc,jproc,mpi_comm,stat,ierr )

      ELSEIF(communication%ctype .EQ. 2)THEN

        CALL MPI_RECV( comm_buffer(icomm)%part_recieve, &
                     & comm_buffer(icomm)%nrecieve, &
                     & mpi_prec_real,iproc,jproc,mpi_comm,stat,ierr )

      ELSE

        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Communication type unknown.')
        GOTO 9999

      END IF
    END IF

  ENDDO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_send

