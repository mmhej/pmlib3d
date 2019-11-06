!---------------------------------------------------------------------------------!
! pmlib_particles_map.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_particles_map(patch,topo,part,dir,nvar_send,ierr)

USE pmlib_mod_communication
USE pmlib_mod_patch
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                         :: patch
  TYPE(class_topology),DIMENSION(:),POINTER :: topo
  TYPE(class_particles)                     :: part
  INTEGER, INTENT(IN)                       :: dir !direction
  INTEGER, INTENT(IN)                       :: nvar_send
  INTEGER, INTENT(OUT)                      :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=19)                     :: caller = 'pmlib_particles_map'
  INTEGER                               :: ncomm, icomm
  INTEGER                               :: iproc, jproc
  INTEGER                               :: i,j,k

  INTEGER, DIMENSION(0:2,0:mpi_nproc-1) :: npart_comm 
  INTEGER, DIMENSION(0:2)               :: npart_dir

  INTEGER                               :: ineigh, neigh
  INTEGER                               :: ipart0, ipart1, ipart2
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin, xmax
  REAL(MK), DIMENSION(pmlib_ndim)       :: dx

  INTEGER                               :: npart_send, npart_recieve
  INTEGER                               :: npart_remove

  INTEGER, DIMENSION(0:mpi_nproc-1)     :: proc_busy 
  INTEGER, DIMENSION(:), POINTER        :: proc_comm
  INTEGER                               :: ncommseq
  INTEGER, DIMENSION(:,:,:), POINTER    :: icommseq 

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find number of particles to communicate
!   index 0 is not moved
!   index 1 is moved in the negative direction 
!   index 2 is moved in the possitive direction
!---------------------------------------------------------------------------------!
  xmin = topo(mpi_rank)%xmin
  xmax = topo(mpi_rank)%xmax 
  npart_comm = 0
  npart_dir  = 0
  npart_remove = 0

  DO i = 1,topo(mpi_rank)%npart
    
    IF( part%pos(dir,i) .LT. xmin(dir) )THEN

      IF( topo(mpi_rank)%bound_cond(2*dir-1) .NE. 0 )THEN
        npart_dir(1) = npart_dir(1) + 1
      ELSE
        npart_remove = npart_remove + 1
      END IF

    ELSEIF( part%pos(dir,i) .GE. xmax(dir) )THEN 

      IF( topo(mpi_rank)%bound_cond(2*dir) .NE. 0 )THEN
        npart_dir(2) = npart_dir(2) + 1
      ELSE
        npart_remove = npart_remove + 1
      END IF

    ELSE
      npart_dir(0) = npart_dir(0) + 1
    END IF

  END DO

!  IF(npart_dir(1) .NE. 0)THEN
!    CALL pmlib_write(mpi_rank,caller,'Communicate particles')
!  ELSEIF(npart_dir(2) .NE. 0)THEN
!    CALL pmlib_write(mpi_rank,caller,'Communicate particles')
!  END IF

!  IF(npart_remove .NE. 0)THEN
!    CALL pmlib_write(mpi_rank,caller,'Particles out of bounce! - removing')
!  END IF

!---------------------------------------------------------------------------------!
! Find indexes of particles to communicate
!---------------------------------------------------------------------------------!
  ALLOCATE( communication%ipart0(npart_dir(0)) )
  ALLOCATE( communication%ipart1(npart_dir(1)) )
  ALLOCATE( communication%ipart2(npart_dir(2)) )

  ipart0 = 0
  ipart1 = 0
  ipart2 = 0
  DO i = 1,topo(mpi_rank)%npart
    
    IF( part%pos(dir,i) .LT. xmin(dir) )THEN
      IF( topo(mpi_rank)%bound_cond(2*dir-1) .NE. 0 )THEN
        ipart1 = ipart1 + 1
        communication%ipart1(ipart1) = i
        IF( part%pos(dir,i) .LT. patch%xmin(dir) )THEN
          part%pos(dir,i) = part%pos(dir,i) &
                        & + ( patch%xmax(dir) - patch%xmin(dir) )
        END IF
      END IF

    ELSEIF( part%pos(dir,i) .GE. xmax(dir) )THEN
      IF( topo(mpi_rank)%bound_cond(2*dir) .NE. 0 )THEN
        ipart2 = ipart2 + 1
        communication%ipart2(ipart2) = i
        IF( part%pos(dir,i) .GE. patch%xmax(dir) )THEN
          part%pos(dir,i) = part%pos(dir,i) &
                        & - ( patch%xmax(dir) - patch%xmin(dir) )
        END IF
      END IF

    ELSE
      ipart0 = ipart0 + 1
      communication%ipart0(ipart0) = i
    END IF
  END DO

  IF( ipart0 .NE. npart_dir(0) .OR. &
    & ipart1 .NE. npart_dir(1) .OR. &
    & ipart2 .NE. npart_dir(2) )THEN
    ierr = -1
    CALL pmlib_write( mpi_rank,caller, &
              & 'Number of communicated particles are not consistent' )
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Gather
!---------------------------------------------------------------------------------!
  CALL MPI_GATHER( npart_dir,3,mpi_prec_int, &
                 & npart_comm,3,mpi_prec_int,0,mpi_comm,ierr )

!---------------------------------------------------------------------------------!
! Setup communication parameters on mpi_rank 0 and distribute
!---------------------------------------------------------------------------------!
  IF (mpi_rank .EQ. 0) THEN
!---------------------------------------------------------------------------------!
! Find number of communications
!---------------------------------------------------------------------------------!
  icomm = 0
  DO i = 1,2
    ineigh = 2*dir - 2 + i
    DO iproc = 0,mpi_nproc-1

      neigh = topo(iproc)%ineigh(ineigh)

      IF( npart_comm(i,iproc) .NE. 0 .AND. neigh .NE. -1)THEN
        icomm = icomm + 1
      END IF

    END DO !iproc
  END DO !ineigh

! add self communications
  icomm = icomm + mpi_nproc

!---------------------------------------------------------------------------------!
! Allocate communication lists
!---------------------------------------------------------------------------------!
  ncomm               = icomm
  communication%ncomm = ncomm
  communication%ctype = 2

  ALLOCATE( communication%proc_send(ncomm) )
  ALLOCATE( communication%proc_recieve(ncomm) )

  ALLOCATE( communication%npart_send(ncomm) )
  ALLOCATE( communication%npart_recieve(ncomm) )
 
  ALLOCATE( communication%comm_dir(ncomm) )

  ALLOCATE( proc_comm(ncomm) )
  proc_comm = 0

!---------------------------------------------------------------------------------!
! Loop through neighbours
!---------------------------------------------------------------------------------!
  icomm = 0
  DO i = 1,2
    ineigh = 2*dir - 2 + i
    DO iproc = 0,mpi_nproc-1

      neigh = topo(iproc)%ineigh(ineigh)

      IF( npart_comm(i,iproc) .NE. 0  .AND. neigh .NE. -1 )THEN

!---------------------------------------------------------------------------------!
! Count communications, tag the inter-processor communications
!---------------------------------------------------------------------------------!
        icomm = icomm + 1
        IF(iproc .NE. neigh)THEN
          proc_comm(icomm) = 1
        END IF
        communication%comm_dir(icomm) = i

!---------------------------------------------------------------------------------!
! Store communicators
!---------------------------------------------------------------------------------!
        communication%proc_send(icomm)     = iproc
        communication%proc_recieve(icomm)  = neigh

        communication%npart_send(icomm)    = npart_comm(i,iproc)
        communication%npart_recieve(icomm) = npart_comm(i,iproc)

      END IF ! send part

    END DO !iproc
  END DO !ineigh

!---------------------------------------------------------------------------------!
! Do self communications
!---------------------------------------------------------------------------------!
  DO iproc = 0,mpi_nproc-1
    icomm = icomm + 1
    proc_comm(icomm) = 0

    communication%comm_dir(icomm) = 0

! Store communicators etc.
    communication%proc_send(icomm)     = iproc
    communication%proc_recieve(icomm)  = iproc

    communication%npart_send(icomm)    = npart_comm(0,iproc)
    communication%npart_recieve(icomm) = npart_comm(0,iproc)
  END DO

  IF(icomm .NE. ncomm)THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
                     & 'Number of communications are inconsistent')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Find communication sequence
!---------------------------------------------------------------------------------!
  ncommseq = 0
  ALLOCATE(icommseq(1:2,0:mpi_nproc-1,ncomm))
  icommseq = 0

  DO WHILE(MAXVAL(proc_comm) .EQ. 1)
    ncommseq = ncommseq + 1 ! next sequence step
    proc_busy = 0           ! free processors

    DO j = 1,ncomm

      iproc = communication%proc_send(j)
      jproc = communication%proc_recieve(j)

!---------------------------------------------------------------------------------!
! Send to/Recieve from (indicated by +/- subdomain)
!---------------------------------------------------------------------------------!
      IF( proc_comm(j) .EQ. 1 .AND. &
        & proc_busy(iproc) .EQ. 0 .AND. &
        & proc_busy(jproc) .EQ. 0 )THEN

        icommseq(1,iproc,ncommseq) =   jproc + 1
        icommseq(1,jproc,ncommseq) = -(iproc + 1)    

        icommseq(2,iproc,ncommseq) = j
        icommseq(2,jproc,ncommseq) = j

        proc_comm(j)     = 0 ! communication is now listed
        proc_busy(iproc) = 1 ! processor busy in this sequence step
        proc_busy(jproc) = 1 ! processor busy in this sequence step

      END IF
    END DO

  END DO

  communication%ncommseq = ncommseq

  ALLOCATE(communication%icommseq(1:2,0:mpi_nproc-1,1:ncommseq))
  communication%icommseq = icommseq(1:2,0:mpi_nproc-1,1:ncommseq) 

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
  DEALLOCATE(proc_comm,stat=ierr)
  DEALLOCATE(icommseq,stat=ierr)

  END IF !mpi_rank0

!---------------------------------------------------------------------------------!
! Broadcast the communication setup
!---------------------------------------------------------------------------------!
  CALL MPI_BCAST(    communication%ncomm,1,mpi_prec_int,0, mpi_comm, ierr)
  CALL MPI_BCAST(    communication%ctype,1,mpi_prec_int,0, mpi_comm, ierr)
  CALL MPI_BCAST( communication%ncommseq,1,mpi_prec_int,0, mpi_comm, ierr)
  
  IF(mpi_rank .NE. 0)THEN
    ncomm = communication%ncomm
    ncommseq = communication%ncommseq

    ALLOCATE( communication%proc_send(ncomm) )
    ALLOCATE( communication%proc_recieve(ncomm) )

    ALLOCATE( communication%comm_dir(ncomm) )

    ALLOCATE( communication%npart_send(ncomm) )
    ALLOCATE( communication%npart_recieve(ncomm) )

    ALLOCATE(communication%icommseq(2,0:mpi_nproc-1,ncommseq))
  END IF

  CALL MPI_BCAST(communication%proc_send, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%proc_recieve, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)

  CALL MPI_BCAST(communication%comm_dir, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)

  CALL MPI_BCAST(communication%npart_send, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%npart_recieve, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)

  CALL MPI_BCAST( communication%icommseq,2*mpi_nproc*ncommseq, &
             & mpi_prec_int, 0, mpi_comm, ierr)

!---------------------------------------------------------------------------------!
! Allocate communication buffer
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(comm_buffer) ) THEN
    DEALLOCATE(comm_buffer)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
      GOTO 9999
    ENDIF
  END IF 
  ALLOCATE( comm_buffer(ncomm) )

  DO j = 1,ncomm
    iproc = communication%proc_send(j)
    jproc = communication%proc_recieve(j)

    npart_send = communication%npart_send(j)
 
    IF(mpi_rank .EQ. iproc)THEN
! Allocate send array for communication j
      ALLOCATE( comm_buffer(j)%part_send(nvar_send,npart_send) )
      comm_buffer(j)%nsend = nvar_send*npart_send
      comm_buffer(j)%npart_send = npart_send
    END IF

    IF(mpi_rank .EQ. jproc)THEN
! Allocate send array for communication j
      ALLOCATE( comm_buffer(j)%part_recieve(nvar_send,npart_send) )
      comm_buffer(j)%nrecieve = nvar_send*npart_send
      comm_buffer(j)%npart_recieve = npart_send
    END IF
  END DO

!---------------------------------------------------------------------------------!
! set counter for number of packed variables
!---------------------------------------------------------------------------------!
  communication%npacked = 0

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_particles_map
