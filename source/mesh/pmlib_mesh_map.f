!---------------------------------------------------------------------------------!
! pmlib_mesh_map.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_mesh_map(topo_send,topo_recieve,nvar_send,ierr,map_ghost)

!---------------------------------------------------------------------------------!
! Sets up the inter topology mesh communication
!  - Finds the overlapping mesh indices
!  - Finds the communication sequence
!  - Allocates the send and recieve buffers
!
! The communication is calculated on processor 0 only
! after which it is broadcasted to the remaining processors.
!---------------------------------------------------------------------------------!

USE pmlib_mod_communication
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_send
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_recieve
  INTEGER, INTENT(IN)                       :: nvar_send
  INTEGER, INTENT(OUT)                      :: ierr
  LOGICAL, INTENT(IN), OPTIONAL             :: map_ghost

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=14)                     :: caller = 'pmlib_mesh_map'
  INTEGER                               :: ncomm, icomm
  INTEGER                               :: iproc, jproc
  INTEGER                               :: i,j,k

  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin_send, xmax_send
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin_recieve, xmax_recieve
  REAL(MK), DIMENSION(pmlib_ndim)       :: dx
  INTEGER,DIMENSION(2*pmlib_ndim)       :: nghost_recieve, nghost_send
  INTEGER,DIMENSION(2*pmlib_ndim)       :: bc_send, bc_recieve
  INTEGER                               :: nghost
  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_send, ncell_recieve

  INTEGER, DIMENSION(0:mpi_nproc-1)     :: proc_busy 
  INTEGER, DIMENSION(:), POINTER        :: proc_comm

  INTEGER                               :: ncommseq
  INTEGER, DIMENSION(:,:,:), POINTER    :: icommseq 

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find number of communications
!---------------------------------------------------------------------------------!
  IF (mpi_rank .EQ. 0) THEN

  icomm = 0
  DO iproc = 0,mpi_nproc-1
    DO jproc = 0,mpi_nproc-1
!---------------------------------------------------------------------------------!
! Extend of subdomains
!---------------------------------------------------------------------------------!
      xmin_send      = topo_send(iproc)%xmin
      xmax_send      = topo_send(iproc)%xmax

      xmin_recieve   = topo_recieve(jproc)%xmin
      xmax_recieve   = topo_recieve(jproc)%xmax

      dx = topo_recieve(jproc)%dx
      IF( dx(1) .NE. topo_send(iproc)%dx(1) .OR. &
        & dx(2) .NE. topo_send(iproc)%dx(2) .OR. &
        & dx(3) .NE. topo_send(iproc)%dx(3))THEN
        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Mesh cell size does not coinside.')
        GOTO 9999
      END IF

!---------------------------------------------------------------------------------!
! Check overlap
!---------------------------------------------------------------------------------!
      IF( NINT( (xmax_recieve(1) - xmin_send(1))/dx(1) ) .GT. 0 .AND. &
        & NINT( (xmax_recieve(2) - xmin_send(2))/dx(2) ) .GT. 0 .AND. &
        & NINT( (xmax_recieve(3) - xmin_send(3))/dx(3) ) .GT. 0 .AND. &
        & NINT( (xmax_send(1) - xmin_recieve(1))/dx(1) ) .GT. 0 .AND. &
        & NINT( (xmax_send(2) - xmin_recieve(2))/dx(2) ) .GT. 0 .AND. &
        & NINT( (xmax_send(3) - xmin_recieve(3))/dx(3) ) .GT. 0 )THEN

        icomm = icomm + 1

      END IF ! intersecting

    END DO !jproc
  END DO !iproc

!---------------------------------------------------------------------------------!
! Allocate communication lists
!---------------------------------------------------------------------------------!
  ncomm               = icomm
  communication%ncomm = ncomm
  communication%ctype = 1

  ALLOCATE( communication%proc_send(ncomm) )
  ALLOCATE( communication%proc_recieve(ncomm) )

  ALLOCATE( communication%imin_send(pmlib_ndim,ncomm) )
  ALLOCATE( communication%imax_send(pmlib_ndim,ncomm) )
  ALLOCATE( communication%imin_recieve(pmlib_ndim,ncomm) )
  ALLOCATE( communication%imax_recieve(pmlib_ndim,ncomm) )

  ALLOCATE( communication%ncell_send(pmlib_ndim,ncomm) )
  ALLOCATE( communication%ncell_recieve(pmlib_ndim,ncomm) )

  ALLOCATE( proc_comm(ncomm) )
  proc_comm = 0

!---------------------------------------------------------------------------------!
! Loop through subdomains
!---------------------------------------------------------------------------------!
  icomm = 0
  DO iproc = 0,mpi_nproc-1
    DO jproc = 0,mpi_nproc-1
!---------------------------------------------------------------------------------!
! Subdomain parameters
!---------------------------------------------------------------------------------!
    xmin_send      = topo_send(iproc)%xmin
    xmax_send      = topo_send(iproc)%xmax

    xmin_recieve   = topo_recieve(jproc)%xmin
    xmax_recieve   = topo_recieve(jproc)%xmax

    ncell_send     = topo_send(iproc)%ncell
    ncell_recieve  = topo_recieve(jproc)%ncell

    nghost_send    = topo_send(iproc)%nghost
    nghost_recieve = topo_recieve(jproc)%nghost

    bc_send        = topo_send(iproc)%bound_cond
    bc_recieve     = topo_recieve(jproc)%bound_cond

    dx = topo_recieve(jproc)%dx

!---------------------------------------------------------------------------------!
! Check overlap
!---------------------------------------------------------------------------------!
    IF( NINT( (xmax_recieve(1) - xmin_send(1))/dx(1) ) .GT. 0 .AND. &
      & NINT( (xmax_recieve(2) - xmin_send(2))/dx(2) ) .GT. 0 .AND. &
      & NINT( (xmax_recieve(3) - xmin_send(3))/dx(3) ) .GT. 0 .AND. &
      & NINT( (xmax_send(1) - xmin_recieve(1))/dx(1) ) .GT. 0 .AND. &
      & NINT( (xmax_send(2) - xmin_recieve(2))/dx(2) ) .GT. 0 .AND. &
      & NINT( (xmax_send(3) - xmin_recieve(3))/dx(3) ) .GT. 0 )THEN

!---------------------------------------------------------------------------------!
! Count communications and tag the inter-processor communications
!---------------------------------------------------------------------------------!
    icomm = icomm + 1
    IF(iproc .EQ. jproc)THEN
      proc_comm(icomm) = 0
    ELSE
      proc_comm(icomm) = 1
    END IF

!---------------------------------------------------------------------------------!
! Store communicators
!---------------------------------------------------------------------------------!
    communication%proc_send(icomm)    = iproc
    communication%proc_recieve(icomm) = jproc

!---------------------------------------------------------------------------------!
! Find indices of intersection
!---------------------------------------------------------------------------------!
    DO i = 1,pmlib_ndim

! min index
      IF( NINT( (xmin_recieve(i) - xmin_send(i))/dx(i) ) .EQ. 0 )THEN

        IF(PRESENT(map_ghost) .AND. map_ghost .AND. & 
          & bc_recieve(2*i-1) .EQ. 0)THEN
          nghost = MIN( nghost_recieve(2*i-1), nghost_send(2*i-1) )
        ELSE
          nghost = 0
        END IF

        communication%imin_send(i,icomm)    = 1 - nghost
        communication%imin_recieve(i,icomm) = 1 - nghost

      ELSEIF( NINT( (xmin_send(i) - xmin_recieve(i))/dx(i) ) .GT. 0 )THEN

        IF(PRESENT(map_ghost) .AND. map_ghost .AND. & 
          & bc_send(2*i-1) .EQ. 0 )THEN
          nghost = MIN( nghost_send(2*i-1), &
                 & NINT((xmin_send(i)-xmin_recieve(i))/dx(i)) &
                 & + nghost_recieve(2*i-1) )
        ELSE
          nghost = 0
        END IF

        communication%imin_send(i,icomm)    = 1 - nghost
        communication%imin_recieve(i,icomm) = 1 & 
              & + NINT((xmin_send(i)-xmin_recieve(i))/dx(i)) - nghost

      ELSEIF( NINT( (xmin_recieve(i) - xmin_send(i))/dx(i) ) .GT. 0 )THEN

        IF( PRESENT(map_ghost) .AND. map_ghost .AND. & 
          & bc_recieve(2*i-1) .EQ. 0 )THEN
          nghost = MIN( nghost_recieve(2*i-1), &
                 & NINT((xmin_recieve(i)-xmin_send(i))/dx(i)) &
                 & + nghost_send(2*i-1) )
        ELSE
          nghost = 0
        END IF

        communication%imin_send(i,icomm)    = 1 & 
              & + NINT((xmin_recieve(i)-xmin_send(i))/dx(i)) - nghost
        communication%imin_recieve(i,icomm) = 1 - nghost

      ELSE
        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Failed to identify intersection indexes')
        GOTO 9999
      ENDIF


! max index
      IF( NINT( (xmax_recieve(i) - xmax_send(i))/dx(i) ) .EQ. 0 )THEN

        IF(PRESENT(map_ghost) .AND. map_ghost .AND. &
          & bc_recieve(2*i) .EQ. 0 )THEN
          nghost = MIN( nghost_recieve(2*i), nghost_send(2*i) )
        ELSE
          nghost = 0
        END IF

        communication%imax_send(i,icomm)    = ncell_send(i) + nghost
        communication%imax_recieve(i,icomm) = ncell_recieve(i) + nghost

      ELSEIF( NINT( (xmax_recieve(i) - xmax_send(i))/dx(i) ) .GT. 0 )THEN

        IF( PRESENT(map_ghost) .AND. map_ghost .AND. &
          & bc_send(2*i) .EQ. 0 )THEN
          nghost = MIN( nghost_send(2*i), &
                 & NINT((xmax_recieve(i)-xmax_send(i))/dx(i)) &
                 & + nghost_recieve(2*i) )
        ELSE
          nghost = 0
        END IF

        communication%imax_send(i,icomm)    = ncell_send(i) + nghost
        communication%imax_recieve(i,icomm) = & 
              & NINT((xmax_send(i)-xmin_recieve(i))/dx(i)) + nghost


      ELSEIF( NINT( (xmax_send(i) - xmax_recieve(i))/dx(i) ) .GT. 0 )THEN

        IF( PRESENT(map_ghost) .AND. map_ghost .AND. &
          & bc_recieve(2*i) .EQ. 0 )THEN
          nghost = MIN( nghost_recieve(2*i), &
                 & NINT((xmax_send(i)-xmax_recieve(i))/dx(i)) &
                 & + nghost_send(2*i) )
        ELSE
          nghost = 0
        END IF

        communication%imax_send(i,icomm) = & 
              & NINT((xmax_recieve(i)-xmin_send(i))/dx(i)) + nghost
        communication%imax_recieve(i,icomm) = ncell_recieve(i) + nghost

      ELSE
        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Failed to identify intersection indexes')
        GOTO 9999
      ENDIF

      communication%ncell_send(i,icomm) = &
                   & communication%imax_send(i,icomm) - &
                   &  communication%imin_send(i,icomm) + 1
      
      communication%ncell_recieve(i,icomm) = &
                   & communication%imax_recieve(i,icomm) - &
                   &  communication%imin_recieve(i,icomm) + 1

! check if the number of sended and recieved cells coinside
      IF(  communication%ncell_send(i,icomm) .NE. &
        &  communication%ncell_recieve(i,icomm))THEN
        ierr = -1
        CALL pmlib_write(mpi_rank,caller,'Number of sended and recieved cells does not correspond')
        GOTO 9999
      END IF

      IF(communication%imin_send(i,icomm) .LT. 1-nghost_send(2*i-1))THEN
        CALL pmlib_write(mpi_rank,caller,'Sending index is too low')
      END IF
      IF(communication%imax_send(i,icomm) .GT. ncell_send(i)+nghost_send(2*i))THEN
        CALL pmlib_write(mpi_rank,caller,'Sending index is too high')
      END IF
      IF(communication%imin_recieve(i,icomm) .LT. 1-nghost_recieve(2*i-1))THEN
        CALL pmlib_write(mpi_rank,caller,'Recieving index is too low')
      END IF
      IF(communication%imax_recieve(i,icomm) .GT. ncell_recieve(i)+nghost_recieve(2*i))THEN
        CALL pmlib_write(mpi_rank,caller,'Recieving index is too high')
      END IF

    END DO ! direction

  END IF ! intersecting

    END DO !jproc
  END DO !iproc

!---------------------------------------------------------------------------------!
! Check number of communications
!---------------------------------------------------------------------------------!
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
      IF( proc_comm(j) .EQ. 1 .AND. proc_busy(iproc) .EQ. 0 .AND. &
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

    ALLOCATE( communication%imin_send(pmlib_ndim,ncomm) )
    ALLOCATE( communication%imax_send(pmlib_ndim,ncomm) )
    ALLOCATE( communication%imin_recieve(pmlib_ndim,ncomm) )
    ALLOCATE( communication%imax_recieve(pmlib_ndim,ncomm) )

    ALLOCATE( communication%ncell_send(pmlib_ndim,ncomm) )
    ALLOCATE( communication%ncell_recieve(pmlib_ndim,ncomm) )

    ALLOCATE(communication%icommseq(2,0:mpi_nproc-1,ncommseq))
  END IF

  CALL MPI_BCAST(communication%proc_send, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%proc_recieve, ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)

  CALL MPI_BCAST(communication%imin_send, pmlib_ndim*ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%imax_send, pmlib_ndim*ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%imin_recieve, pmlib_ndim*ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%imax_recieve, pmlib_ndim*ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)

  CALL MPI_BCAST(communication%ncell_send, pmlib_ndim*ncomm, & 
             & mpi_prec_int, 0, mpi_comm, ierr)
  CALL MPI_BCAST(communication%ncell_recieve, pmlib_ndim*ncomm, & 
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

    ncell_send(1) = communication%ncell_send(1,j)  
    ncell_send(2) = communication%ncell_send(2,j)  
    ncell_send(3) = communication%ncell_send(3,j) 
 
    IF(mpi_rank .EQ. iproc)THEN
! Allocate send array for communication j
      ALLOCATE( comm_buffer(j)%mesh_send(nvar_send, &
              & ncell_send(1),ncell_send(2),ncell_send(3)) )
      comm_buffer(j)%nsend = nvar_send*ncell_send(1) &
                         & * ncell_send(2)*ncell_send(3)
      comm_buffer(j)%ncell_send = ncell_send
    END IF

    IF(mpi_rank .EQ. jproc)THEN
! Allocate send array for communication j
      ALLOCATE( comm_buffer(j)%mesh_recieve(nvar_send, &
              & ncell_send(1),ncell_send(2),ncell_send(3)) )
      comm_buffer(j)%nrecieve = nvar_send*ncell_send(1) &
                            & * ncell_send(2)*ncell_send(3)
      comm_buffer(j)%ncell_recieve = ncell_send
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
END SUBROUTINE pmlib_mesh_map
