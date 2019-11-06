!---------------------------------------------------------------------------------!
! pmlib_mesh_map_ghost.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_mesh_map_ghost(topo,dir,nvar_send,ierr,incl_edges)

!---------------------------------------------------------------------------------!
! Sets up the mesh boundary communication
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
  TYPE(class_topology),DIMENSION(:),POINTER :: topo
  INTEGER, INTENT(IN)                       :: dir
  INTEGER, INTENT(IN)                       :: nvar_send
  INTEGER, INTENT(OUT)                      :: ierr
  LOGICAL, INTENT(IN), OPTIONAL             :: incl_edges

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                     :: caller = 'pmlib_mesh_map_ghost'
  INTEGER                               :: ncomm, icomm
  INTEGER                               :: iproc, jproc
  INTEGER                               :: i,j,k

  INTEGER                               :: ineigh, neigh
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin_send, xmax_send
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin_recieve, xmax_recieve
  REAL(MK), DIMENSION(pmlib_ndim)       :: dx

  INTEGER,DIMENSION(2*pmlib_ndim)       :: nghost_recieve, nghost_send
  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_send, ncell_recieve
  INTEGER                               :: nghost, ncell

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
  DO ineigh = 2*dir-1,2*dir
    DO iproc = 0,mpi_nproc-1
      neigh = topo(iproc)%ineigh(ineigh)

      IF( neigh .NE. -1 )THEN

        IF(PRESENT(incl_edges) .AND. incl_edges)THEN
          nghost_send    = topo(iproc)%nghost
          nghost_recieve = topo(neigh)%nghost

          IF( ( nghost_send(ineigh) .GT. 0 ) .OR. &
            & ( ineigh .EQ. 2*dir-1 .AND.  nghost_recieve(2*dir) .GT. 0 ) .OR. &
            & ( ineigh .EQ. 2*dir   .AND. nghost_recieve(2*dir-1) .GT. 0 ) )THEN
            icomm = icomm + 1
          END IF

        ELSE
          nghost_recieve = topo(neigh)%nghost

          IF( ( ineigh .EQ. 2*dir-1 .AND. nghost_recieve(2*dir) .GT. 0 ) .OR. &
            & ( ineigh .EQ. 2*dir .AND. nghost_recieve(2*dir-1) .GT. 0 ) )THEN
            icomm = icomm + 1
          END IF
        END IF

      END IF ! neighbours

    END DO !iproc
  END DO !ineigh

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
! Loop through neighbours
!---------------------------------------------------------------------------------!
  icomm = 0
  DO ineigh = 2*dir-1,2*dir
    DO iproc = 0,mpi_nproc-1
    
  neigh = topo(iproc)%ineigh(ineigh)

  IF( neigh .NE. -1 )THEN
    nghost_send    = topo(iproc)%nghost
    nghost_recieve = topo(neigh)%nghost

      IF(PRESENT(incl_edges) .AND. incl_edges)THEN
        IF( ( nghost_send(ineigh) .GT. 0 ) .OR. &
          & ( ineigh .EQ. 2*dir-1 .AND.   nghost_recieve(2*dir) .GT. 0 ) .OR. &
          & ( ineigh .EQ. 2*dir   .AND. nghost_recieve(2*dir-1) .GT. 0 ) )THEN
!---------------------------------------------------------------------------------!
! Count communications, tag the inter-processor communications
!---------------------------------------------------------------------------------!
          icomm = icomm + 1
          IF(iproc .NE. neigh)THEN
            proc_comm(icomm) = 1
          END IF

!---------------------------------------------------------------------------------!
! Store communicators
!---------------------------------------------------------------------------------!
          communication%proc_send(icomm)    = iproc
          communication%proc_recieve(icomm) = neigh

!---------------------------------------------------------------------------------!
! Store number of mesh cells in the topologies
!---------------------------------------------------------------------------------!
          ncell_send     = topo(iproc)%ncell
          ncell_recieve  = topo(neigh)%ncell

!---------------------------------------------------------------------------------!
! Find indices of ghost cells
!---------------------------------------------------------------------------------!
          DO i = 1,pmlib_ndim

            IF(2*i-1 .EQ. ineigh)THEN

              communication%imin_send(i,icomm)    = 1 - nghost_send(2*i-1)
              communication%imin_recieve(i,icomm) = ncell_recieve(i) + 1 &
                                                & - nghost_send(2*i-1)

              communication%imax_send(i,icomm)    = nghost_recieve(2*i)
              communication%imax_recieve(i,icomm) = ncell_recieve(i) &
                                                & + nghost_recieve(2*i)

            ELSEIF(2*i .EQ. ineigh)THEN

              communication%imin_send(i,icomm)    = ncell_send(i) &
                                                & - nghost_recieve(2*i-1) + 1
              communication%imin_recieve(i,icomm) = 1 - nghost_recieve(2*i-1)

              communication%imax_send(i,icomm)    = ncell_send(i) &
                                                & + nghost_send(2*i)
              communication%imax_recieve(i,icomm) = nghost_send(2*i)

            ELSE
              ncell  = MIN(ncell_send(i),ncell_recieve(i))

              nghost = MIN(nghost_send(2*i-1),nghost_recieve(2*i-1))
              communication%imin_send(i,icomm)    = 1 - nghost
              communication%imin_recieve(i,icomm) = 1 - nghost

              nghost = MIN(nghost_send(2*i),nghost_recieve(2*i))
              communication%imax_send(i,icomm)    = ncell + nghost
              communication%imax_recieve(i,icomm) = ncell + nghost
            END IF

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
              CALL pmlib_write( mpi_rank,caller, &
                  & 'Number of sended and recieved cells does not correspond')
              GOTO 9999
            END IF

          END DO ! ndim
        END IF ! ineigh

      ELSE
        IF( ( ineigh .EQ. 2*dir-1 .AND. nghost_recieve(2*dir) .GT. 0 ) .OR. &
          & ( ineigh .EQ. 2*dir   .AND. nghost_recieve(2*dir-1) .GT. 0 ) )THEN

!---------------------------------------------------------------------------------!
! Count communications
!---------------------------------------------------------------------------------!
          icomm = icomm + 1
          IF(iproc .NE. neigh)THEN
            proc_comm(icomm) = 1
          END IF

!---------------------------------------------------------------------------------!
! Store communicators
!---------------------------------------------------------------------------------!
          communication%proc_send(icomm)    = iproc
          communication%proc_recieve(icomm) = neigh

!---------------------------------------------------------------------------------!
! Store number of mesh cells in the topologies
!---------------------------------------------------------------------------------!
          ncell_send    = topo(iproc)%ncell
          ncell_recieve = topo(neigh)%ncell

!---------------------------------------------------------------------------------!
! Find indices of ghost cells
!---------------------------------------------------------------------------------!
          DO i = 1,pmlib_ndim

            IF(2*i-1 .EQ. ineigh)THEN

              communication%imin_send(i,icomm)    = 1
              communication%imin_recieve(i,icomm) = ncell_recieve(i) + 1

              communication%imax_send(i,icomm)    = nghost_recieve(2*i)
              communication%imax_recieve(i,icomm) = ncell_recieve(i) &
                                                & + nghost_recieve(2*i)

            ELSEIF(2*i .EQ. ineigh)THEN

              communication%imin_send(i,icomm)    = ncell_send(i) &
                                                & - nghost_recieve(2*i-1) + 1
              communication%imin_recieve(i,icomm) = 1 - nghost_recieve(2*i-1)

              communication%imax_send(i,icomm)    = ncell_send(i)
              communication%imax_recieve(i,icomm) = 0

            ELSE
              ncell  = MIN(ncell_send(i),ncell_recieve(i))

              nghost = MIN(nghost_send(2*i-1),nghost_recieve(2*i-1))
              communication%imin_send(i,icomm)    = 1 - nghost
              communication%imin_recieve(i,icomm) = 1 - nghost

              nghost = MIN(nghost_send(2*i),nghost_recieve(2*i))
              communication%imax_send(i,icomm)    = ncell + nghost
              communication%imax_recieve(i,icomm) = ncell + nghost
            END IF

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
              CALL pmlib_write( mpi_rank,caller, &
                  & 'Number of sended and recieved cells does not correspond')
              GOTO 9999
            END IF

          END DO ! ndim
        END IF ! ineigh

      END IF ! incl_edges

    END IF ! neighbours

    END DO !iproc
  END DO !ineigh

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
END SUBROUTINE pmlib_mesh_map_ghost
