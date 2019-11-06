!---------------------------------------------------------------------------------!
! pmlib_fourier_ifft_mkl.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_ifft( topo_all,topo_in,mesh_fft,topo_out,mesh,ifft_seq, &
                     & ierr,fft_shift )

USE MKL_DFTI

  IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology_all)                        :: topo_all
  TYPE(class_topology),DIMENSION(:),POINTER       :: topo_in
  TYPE(class_topology),DIMENSION(:),POINTER       :: topo_out
  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER         :: mesh_fft
  INTEGER, DIMENSION(pmlib_ndim), INTENT(IN)      :: ifft_seq
  LOGICAL, OPTIONAL, INTENT(IN)                   :: fft_shift
  REAL(MK),DIMENSION(:,:,:,:),POINTER             :: mesh
  INTEGER, INTENT(OUT)                            :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=10)                 :: caller = 'pmlib_ifft'
  INTEGER                           :: i,j,k,l,m

  INTEGER, DIMENSION(pmlib_ndim)    :: ncell
  INTEGER, DIMENSION(2*pmlib_ndim)  :: nghost
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin, xmax, dx 

  INTEGER                           :: nvar
  INTEGER                           :: nfft

  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft_inout  => NULL()

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_send => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_recieve => NULL()

  type(DFTI_DESCRIPTOR), POINTER    :: mkl_desc_handle

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  nvar = SIZE(mesh,1)

!---------------------------------------------------------------------------------!
! Allocate auxilirary topologies
!---------------------------------------------------------------------------------!
  ALLOCATE(topo_send(0:mpi_nproc-1), topo_recieve(0:mpi_nproc-1))

  DO m = 1,pmlib_ndim
     IF(ifft_seq(m) .EQ. 0)THEN
        GOTO 100
     END IF

!---------------------------------------------------------------------------------!
! Map topology to the fft direction
!---------------------------------------------------------------------------------!
     IF(m .EQ. 1)THEN
        DO i = 0,mpi_nproc-1
           topo_send(i) = topo_in(i)
        END DO
     END IF

     IF(ifft_seq(m) .EQ. 1)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%xpencil(i)
        END DO
     ELSEIF(ifft_seq(m) .EQ. 2)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%ypencil(i)
        END DO
     ELSEIF(ifft_seq(m) .EQ. 3)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%zpencil(i)
        END DO
     ELSE
        ierr = -1
        CALL pmlib_write(mpi_rank,caller, 'ifft direction unknown.')
        GOTO 9999
     END IF

     CALL pmlib_mesh_map(topo_send,topo_recieve,2*nvar,ierr,map_ghost=.TRUE.)
     CALL pmlib_comm_pack(mesh_fft,ierr)
     CALL pmlib_comm_send(ierr)
     CALL pmlib_comm_unpack(topo_recieve,mesh_fft,0,ierr,clear=.TRUE.)
     CALL pmlib_comm_finalise(ierr)

!---------------------------------------------------------------------------------!
! Do IFFTs
!---------------------------------------------------------------------------------!
     IF(ifft_seq(m) .EQ. 1)THEN! IFFT x-pencils

        dx    = topo_recieve(mpi_rank)%dx
        ncell = topo_recieve(mpi_rank)%ncell
        nfft  = ncell(1)

        IF(MOD(nfft,2) .NE. 0)THEN
           ierr = -1
           CALL pmlib_write(mpi_rank,caller, 'nfft is an odd number.')
           GOTO 9999
        END IF

! Setup a complex to complex in-place transform
#ifdef __single
! Single precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_SINGLE, &
          & DFTI_COMPLEX, 1, nfft) 
#elif __double
! Double precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_DOUBLE, &
         & DFTI_COMPLEX, 1, nfft) 
#endif
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to create fft.')
           GOTO 9999
        ENDIF
        ierr = DftiCommitDescriptor(mkl_desc_handle)
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to commit fft.')
           GOTO 9999
        ENDIF

        ALLOCATE(ifft_inout(nfft))

        DO k = 1, ncell(3)
           DO j = 1, ncell(2)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO i = 1, nfft/2
                       ifft_inout(i) = mesh_fft(l,i + nfft/2,j,k)
                    END DO
                    DO i = nfft/2 + 1, nfft
                       ifft_inout(i) = mesh_fft(l,i - nfft/2,j,k)
                    END DO
                 ELSE
                    DO i = 1, nfft
                       ifft_inout(i) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier inverse transform
                 ierr = DftiComputeBackward(mkl_desc_handle, ifft_inout)! overwrite fft_inout
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO i = 1, ncell(1)
                    mesh_fft(l,i,j,k) = ifft_inout(i)/REAL(nfft,MK)
                 END DO

              END DO!l
           END DO!j
        END DO!k

     ELSEIF(ifft_seq(m) .EQ. 2)THEN! IFFT y-pencils

        dx     = topo_recieve(mpi_rank)%dx
        ncell  = topo_recieve(mpi_rank)%ncell
        nfft = ncell(2)

        IF(MOD(nfft,2) .NE. 0)THEN
           ierr = -1
           CALL pmlib_write(mpi_rank,caller, 'nfft is an odd number.')
           GOTO 9999
        END IF

! Setup a complex to complex in-place transform
#ifdef __single
! Single precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_SINGLE, &
          & DFTI_COMPLEX, 1, nfft) 
#elif __double
! Double precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_DOUBLE, &
         & DFTI_COMPLEX, 1, nfft) 
#endif
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to create fft.')
           GOTO 9999
        ENDIF
        ierr = DftiCommitDescriptor(mkl_desc_handle)
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to commit fft.')
           GOTO 9999
        ENDIF

        ALLOCATE(ifft_inout(nfft))

        DO k = 1, ncell(3)
           DO i = 1, ncell(1)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO j = 1, nfft/2
                       ifft_inout(j) = mesh_fft(l,i,j + nfft/2,k)
                    END DO
                    DO j = nfft/2 + 1, nfft
                       ifft_inout(j) = mesh_fft(l,i,j - nfft/2,k)
                    END DO
                 ELSE
                    DO j = 1, nfft
                       ifft_inout(j) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier inverse transform
                 ierr = DftiComputeBackward(mkl_desc_handle, ifft_inout)! overwrite fft_inout
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO j = 1, ncell(2)
                    mesh_fft(l,i,j,k) = ifft_inout(j)/REAL(nfft,MK)
                 END DO

              END DO!l
           END DO!i
        END DO!k

     ELSEIF(ifft_seq(m) .EQ. 3)THEN! FFT z-pencils 

        dx     = topo_recieve(mpi_rank)%dx
        ncell  = topo_recieve(mpi_rank)%ncell
        nfft = ncell(3)

        IF(MOD(nfft,2) .NE. 0)THEN
           ierr = -1
           CALL pmlib_write(mpi_rank,caller, 'nfft is an odd number.')
           GOTO 9999
        END IF

! Setup a complex to complex in-place transform
#ifdef __single
! Single precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_SINGLE, &
          & DFTI_COMPLEX, 1, nfft) 
#elif __double
! Double precision
    ierr = DftiCreateDescriptor(mkl_desc_handle, DFTI_DOUBLE, &
         & DFTI_COMPLEX, 1, nfft) 
#endif
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to create fft.')
           GOTO 9999
        ENDIF
        ierr = DftiCommitDescriptor(mkl_desc_handle)
        IF (ierr .NE. 0) THEN
           CALL pmlib_write(mpi_rank,caller,'Failed to commit fft.')
           GOTO 9999
        ENDIF

        ALLOCATE(ifft_inout(nfft))

        DO j = 1, ncell(2)
           DO i = 1, ncell(1)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO k = 1, nfft/2
                       ifft_inout(k) = mesh_fft(l,i,j,k + nfft/2)
                    END DO
                    DO k = nfft/2 + 1, nfft
                       ifft_inout(k) = mesh_fft(l,i,j,k - nfft/2)
                    END DO
                 ELSE
                    DO k = 1, nfft
                       ifft_inout(k) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier inverse transform
                 ierr = DftiComputeBackward(mkl_desc_handle, ifft_inout)! overwrite fft_inout
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO k = 1, ncell(3)
                    mesh_fft(l,i,j,k) = ifft_inout(k)/REAL(nfft,MK)
                 END DO

              END DO!l
           END DO!i
        END DO!j

     END IF! fft direction

! Clean up mkl
     ierr = DftiFreeDescriptor(mkl_desc_handle)
     IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to free fft.')
        GOTO 9999
     ENDIF

     DEALLOCATE(ifft_inout)

!---------------------------------------------------------------------------------!
! Shift topology reciever to sender
!---------------------------------------------------------------------------------!
     DO i = 0,mpi_nproc-1
        topo_send(i) = topo_recieve(i)
     END DO

  END DO! fft_seq
100 CONTINUE

!---------------------------------------------------------------------------------!
! Map to the output topology
!---------------------------------------------------------------------------------!
  CALL pmlib_mesh_map(topo_send,topo_out,2*nvar,ierr,map_ghost=.TRUE.)
  CALL pmlib_comm_pack(mesh_fft,ierr)
  CALL pmlib_comm_send(ierr)
  CALL pmlib_comm_unpack(topo_out,mesh_fft,0,ierr,clear=.TRUE.)
  CALL pmlib_comm_finalise(ierr)

!---------------------------------------------------------------------------------!
! Map boundaries
!---------------------------------------------------------------------------------!
  DO i = 1,pmlib_ndim
     CALL pmlib_mesh_map_ghost(topo_out,i,2*nvar,ierr,incl_edges=.FALSE.)
     CALL pmlib_comm_pack(mesh_fft,ierr)
     CALL pmlib_comm_send(ierr)
     CALL pmlib_comm_unpack(topo_out,mesh_fft,0,ierr,clear=.FALSE.)
     CALL pmlib_comm_finalise(ierr)
  END DO

  DEALLOCATE(topo_send, topo_recieve)
!---------------------------------------------------------------------------------!
! Allocate the real mesh and extract the real part of the complex mesh
!---------------------------------------------------------------------------------!
  ncell  = topo_out(mpi_rank)%ncell
  nghost = topo_out(mpi_rank)%nghost

  IF( ASSOCIATED(mesh) ) THEN
     DEALLOCATE(mesh,stat=ierr)
     IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
        GOTO 9999
     ENDIF
  END IF

  ALLOCATE(mesh(nvar, &
       & 1-nghost(1):ncell(1)+nghost(2), &
       & 1-nghost(3):ncell(2)+nghost(4), &
       & 1-nghost(5):ncell(3)+nghost(6)), &
       & stat=ierr)
  IF (ierr .NE. 0) THEN
     CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
     GOTO 9999
  ENDIF

! Extract the real part
  DO k = 1 - nghost(5), ncell(3) + nghost(6)
     DO j = 1 - nghost(3), ncell(2) + nghost(4)
        DO i = 1 - nghost(1), ncell(1) + nghost(2)
           DO l = 1,nvar
              mesh(l,i,j,k) = REAL(mesh_fft(l,i,j,k))
           END DO
        END DO
     END DO
  END DO

!---------------------------------------------------------------------------------!
! Deallocate the Fourier coefficients
!---------------------------------------------------------------------------------!
  DEALLOCATE(mesh_fft)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
9999 CONTINUE
  RETURN
END SUBROUTINE pmlib_ifft
