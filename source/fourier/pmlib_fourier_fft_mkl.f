!---------------------------------------------------------------------------------!
! pmlib_fourier_fft_mkl.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_fft( topo_all,topo_in,mesh,topo_out,mesh_fft,fft_seq, &
                    & ierr,fft_shift )

USE MKL_DFTI

  IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology_all)                        :: topo_all
  TYPE(class_topology),DIMENSION(:),POINTER       :: topo_in
  TYPE(class_topology),DIMENSION(:),POINTER       :: topo_out
  REAL(MK),DIMENSION(:,:,:,:),POINTER             :: mesh
  INTEGER, DIMENSION(pmlib_ndim), INTENT(IN)      :: fft_seq
  LOGICAL, OPTIONAL, INTENT(IN)                   :: fft_shift
  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER         :: mesh_fft
  INTEGER, INTENT(OUT)                            :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=10)                 :: caller = 'pmlib_fft'
  INTEGER                           :: i,j,k,l,m

  INTEGER, DIMENSION(pmlib_ndim)    :: ncell
  INTEGER, DIMENSION(2*pmlib_ndim)  :: nghost
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin, xmax, dx 

  INTEGER                           :: nvar
  INTEGER                           :: nfft

  COMPLEX(MKC),DIMENSION(:),POINTER :: fft_inout  => NULL()

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_send => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_recieve => NULL()

  type(DFTI_DESCRIPTOR), POINTER    :: mkl_desc_handle

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  nvar = SIZE(mesh,1)

!---------------------------------------------------------------------------------!
! Move the real mesh to the complex fft mesh
!---------------------------------------------------------------------------------!
  ncell  = topo_in(mpi_rank)%ncell
  nghost = topo_in(mpi_rank)%nghost

  IF( ASSOCIATED(mesh_fft) ) THEN
     DEALLOCATE(mesh_fft,stat=ierr)
     IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
        GOTO 9999
     ENDIF
  END IF

  ALLOCATE(mesh_fft(nvar, &
       & 1-nghost(1):ncell(1)+nghost(2), &
       & 1-nghost(3):ncell(2)+nghost(4), &
       & 1-nghost(5):ncell(3)+nghost(6)), &
       & stat=ierr)
  IF (ierr .NE. 0) THEN
     CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
     GOTO 9999
  ENDIF

  DO k = 1 - nghost(5), ncell(3) + nghost(6)
     DO j = 1 - nghost(3), ncell(2) + nghost(4)
        DO i = 1 - nghost(1), ncell(1) + nghost(2)
           DO l = 1,nvar
              mesh_fft(l,i,j,k) = CMPLX(mesh(l,i,j,k),0.0_MK,MKC)
           END DO
        END DO
     END DO
  END DO

  ALLOCATE(topo_send(0:mpi_nproc-1), topo_recieve(0:mpi_nproc-1))

  DO m = 1,pmlib_ndim
     IF(fft_seq(m) .EQ. 0)THEN
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

     IF(fft_seq(m) .EQ. 1)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%xpencil(i)
        END DO
     ELSEIF(fft_seq(m) .EQ. 2)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%ypencil(i)
        END DO
     ELSEIF(fft_seq(m) .EQ. 3)THEN
        DO i = 0,mpi_nproc-1
           topo_recieve(i) = topo_all%zpencil(i)
        END DO
     ELSE
        ierr = -1
        CALL pmlib_write(mpi_rank,caller, 'fft direction unknown.')
        GOTO 9999
     END IF

     CALL pmlib_mesh_map(topo_send,topo_recieve,2*nvar,ierr,map_ghost=.TRUE.)
     CALL pmlib_comm_pack(mesh_fft,ierr)
     CALL pmlib_comm_send(ierr)
     CALL pmlib_comm_unpack(topo_recieve,mesh_fft,0,ierr,clear=.TRUE.)
     CALL pmlib_comm_finalise(ierr)

!---------------------------------------------------------------------------------!
! Do FFTs
!---------------------------------------------------------------------------------!
     IF(fft_seq(m) .EQ. 1)THEN ! FFT x-pencils

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

        ALLOCATE(fft_inout(nfft))

        DO k = 1, ncell(3)
           DO j = 1, ncell(2)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO i = 1, nfft/2
                       fft_inout(i) = mesh_fft(l,i + nfft/2,j,k)
                    END DO
                    DO i = nfft/2 + 1, nfft
                       fft_inout(i) = mesh_fft(l,i - nfft/2,j,k)
                    END DO
                 ELSE
                    DO i = 1, nfft
                       fft_inout(i) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier transform (overwrites fft_inout)
                 ierr = DftiComputeForward(mkl_desc_handle, fft_inout)
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO i = 1, ncell(1)
                    mesh_fft(l,i,j,k) = fft_inout(i)
                 END DO

              END DO!l
           END DO!j
        END DO!k

     ELSEIF(fft_seq(m) .EQ. 2)THEN ! FFT y-pencils

        ncell = topo_recieve(mpi_rank)%ncell
        nfft  = ncell(2)

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

        ALLOCATE(fft_inout(nfft))

        DO k = 1, ncell(3)
           DO i = 1, ncell(1)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO j = 1, nfft/2
                       fft_inout(j) = mesh_fft(l,i,j + nfft/2,k)
                    END DO
                    DO j = nfft/2 + 1, nfft
                       fft_inout(j) = mesh_fft(l,i,j - nfft/2,k)
                    END DO
                 ELSE
                    DO j = 1, nfft
                       fft_inout(j) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier transform (overwrites fft_inout)
                 ierr = DftiComputeForward(mkl_desc_handle, fft_inout) 
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO j = 1, ncell(2)
                    mesh_fft(l,i,j,k) = fft_inout(j)
                 END DO

              END DO!l
           END DO!i
        END DO!k

     ELSEIF(fft_seq(m) .EQ. 3)THEN ! FFT z-pencils 

        ncell = topo_recieve(mpi_rank)%ncell
        nfft  = ncell(3)

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

        ALLOCATE(fft_inout(nfft))

        DO j = 1, ncell(2)
           DO i = 1, ncell(1)
              DO l = 1,nvar

! Store fft pencil
                 IF(PRESENT(fft_shift) .AND. fft_shift)THEN
                    DO k = 1, nfft/2
                       fft_inout(k) = mesh_fft(l,i,j,k + nfft/2)
                    END DO
                    DO k = nfft/2 + 1, nfft
                       fft_inout(k) = mesh_fft(l,i,j,k - nfft/2)
                    END DO
                 ELSE
                    DO k = 1, nfft
                       fft_inout(k) = mesh_fft(l,i,j,k)
                    END DO
                 END IF

! Fourier transform (overwrites fft_inout)
                 ierr = DftiComputeForward(mkl_desc_handle, fft_inout)
                 IF (ierr .NE. 0) THEN
                    CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                    GOTO 9999
                 ENDIF

! Store kernel
                 DO k = 1, ncell(3)
                    mesh_fft(l,i,j,k) = fft_inout(k)
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

     DEALLOCATE(fft_inout)

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
! Deallocate local pointers
!---------------------------------------------------------------------------------!
  DEALLOCATE(topo_send, topo_recieve)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
9999 CONTINUE
  RETURN
END SUBROUTINE pmlib_fft
