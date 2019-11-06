!---------------------------------------------------------------------------------!
! pmlib_poisson_solve_mkl.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_poisson_solve(patch,topo_all,mesh,ierr,reg_vort,reproj)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_fourier
USE pmlib_mod_communication

USE MKL_DFTI

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                             :: patch
  TYPE(class_topology_all)                      :: topo_all
  TYPE(class_mesh)                              :: mesh
  INTEGER, INTENT(OUT)                          :: ierr
  LOGICAL, OPTIONAL, INTENT(IN)                 :: reg_vort
  LOGICAL, OPTIONAL, INTENT(IN)                 :: reproj

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                 :: caller = 'pmlib_poisson_solve'
  REAL(MK)                          :: PI
  INTEGER                           :: nreproj = 2

  INTEGER                           :: i,j,k,l

  REAL(MK)                          :: px, py, pz
  REAL(MK)                          :: rho

  INTEGER, DIMENSION(pmlib_ndim)    :: ncell
  INTEGER, DIMENSION(2*pmlib_ndim)  :: nghost
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin
  REAL(MK), DIMENSION(pmlib_ndim)   :: dx 

  INTEGER                           :: i_ptc, j_ptc, k_ptc
  REAL(MK), DIMENSION(pmlib_ndim)   :: l_ptc, ks, dk, kmin, ksmp
  INTEGER, DIMENSION(pmlib_ndim)    :: ncell_ptc
  REAL(MK)                          :: ki, kj, kk, s, kmag
  REAL(MK)                          :: eps
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin_ptc, xmax_ptc
  INTEGER, DIMENSION(pmlib_ndim)    :: bc_ptc

  REAL(MK)                          :: c1, c2, c3, c4, c5

  INTEGER, DIMENSION(pmlib_ndim)    :: fft_seq, ifft_seq
  LOGICAL                           :: ifft_vort

  COMPLEX(MKC)                      :: G, Kx, Ky, Kz
  COMPLEX(MKC)                      :: zeta_fft

  COMPLEX(MKC)                      :: psi1_fft, psi2_fft, psi3_fft
  COMPLEX(MKC)                      :: div_psi, div_vort

  INTEGER                           :: nfft

  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER :: vel_fft => NULL()
  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER :: vort_fft => NULL() 

  COMPLEX(MKC),DIMENSION(:),POINTER :: fft1_inout   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft1_inout  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: fft2_inout   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft2_inout  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: fft3_inout   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft3_inout  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft4_inout => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft5_inout => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft6_inout => NULL()

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_third => NULL()

  type(DFTI_DESCRIPTOR), POINTER    :: mkl_desc_handle

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

  ifft_vort = .FALSE.
!  IF(PRESENT(reg_vort) .AND. reg_vort .AND. patch%level .EQ. 1)THEN
  IF(PRESENT(reg_vort) .AND. reg_vort)THEN
     ifft_vort = .TRUE.
  END IF
!  IF(PRESENT(reproj) .AND. reproj .AND. patch%level .EQ. 1)THEN
  IF(PRESENT(reproj) .AND. reproj)THEN
     ifft_vort = .TRUE.
  END IF

!---------------------------------------------------------------------------------!
! Determine the fft sequence
!---------------------------------------------------------------------------------!
  ALLOCATE( topo_third(0:mpi_nproc-1) )

  fft_seq = patch%fft_seq

  IF(fft_seq(3) .EQ. 1)THEN
     DO i = 0,mpi_nproc-1
        topo_third(i) = topo_all%xpencil(i)
     END DO
  ELSEIF(fft_seq(3) .EQ. 2)THEN
     DO i = 0,mpi_nproc-1
        topo_third(i) = topo_all%ypencil(i)
     END DO
  ELSEIF(fft_seq(3) .EQ. 3)THEN
     DO i = 0,mpi_nproc-1
        topo_third(i) = topo_all%zpencil(i)
     END DO
  ELSE
     ierr = -1
     CALL pmlib_write(mpi_rank,caller, 'fft direction unknown.')
     GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Allocate Fourier coefficients for velocity field
!---------------------------------------------------------------------------------!
  ncell   = topo_third(mpi_rank)%ncell
  ALLOCATE( vel_fft(pmlib_ndim,ncell(1),ncell(2),ncell(3)) )

!---------------------------------------------------------------------------------!
! FFT the field in the first two directions and map to the third direction
!---------------------------------------------------------------------------------!
  fft_seq(3) = 0
  CALL pmlib_fft( topo_all,topo_all%cuboid,mesh%vort,topo_third,vort_fft, &
       & fft_seq,ierr,fft_shift = .FALSE.)

!---------------------------------------------------------------------------------!
! Setup the third direction
!---------------------------------------------------------------------------------!
  fft_seq  = patch%fft_seq
  ncell    = topo_third(mpi_rank)%ncell
  dx       = topo_third(mpi_rank)%dx
  xmin     = topo_third(mpi_rank)%xmin

  DO i = 1,pmlib_ndim
    IF( patch%bound_cond(i) .NE. 0 )THEN
      bc_ptc(i) = 1
    ELSE
      bc_ptc(i) = 0
    END IF
  END DO

! Include/exclude ghost cells
  xmin_ptc = patch%xmin - &
           & REAL(( 1 - bc_ptc ) * patch%nghost,MK)*patch%dx

!  IF(patch%level .EQ. 1)THEN
    ncell_ptc = (2 - bc_ptc) * (patch%ncell + 2 * ( 1 - bc_ptc ) * patch%nghost)
!  ELSE
!    ncell_ptc = patch%ncell + (1 - bc_ptc)*pmlib_patch_ext &
!            & + 2 * ( 1 - bc_ptc ) * patch%nghost)
!  END IF

  l_ptc = REAL(ncell_ptc,MK) * dx

!---------------------------------------------------------------------------------!
! Setup wavenumber parameters
!---------------------------------------------------------------------------------!
  dk(1) = 1.0_MK/l_ptc(1)
  dk(2) = 1.0_MK/l_ptc(2)
  dk(3) = 1.0_MK/l_ptc(3)

! patch index
  i_ptc = NINT((xmin(1) - xmin_ptc(1))/dx(1))
  j_ptc = NINT((xmin(2) - xmin_ptc(2))/dx(2))
  k_ptc = NINT((xmin(3) - xmin_ptc(3))/dx(3))

  ksmp(1) = 1.0_MK/(2.0_MK * dx(1))
  ksmp(2) = 1.0_MK/(2.0_MK * dx(2))
  ksmp(3) = 1.0_MK/(2.0_MK * dx(3))

  kmin(1) = - ksmp(1) + REAL(i_ptc,MK) * dk(1)
  kmin(2) = - ksmp(2) + REAL(j_ptc,MK) * dk(2)
  kmin(3) = - ksmp(3) + REAL(k_ptc,MK) * dk(3)

  eps = pmlib_regularisation_radius*MAXVAL(dx)

!---------------------------------------------------------------------------------!
! Setup coefficients for the regularisation filter
!---------------------------------------------------------------------------------!
  c1 = 1.0_MK
  c2 = 0.5_MK
  c3 = 0.125_MK
  c4 = 1.0_MK/48.0_MK
  c5 = 1.0_MK/384.0_MK

!---------------------------------------------------------------------------------!
! FFT the field in the third direction, convolve and IFFT
!---------------------------------------------------------------------------------!
  IF(fft_seq(3) .EQ. 1)THEN! Convolve on x-pencils

    IF(bc_ptc(1) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(1)
!      ELSE
!        nfft = ncell(1) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(1)
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

! Allocate fft pencils
     ALLOCATE( fft1_inout(nfft), ifft1_inout(nfft) )
     ALLOCATE( fft2_inout(nfft), ifft2_inout(nfft) )
     ALLOCATE( fft3_inout(nfft), ifft3_inout(nfft) )
     IF( ifft_vort )THEN
        ALLOCATE( ifft4_inout(nfft), ifft5_inout(nfft), ifft6_inout(nfft) )
     END IF

     DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
        IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
           kk = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
        ELSE
           kk = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
        END IF

        DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
           IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
              kj = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
           ELSE
              kj = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
           END IF

! Store fft pencil
           DO i = 1, nfft
              IF( i .LE. ncell(1) )THEN
                 fft1_inout(i) = vort_fft(1,i,j,k)
                 fft2_inout(i) = vort_fft(2,i,j,k)
                 fft3_inout(i) = vort_fft(3,i,j,k)
              ELSE! zero padding
                 fft1_inout(i) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft2_inout(i) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft3_inout(i) = CMPLX(0.0_MK,0.0_MK,MKC)
              END IF
           END DO

! Fourier transform (overwrites fft_inout)
           ierr = DftiComputeForward(mkl_desc_handle, fft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF

           DO i = 1,nfft
! Calculate wave number in x-direction (FFT-shifted)
              IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
                 ki = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
              ELSE
                 ki = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
              END IF

! Regularise output the vorticity field
              IF( ifft_vort )THEN
                 s = eps * 2.0_MK * PI * SQRT( ki**2 + kj**2 + kk**2 )

                 IF(pmlib_poisson_order .EQ. 0)THEN             
                    zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 2)THEN
                    zeta_fft = CMPLX( c1  &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 4)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 6)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 8)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 10)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                         & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 END IF

 ! Reguralise
                 ifft4_inout(i) = fft1_inout(i) * zeta_fft
                 ifft5_inout(i) = fft2_inout(i) * zeta_fft
                 ifft6_inout(i) = fft3_inout(i) * zeta_fft
              ELSE
                 zeta_fft = 1.0_MK
              END IF

              IF(pmlib_poisson_kernel .EQ. 1)THEN
! Get the Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                G = zeta_fft/(4.0_MK * PI**2 * kmag )
              ELSE
                G = 0.0_MK
              END IF

            ELSE

              G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kj,MKC) * fft3_inout(i) &
                                &                   - CMPLX(0.0_MK,kk,MKC) * fft2_inout(i) )
                 ifft2_inout(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kk,MKC) * fft1_inout(i) &
                                &                   - CMPLX(0.0_MK,ki,MKC) * fft3_inout(i) )
                 ifft3_inout(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC) * fft2_inout(i) &
                                &                   - CMPLX(0.0_MK,kj,MKC) * fft1_inout(i) )

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(i) &
                            &                   + CMPLX(0.0_MK,kj,MKC)*fft2_inout(i) &
                            &                   + CMPLX(0.0_MK,kk,MKC)*fft3_inout(i) )

                    ifft4_inout(i) = ifft4_inout(i) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
                    ifft5_inout(i) = ifft5_inout(i) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
                    ifft6_inout(i) = ifft6_inout(i) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kk,MKC) * div_psi
                 END IF
              ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
 ! Get the derivative Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                Kx = zeta_fft * CMPLX(0.0_MK,ki,MKC) /(2.0_MK * PI * kmag )
                Ky = zeta_fft * CMPLX(0.0_MK,kj,MKC) /(2.0_MK * PI * kmag )
                Kz = zeta_fft * CMPLX(0.0_MK,kk,MKC) /(2.0_MK * PI * kmag )
              ELSE
                Kx = 0.0_MK
                Ky = 0.0_MK
                Kz = 0.0_MK
              END IF

            ELSE

              Kx = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
              Ky = mesh%pkernel_fft(2,i,j,k) *dx(1)*dx(2)*dx(3)
              Kz = mesh%pkernel_fft(3,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(i) = Ky * fft3_inout(i) - Kz * fft2_inout(i)
                 ifft2_inout(i) = Kz * fft1_inout(i) - Kx * fft3_inout(i)
                 ifft3_inout(i) = Kx * fft2_inout(i) - Ky * fft1_inout(i)

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(i) &
                             &               + CMPLX(0.0_MK,kj,MKC)*fft2_inout(i) &
                             &               + CMPLX(0.0_MK,kk,MKC)*fft3_inout(i) )

                    ifft4_inout(i) = ifft4_inout(i) + Kx * div_vort
                    ifft5_inout(i) = ifft5_inout(i) + Ky * div_vort
                    ifft6_inout(i) = ifft6_inout(i) + Kz * div_vort

                 END IF
              END IF

           END DO

! Inverse fourier transform
           ierr = DftiComputeBackward(mkl_desc_handle, ifft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           IF( ifft_vort )THEN
              ierr = DftiComputeBackward(mkl_desc_handle, ifft4_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft5_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft6_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
           ENDIF

! Store fourier coefficients discard zero-padding
           DO i = 1,ncell(1)
              vel_fft(1,i,j,k) = ifft1_inout(i)/REAL(nfft,MK)
              vel_fft(2,i,j,k) = ifft2_inout(i)/REAL(nfft,MK)
              vel_fft(3,i,j,k) = ifft3_inout(i)/REAL(nfft,MK)
              IF( ifft_vort )THEN
                 vort_fft(1,i,j,k) = ifft4_inout(i)/REAL(nfft,MK)
                 vort_fft(2,i,j,k) = ifft5_inout(i)/REAL(nfft,MK)
                 vort_fft(3,i,j,k) = ifft6_inout(i)/REAL(nfft,MK)
              END IF
           END DO

        END DO!j
     END DO!k

  ELSEIF(fft_seq(3) .EQ. 2)THEN! Convolve on y-pencils

    IF(bc_ptc(2) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(2)
!      ELSE
!        nfft = ncell(2) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(2)
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

! Allocate fft pencils
     ALLOCATE( fft1_inout(nfft), ifft1_inout(nfft) )
     ALLOCATE( fft2_inout(nfft), ifft2_inout(nfft) )
     ALLOCATE( fft3_inout(nfft), ifft3_inout(nfft) )
     IF( ifft_vort )THEN
        ALLOCATE( ifft4_inout(nfft), ifft5_inout(nfft), ifft6_inout(nfft) )
     END IF

     DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
        IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
           kk = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
        ELSE
           kk = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
        END IF

        DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
           IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
              ki = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
           ELSE
              ki = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
           END IF

! Store fft pencil
           DO j = 1, nfft
              IF( j .LE. ncell(2) )THEN
                 fft1_inout(j) = vort_fft(1,i,j,k)
                 fft2_inout(j) = vort_fft(2,i,j,k)
                 fft3_inout(j) = vort_fft(3,i,j,k)
              ELSE! zero padding
                 fft1_inout(j) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft2_inout(j) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft3_inout(j) = CMPLX(0.0_MK,0.0_MK,MKC)
              END IF
           END DO

! Fourier transform (overwrites fft_inout)
           ierr = DftiComputeForward(mkl_desc_handle, fft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF

           DO j = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
              IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
                 kj = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
              ELSE
                 kj = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
              END IF

! regularise the output vorticity field
              IF( ifft_vort )THEN
                 s = eps * 2.0_MK * PI * SQRT( ki**2 + kj**2 + kk**2 )

                 IF(pmlib_poisson_order .EQ. 0)THEN             
                    zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 2)THEN
                    zeta_fft = CMPLX( c1  &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 4)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 6)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 8)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 10)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                         & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 END IF

 ! Reguralise
                 ifft4_inout(j) = fft1_inout(j) * zeta_fft
                 ifft5_inout(j) = fft2_inout(j) * zeta_fft
                 ifft6_inout(j) = fft3_inout(j) * zeta_fft
              ELSE
                 zeta_fft = 1.0_MK
              END IF

              IF(pmlib_poisson_kernel .EQ. 1)THEN
! Get the Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                G = zeta_fft/(4.0_MK * PI**2 * kmag )
              ELSE
                G = 0.0_MK
              END IF

            ELSE

              G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(j) = G * ( CMPLX(0.0_MK,kj,MKC) * fft3_inout(j) &
                      &     - CMPLX(0.0_MK,kk,MKC) * fft2_inout(j) )
                 ifft2_inout(j) = G * ( CMPLX(0.0_MK,kk,MKC) * fft1_inout(j) &
                      &     - CMPLX(0.0_MK,ki,MKC) * fft3_inout(j) )
                 ifft3_inout(j) = G * ( CMPLX(0.0_MK,ki,MKC) * fft2_inout(j) &
                      &     - CMPLX(0.0_MK,kj,MKC) * fft1_inout(j) )

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(j) &
                            &                   + CMPLX(0.0_MK,kj,MKC)*fft2_inout(j) &
                            &                   + CMPLX(0.0_MK,kk,MKC)*fft3_inout(j) )

                    ifft4_inout(j) = ifft4_inout(j) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
                    ifft5_inout(j) = ifft5_inout(j) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
                    ifft6_inout(j) = ifft6_inout(j) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kk,MKC) * div_psi
                 END IF
              ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
 ! Get the derivative Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                Kx = zeta_fft * CMPLX(0.0_MK,ki,MKC) /(2.0_MK * PI * kmag )
                Ky = zeta_fft * CMPLX(0.0_MK,kj,MKC) /(2.0_MK * PI * kmag )
                Kz = zeta_fft * CMPLX(0.0_MK,kk,MKC) /(2.0_MK * PI * kmag )
              ELSE
                Kx = 0.0_MK
                Ky = 0.0_MK
                Kz = 0.0_MK
              END IF

            ELSE

              Kx = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
              Ky = mesh%pkernel_fft(2,i,j,k) *dx(1)*dx(2)*dx(3)
              Kz = mesh%pkernel_fft(3,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(j) = Ky * fft3_inout(j) - Kz * fft2_inout(j)
                 ifft2_inout(j) = Kz * fft1_inout(j) - Kx * fft3_inout(j)
                 ifft3_inout(j) = Kx * fft2_inout(j) - Ky * fft1_inout(j)

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(j) &
                             &               + CMPLX(0.0_MK,kj,MKC)*fft2_inout(j) &
                             &               + CMPLX(0.0_MK,kk,MKC)*fft3_inout(j) )

                    ifft4_inout(j) = ifft4_inout(j) + Kx * div_vort
                    ifft5_inout(j) = ifft5_inout(j) + Ky * div_vort
                    ifft6_inout(j) = ifft6_inout(j) + Kz * div_vort
                 END IF
              END IF

           END DO

! Inverse fourier transform
           ierr = DftiComputeBackward(mkl_desc_handle, ifft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           IF( ifft_vort )THEN
              ierr = DftiComputeBackward(mkl_desc_handle, ifft4_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft5_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft6_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
           ENDIF

! Store fourier coefficients discard zero-padding
           DO j = 1,ncell(2)
              vel_fft(1,i,j,k) = ifft1_inout(j)/REAL(nfft,MK)
              vel_fft(2,i,j,k) = ifft2_inout(j)/REAL(nfft,MK)
              vel_fft(3,i,j,k) = ifft3_inout(j)/REAL(nfft,MK)
              IF( ifft_vort )THEN
                 vort_fft(1,i,j,k) = ifft4_inout(j)/REAL(nfft,MK)
                 vort_fft(2,i,j,k) = ifft5_inout(j)/REAL(nfft,MK)
                 vort_fft(3,i,j,k) = ifft6_inout(j)/REAL(nfft,MK)
              END IF
           END DO

        END DO!i
     END DO!k


  ELSEIF(fft_seq(3) .EQ. 3)THEN! Convolve on z-pencils

    IF(bc_ptc(3) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(3)
!      ELSE
!        nfft = ncell(3) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(3)
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

! Allocate fft pencils
     ALLOCATE( fft1_inout(nfft), ifft1_inout(nfft) )
     ALLOCATE( fft2_inout(nfft), ifft2_inout(nfft) )
     ALLOCATE( fft3_inout(nfft), ifft3_inout(nfft) )
     IF( ifft_vort )THEN
        ALLOCATE( ifft4_inout(nfft), ifft5_inout(nfft), ifft6_inout(nfft) )
     END IF

     DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
        IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
           kj = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
        ELSE
           kj = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
        END IF

        DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
           IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
              ki = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
           ELSE
              ki = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
           END IF

! Store fft pencil
           DO k = 1, nfft
              IF( k .LE. ncell(3) )THEN
                 fft1_inout(k) = vort_fft(1,i,j,k)
                 fft2_inout(k) = vort_fft(2,i,j,k)
                 fft3_inout(k) = vort_fft(3,i,j,k)
              ELSE! zero padding
                 fft1_inout(k) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft2_inout(k) = CMPLX(0.0_MK,0.0_MK,MKC)
                 fft3_inout(k) = CMPLX(0.0_MK,0.0_MK,MKC)
              END IF
           END DO

! Fourier transform (overwrites fft_inout)
           ierr = DftiComputeForward(mkl_desc_handle, fft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeForward(mkl_desc_handle, fft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF

           DO k = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
              IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
                 kk = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
              ELSE
                 kk = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
              END IF

! Regularise the output vorticity field
              IF( ifft_vort )THEN
                 s = eps * 2.0_MK * PI * SQRT( ki**2 + kj**2 + kk**2 )

                 IF(pmlib_poisson_order .EQ. 0)THEN             
                    zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 2)THEN
                    zeta_fft = CMPLX( c1  &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 4)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 6)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 8)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                         & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 ELSEIF(pmlib_poisson_order .EQ. 10)THEN
                    zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                         & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
                 END IF

 ! Reguralise
                 ifft4_inout(k) = fft1_inout(k) * zeta_fft
                 ifft5_inout(k) = fft2_inout(k) * zeta_fft
                 ifft6_inout(k) = fft3_inout(k) * zeta_fft
              ELSE
                 zeta_fft = 1.0_MK
              END IF

              IF(pmlib_poisson_kernel .EQ. 1)THEN
! Get the Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                G = zeta_fft/(4.0_MK * PI**2 * kmag )
              ELSE
                G = 0.0_MK
              END IF

            ELSE

              G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kj,MKC) * fft3_inout(k) &
                                &                   - CMPLX(0.0_MK,kk,MKC) * fft2_inout(k) )
                 ifft2_inout(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kk,MKC) * fft1_inout(k) &
                                &                   - CMPLX(0.0_MK,ki,MKC) * fft3_inout(k) )
                 ifft3_inout(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC) * fft2_inout(k) &
                                &                   - CMPLX(0.0_MK,kj,MKC) * fft1_inout(k) )

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(k) &
                            &                   + CMPLX(0.0_MK,kj,MKC)*fft2_inout(k) &
                            &                   + CMPLX(0.0_MK,kk,MKC)*fft3_inout(k) )

                    ifft4_inout(k) = ifft4_inout(k) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
                    ifft5_inout(k) = ifft5_inout(k) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
                    ifft6_inout(k) = ifft6_inout(k) &
                                 & + 2.0_MK * PI * CMPLX(0.0_MK,kk,MKC) * div_psi
                 END IF
              ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
 ! Get the derivative Poisson integration kernel
            IF( patch%bound_cond(1) .EQ. 1 .AND. &
              & patch%bound_cond(2) .EQ. 1 .AND. &
              & patch%bound_cond(3) .EQ. 1 )THEN

              kmag = ki**2 + kj**2 + kk**2
              IF(kmag .GT. 0.0_MK)THEN
                Kx = zeta_fft * CMPLX(0.0_MK,ki,MKC) /(2.0_MK * PI * kmag )
                Ky = zeta_fft * CMPLX(0.0_MK,kj,MKC) /(2.0_MK * PI * kmag )
                Kz = zeta_fft * CMPLX(0.0_MK,kk,MKC) /(2.0_MK * PI * kmag )
              ELSE
                Kx = 0.0_MK
                Ky = 0.0_MK
                Kz = 0.0_MK
              END IF

            ELSE

              Kx = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
              Ky = mesh%pkernel_fft(2,i,j,k) *dx(1)*dx(2)*dx(3)
              Kz = mesh%pkernel_fft(3,i,j,k) *dx(1)*dx(2)*dx(3)

            END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
                 ifft1_inout(k) = Ky * fft3_inout(k) - Kz * fft2_inout(k)
                 ifft2_inout(k) = Kz * fft1_inout(k) - Kx * fft3_inout(k)
                 ifft3_inout(k) = Kx * fft2_inout(k) - Ky * fft1_inout(k)

! Correct the output vorticity by the divergence of the stream function
                 IF(PRESENT(reproj) .AND. reproj)THEN
                    div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_inout(k) &
                             &               + CMPLX(0.0_MK,kj,MKC)*fft2_inout(k) &
                             &               + CMPLX(0.0_MK,kk,MKC)*fft3_inout(k) )

                    ifft4_inout(k) = ifft4_inout(k) + Kx * div_vort
                    ifft5_inout(k) = ifft5_inout(k) + Ky * div_vort
                    ifft6_inout(k) = ifft6_inout(k) + Kz * div_vort
                 END IF
              END IF

           END DO

! Inverse fourier transform
           ierr = DftiComputeBackward(mkl_desc_handle, ifft1_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft2_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           ierr = DftiComputeBackward(mkl_desc_handle, ifft3_inout)
           IF (ierr .NE. 0) THEN
              CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
              GOTO 9999
           ENDIF
           IF( ifft_vort )THEN
              ierr = DftiComputeBackward(mkl_desc_handle, ifft4_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft5_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
              ierr = DftiComputeBackward(mkl_desc_handle, ifft6_inout)
              IF (ierr .NE. 0) THEN
                 CALL pmlib_write(mpi_rank,caller,'Failed to fft.')
                 GOTO 9999
              ENDIF
           ENDIF

! Store fourier coefficients discard zero-padding
           DO k = 1,ncell(3)
              vel_fft(1,i,j,k) = ifft1_inout(k)/REAL(nfft,MK)
              vel_fft(2,i,j,k) = ifft2_inout(k)/REAL(nfft,MK)
              vel_fft(3,i,j,k) = ifft3_inout(k)/REAL(nfft,MK)
              IF( ifft_vort )THEN
                 vort_fft(1,i,j,k) = ifft4_inout(k)/REAL(nfft,MK)
                 vort_fft(2,i,j,k) = ifft5_inout(k)/REAL(nfft,MK)
                 vort_fft(3,i,j,k) = ifft6_inout(k)/REAL(nfft,MK)
              END IF
           END DO

        END DO!i
     END DO!j

  END IF

!---------------------------------------------------------------------------------!
! Clean up mkl
!---------------------------------------------------------------------------------!
  ierr = DftiFreeDescriptor(mkl_desc_handle)
  IF (ierr .NE. 0) THEN
     CALL pmlib_write(mpi_rank,caller,'Failed to free fft.')
     GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! De-allocate fft arrays
!---------------------------------------------------------------------------------!
  DEALLOCATE( fft1_inout, ifft1_inout )
  DEALLOCATE( fft2_inout, ifft2_inout )
  DEALLOCATE( fft3_inout, ifft3_inout )

  IF( ifft_vort )THEN
     DEALLOCATE( ifft4_inout, ifft5_inout, ifft6_inout )
  END IF

!---------------------------------------------------------------------------------!
! IFFT velocity field in the first direction
!---------------------------------------------------------------------------------!
  ifft_seq(1) = patch%fft_seq(2)
  ifft_seq(2) = patch%fft_seq(1)
  ifft_seq(3) = 0
  CALL pmlib_ifft( topo_all,topo_third,vel_fft,topo_all%cuboid,mesh%vel, &
       & ifft_seq,ierr,fft_shift = .FALSE.)

  IF( ifft_vort )THEN
     CALL pmlib_ifft( topo_all,topo_third,vort_fft,topo_all%cuboid,mesh%vort, &
          & ifft_seq,ierr,fft_shift = .FALSE.)
  END IF

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(vort_fft) ) THEN
     DEALLOCATE(vort_fft)
  END IF
  IF( ASSOCIATED(vel_fft) ) THEN
     DEALLOCATE(vel_fft)
  END IF

!---------------------------------------------------------------------------------!
! Deallocate
!---------------------------------------------------------------------------------!
  DEALLOCATE(topo_third,stat=ierr)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
9999 CONTINUE
  RETURN
END SUBROUTINE pmlib_poisson_solve
