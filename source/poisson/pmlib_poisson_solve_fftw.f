!---------------------------------------------------------------------------------!
! pmlib_poisson_solve_fftw.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_poisson_solve(patch,topo_all,mesh,ierr,reg_vort,reproj)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_fourier
USE pmlib_mod_communication

IMPLICIT NONE

INCLUDE 'fftw3.f'
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
  REAL(MK)                          :: ki, kj, kk, s
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

  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER :: vel_fft, vort_fft 

  COMPLEX(MKC),DIMENSION(:),POINTER :: fft1_in   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: fft1_out  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft1_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft1_out => NULL()

  COMPLEX(MKC),DIMENSION(:),POINTER :: fft2_in   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: fft2_out  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft2_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft2_out => NULL()

  COMPLEX(MKC),DIMENSION(:),POINTER :: fft3_in   => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: fft3_out  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft3_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft3_out => NULL()

  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft4_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft4_out => NULL()

  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft5_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft5_out => NULL()

  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft6_in  => NULL()
  COMPLEX(MKC),DIMENSION(:),POINTER :: ifft6_out => NULL()

  INTEGER                           :: fft_plan, ifft_plan

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_third => NULL()

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
  IF(fft_seq(3) .EQ. 1)THEN ! Convolve on x-pencils

    IF(bc_ptc(1) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(1)
!      ELSE
!        nfft = ncell(1) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(1)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )

    IF( ifft_vort )THEN
      ALLOCATE( ifft4_in(nfft), ifft5_in(nfft), ifft6_in(nfft) )
      ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )
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
            fft1_in(i) = vort_fft(1,i,j,k)
            fft2_in(i) = vort_fft(2,i,j,k)
            fft3_in(i) = vort_fft(3,i,j,k)
          ELSE ! zero padding
            fft1_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft2_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft3_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
          END IF
        END DO


! Fourier transform
        CALL dfftw_plan_dft_1d( fft_plan, nfft, fft1_in, fft1_out, & 
                              & FFTW_FORWARD, FFTW_ESTIMATE)
        CALL dfftw_execute_dft( fft_plan, fft1_in, fft1_out )
        CALL dfftw_execute_dft( fft_plan, fft2_in, fft2_out )
        CALL dfftw_execute_dft( fft_plan, fft3_in, fft3_out )
        CALL dfftw_destroy_plan( fft_plan )

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

! Regularise
            ifft4_in(i) = fft1_out(i) * zeta_fft
            ifft5_in(i) = fft2_out(i) * zeta_fft
            ifft6_in(i) = fft3_out(i) * zeta_fft
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
            ifft1_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kj,MKC) * fft3_out(i) &
                        &                   - CMPLX(0.0_MK,kk,MKC) * fft2_out(i) )
            ifft2_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kk,MKC) * fft1_out(i) &
                        &                   - CMPLX(0.0_MK,ki,MKC) * fft3_out(i) )
            ifft3_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC) * fft2_out(i) &
                        &                   - CMPLX(0.0_MK,kj,MKC) * fft1_out(i) )

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(i) &
                      &                   + CMPLX(0.0_MK,kj,MKC)*fft2_out(i) &
                      &                   + CMPLX(0.0_MK,kk,MKC)*fft3_out(i) )

              ifft4_in(i) = ifft4_in(i) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
              ifft5_in(i) = ifft5_in(i) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
              ifft6_in(i) = ifft6_in(i) &
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
            ifft1_in(i) = Ky * fft3_out(i) - Kz * fft2_out(i)
            ifft2_in(i) = Kz * fft1_out(i) - Kx * fft3_out(i)
            ifft3_in(i) = Kx * fft2_out(i) - Ky * fft1_out(i)

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(i) &
                       &               + CMPLX(0.0_MK,kj,MKC)*fft2_out(i) &
                       &               + CMPLX(0.0_MK,kk,MKC)*fft3_out(i) )

              ifft4_in(i) = ifft4_in(i) + Kx * div_vort
              ifft5_in(i) = ifft5_in(i) + Ky * div_vort
              ifft6_in(i) = ifft6_in(i) + Kz * div_vort

            END IF
          END IF

        END DO

! Inverse Fourier transform
        CALL dfftw_plan_dft_1d( ifft_plan, nfft, ifft1_in, ifft1_out, &
                              & FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL dfftw_execute_dft( ifft_plan, ifft1_in, ifft1_out )
        CALL dfftw_execute_dft( ifft_plan, ifft2_in, ifft2_out )
        CALL dfftw_execute_dft( ifft_plan, ifft3_in, ifft3_out )
        IF( ifft_vort )THEN
          CALL dfftw_execute_dft( ifft_plan, ifft4_in, ifft4_out )
          CALL dfftw_execute_dft( ifft_plan, ifft5_in, ifft5_out )
          CALL dfftw_execute_dft( ifft_plan, ifft6_in, ifft6_out )
        END IF
        CALL dfftw_destroy_plan( ifft_plan )

! Store fourier coefficients discard zero-padding
        DO i = 1,ncell(1)
          vel_fft(1,i,j,k) = ifft1_out(i)/REAL(nfft,MK)
          vel_fft(2,i,j,k) = ifft2_out(i)/REAL(nfft,MK)
          vel_fft(3,i,j,k) = ifft3_out(i)/REAL(nfft,MK)
          IF( ifft_vort )THEN
            vort_fft(1,i,j,k) = ifft4_out(i)/REAL(nfft,MK)
            vort_fft(2,i,j,k) = ifft5_out(i)/REAL(nfft,MK)
            vort_fft(3,i,j,k) = ifft6_out(i)/REAL(nfft,MK)
          END IF
        END DO

      END DO !j
    END DO !k

! Clean fftw
    CALL dfftw_cleanup()

  ELSEIF(fft_seq(3) .EQ. 2)THEN ! Convolve on y-pencils

    IF(bc_ptc(2) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(2)
!      ELSE
!        nfft = ncell(2) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(2)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )

    IF( ifft_vort )THEN
      ALLOCATE( ifft4_in(nfft), ifft5_in(nfft), ifft6_in(nfft) )
      ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )
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
            fft1_in(j) = vort_fft(1,i,j,k)
            fft2_in(j) = vort_fft(2,i,j,k)
            fft3_in(j) = vort_fft(3,i,j,k)
          ELSE ! zero padding
            fft1_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft2_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft3_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
          END IF
        END DO

! Fourier transform
        CALL dfftw_plan_dft_1d( fft_plan, nfft, fft1_in, fft1_out, & 
                              & FFTW_FORWARD, FFTW_ESTIMATE)
        CALL dfftw_execute_dft( fft_plan, fft1_in, fft1_out )
        CALL dfftw_execute_dft( fft_plan, fft2_in, fft2_out )
        CALL dfftw_execute_dft( fft_plan, fft3_in, fft3_out )
        CALL dfftw_destroy_plan( fft_plan )

        DO j = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
          IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
            kj = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
          ELSE
            kj = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
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

! Regularise
            ifft4_in(j) = fft1_out(j) * zeta_fft
            ifft5_in(j) = fft2_out(j) * zeta_fft
            ifft6_in(j) = fft3_out(j) * zeta_fft
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
            ifft1_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kj,MKC) * fft3_out(j) &
                        &                   - CMPLX(0.0_MK,kk,MKC) * fft2_out(j) )
            ifft2_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kk,MKC) * fft1_out(j) &
                        &                   - CMPLX(0.0_MK,ki,MKC) * fft3_out(j) )
            ifft3_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC) * fft2_out(j) &
                        &                   - CMPLX(0.0_MK,kj,MKC) * fft1_out(j) )

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(j) &
                      &                   + CMPLX(0.0_MK,kj,MKC)*fft2_out(j) &
                      &                   + CMPLX(0.0_MK,kk,MKC)*fft3_out(j) )

              ifft4_in(j) = ifft4_in(j) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
              ifft5_in(j) = ifft5_in(j) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
              ifft6_in(j) = ifft6_in(j) &
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
            ifft1_in(j) = Ky * fft3_out(j) - Kz * fft2_out(j)
            ifft2_in(j) = Kz * fft1_out(j) - Kx * fft3_out(j)
            ifft3_in(j) = Kx * fft2_out(j) - Ky * fft1_out(j)

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(j) &
                       &               + CMPLX(0.0_MK,kj,MKC)*fft2_out(j) &
                       &               + CMPLX(0.0_MK,kk,MKC)*fft3_out(j) )

              ifft4_in(j) = ifft4_in(j) + Kx * div_vort
              ifft5_in(j) = ifft5_in(j) + Ky * div_vort
              ifft6_in(j) = ifft6_in(j) + Kz * div_vort
            END IF
          END IF

        END DO

! Inverse Fourier transform
        CALL dfftw_plan_dft_1d( ifft_plan, nfft, ifft1_in, ifft1_out, &
                              & FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL dfftw_execute_dft( ifft_plan, ifft1_in, ifft1_out )
        CALL dfftw_execute_dft( ifft_plan, ifft2_in, ifft2_out )
        CALL dfftw_execute_dft( ifft_plan, ifft3_in, ifft3_out )
        IF( ifft_vort )THEN
          CALL dfftw_execute_dft( ifft_plan, ifft4_in, ifft4_out )
          CALL dfftw_execute_dft( ifft_plan, ifft5_in, ifft5_out )
          CALL dfftw_execute_dft( ifft_plan, ifft6_in, ifft6_out )
        END IF
        CALL dfftw_destroy_plan( ifft_plan )

! Store fourier coefficients discard zero-padding
        DO j = 1,ncell(2)
          vel_fft(1,i,j,k) = ifft1_out(j)/REAL(nfft,MK)
          vel_fft(2,i,j,k) = ifft2_out(j)/REAL(nfft,MK)
          vel_fft(3,i,j,k) = ifft3_out(j)/REAL(nfft,MK)
          IF( ifft_vort )THEN
            vort_fft(1,i,j,k) = ifft4_out(j)/REAL(nfft,MK)
            vort_fft(2,i,j,k) = ifft5_out(j)/REAL(nfft,MK)
            vort_fft(3,i,j,k) = ifft6_out(j)/REAL(nfft,MK)
          END IF
        END DO

      END DO !i
    END DO !k

! Clean fftw
    CALL dfftw_cleanup()

  ELSEIF(fft_seq(3) .EQ. 3)THEN ! Convolve on z-pencils

    IF(bc_ptc(3) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        nfft = 2*ncell(3)
!      ELSE
!        nfft = ncell(3) + pmlib_patch_ext
!      END IF
    ELSE
      nfft = ncell(3)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )

    IF( ifft_vort )THEN
      ALLOCATE( ifft4_in(nfft), ifft5_in(nfft), ifft6_in(nfft) )
      ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )
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
            fft1_in(k) = vort_fft(1,i,j,k)
            fft2_in(k) = vort_fft(2,i,j,k)
            fft3_in(k) = vort_fft(3,i,j,k)
          ELSE ! zero padding
            fft1_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft2_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
            fft3_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
          END IF
        END DO

! Fourier transform
        CALL dfftw_plan_dft_1d( fft_plan, nfft, fft1_in, fft1_out, & 
                              & FFTW_FORWARD, FFTW_ESTIMATE )
        CALL dfftw_execute_dft( fft_plan, fft1_in, fft1_out )
        CALL dfftw_execute_dft( fft_plan, fft2_in, fft2_out )
        CALL dfftw_execute_dft( fft_plan, fft3_in, fft3_out )
        CALL dfftw_destroy_plan( fft_plan )

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

! Regularise
            ifft4_in(k) = fft1_out(k) * zeta_fft
            ifft5_in(k) = fft2_out(k) * zeta_fft
            ifft6_in(k) = fft3_out(k) * zeta_fft
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
            ifft1_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kj,MKC) * fft3_out(k) &
                        &                   - CMPLX(0.0_MK,kk,MKC) * fft2_out(k) )
            ifft2_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,kk,MKC) * fft1_out(k) &
                        &                   - CMPLX(0.0_MK,ki,MKC) * fft3_out(k) )
            ifft3_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC) * fft2_out(k) &
                        &                   - CMPLX(0.0_MK,kj,MKC) * fft1_out(k) )

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(k) &
                      &                   + CMPLX(0.0_MK,kj,MKC)*fft2_out(k) &
                      &                   + CMPLX(0.0_MK,kk,MKC)*fft3_out(k) )

              ifft4_in(k) = ifft4_in(k) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,ki,MKC) * div_psi
              ifft5_in(k) = ifft5_in(k) &
                        & + 2.0_MK * PI * CMPLX(0.0_MK,kj,MKC) * div_psi
              ifft6_in(k) = ifft6_in(k) &
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

! convolve and do spectral curl to get the FFT coefficients of the velocity
            ifft1_in(k) = Ky * fft3_out(k) - Kz * fft2_out(k)
            ifft2_in(k) = Kz * fft1_out(k) - Kx * fft3_out(k)
            ifft3_in(k) = Kx * fft2_out(k) - Ky * fft1_out(k)

! Correct the output vorticity by the divergence of the stream function
            IF(PRESENT(reproj) .AND. reproj)THEN
              div_vort = 2.0_MK * PI * ( CMPLX(0.0_MK,ki,MKC)*fft1_out(k) &
                       &               + CMPLX(0.0_MK,kj,MKC)*fft2_out(k) &
                       &               + CMPLX(0.0_MK,kk,MKC)*fft3_out(k) )

              ifft4_in(k) = ifft4_in(k) + Kx * div_vort
              ifft5_in(k) = ifft5_in(k) + Ky * div_vort
              ifft6_in(k) = ifft6_in(k) + Kz * div_vort
            END IF
          END IF

        END DO

! Inverse Fourier transform
        CALL dfftw_plan_dft_1d( ifft_plan, nfft, ifft1_in, ifft1_out, &
                              & FFTW_BACKWARD, FFTW_ESTIMATE)
        CALL dfftw_execute_dft( ifft_plan, ifft1_in, ifft1_out )
        CALL dfftw_execute_dft( ifft_plan, ifft2_in, ifft2_out )
        CALL dfftw_execute_dft( ifft_plan, ifft3_in, ifft3_out )
        IF( ifft_vort )THEN
          CALL dfftw_execute_dft( ifft_plan, ifft4_in, ifft4_out )
          CALL dfftw_execute_dft( ifft_plan, ifft5_in, ifft5_out )
          CALL dfftw_execute_dft( ifft_plan, ifft6_in, ifft6_out )
        END IF
        CALL dfftw_destroy_plan( ifft_plan )

! Store fourier coefficients discard zero-padding
        DO k = 1,ncell(3)
          vel_fft(1,i,j,k) = ifft1_out(k)/REAL(nfft,MK)
          vel_fft(2,i,j,k) = ifft2_out(k)/REAL(nfft,MK)
          vel_fft(3,i,j,k) = ifft3_out(k)/REAL(nfft,MK)
          IF( ifft_vort )THEN
            vort_fft(1,i,j,k) = ifft4_out(k)/REAL(nfft,MK)
            vort_fft(2,i,j,k) = ifft5_out(k)/REAL(nfft,MK)
            vort_fft(3,i,j,k) = ifft6_out(k)/REAL(nfft,MK)
          END IF
        END DO

      END DO !i
    END DO !j

! Clean fftw
    CALL dfftw_cleanup()

  END IF

!---------------------------------------------------------------------------------!
! De-allocate fft arrays
!---------------------------------------------------------------------------------!
  DEALLOCATE( fft1_in, fft1_out, ifft1_in, ifft1_out )
  DEALLOCATE( fft2_in, fft2_out, ifft2_in, ifft2_out )
  DEALLOCATE( fft3_in, fft3_out, ifft3_in, ifft3_out )

  IF( ifft_vort )THEN
    DEALLOCATE( ifft4_in, ifft5_in, ifft6_in )
    DEALLOCATE( ifft4_out, ifft5_out, ifft6_out )
  END IF

!---------------------------------------------------------------------------------!
! IFFT velocity field in the first direction
!---------------------------------------------------------------------------------!
  ifft_seq(1) = patch%fft_seq(2)
  ifft_seq(2) = patch%fft_seq(1)
  ifft_seq(3) = 0
  CALL pmlib_ifft( topo_all,topo_third,vel_fft,topo_all%cuboid, &
                 & mesh%vel,ifft_seq,ierr,fft_shift = .FALSE.)

  IF( ifft_vort )THEN
    CALL pmlib_ifft( topo_all,topo_third,vort_fft,topo_all%cuboid, & 
                   & mesh%vort,ifft_seq,ierr,fft_shift = .FALSE.)
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
