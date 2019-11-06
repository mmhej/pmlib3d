!---------------------------------------------------------------------------------!
! flowcases_turbulence.f
!---------------------------------------------------------------------------------!
! Generates a stochastically synthesised turbulent vorticity field
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_turbulence(patch,topo_all,mesh,ierr,reg_vort,reproj)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_fourier
USE pmlib_mod_regularise
!USE pmlib_mod_communication

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                            :: patch
  TYPE(class_topology_all)                     :: topo_all
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr
  LOGICAL, OPTIONAL, INTENT(IN)                :: reg_vort
  LOGICAL, OPTIONAL, INTENT(IN)                :: reproj

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                 :: caller = 'flowcases_turbulence'
  REAL(MK)                                     :: PI
  INTEGER,DIMENSION(8)                         :: timedate
  INTEGER                                      :: i,j,k,l, ii, jj, kk
  REAL(MK)                                     :: px,py,pz

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  REAL(MK)                                     :: sigma, rho, scaling

  REAL(MK), DIMENSION(:), POINTER              :: aux
  REAL(MK), DIMENSION(3)                       :: u
  INTEGER, DIMENSION(:), ALLOCATABLE           :: seed
  INTEGER                                      :: nrand

  INTEGER                           :: i_ptc, j_ptc, k_ptc
  REAL(MK), DIMENSION(pmlib_ndim)   :: l_ptc, ks, dk, kmin, ksmp, ki
  INTEGER, DIMENSION(pmlib_ndim)    :: ncell_ptc
  REAL(MK)                          :: kmag, s, kolm_len

  REAL(MK)                          :: eps
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin_ptc, xmax_ptc
  INTEGER, DIMENSION(pmlib_ndim)    :: bc_ptc

  REAL(MK), DIMENSION(pmlib_ndim)   :: sum_std_vel, std_vel
  REAL(MK), DIMENSION(pmlib_ndim)   :: sum_std_vort, std_vort
  REAL(MK)                          :: max_vort

  REAL(MK)                          :: c1, c2, c3, c4, c5, c17_6

  INTEGER, DIMENSION(pmlib_ndim)    :: fft_seq, ifft_seq

  COMPLEX(MKC)                      :: G, zeta_fft
  
  REAL(MK)                          :: E
  REAL(MK), DIMENSION(pmlib_ndim,pmlib_ndim)  :: C, H

  COMPLEX(MKC)                      :: psi1_fft, psi2_fft, psi3_fft
  COMPLEX(MKC)                      :: div_psi, div_vort

  INTEGER                           :: nfft

  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER :: vort_fft, vel_fft

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

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_third => NULL()

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Construct a random vorticity field
!---------------------------------------------------------------------------------!
  length   = 0.10_MK
  sigma    = 0.10_MK
  rho      = 1.0_MK

  dx     = topo_all%cuboid(rank)%dx
  xmin   = topo_all%cuboid(rank)%xmin
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel  = 0.0_MK

  nrand = nvort * ncell(1) * ncell(2) * ncell(3)
  CALL DATE_AND_TIME(values = timedate)

  CALL RANDOM_SEED(SIZE = nrand)
  ALLOCATE(seed(nrand))
  DO i = 1,nrand
    seed(i) = timedate(8) + rank*timedate(7)
  END DO
  CALL RANDOM_SEED(PUT=seed)

  l = 0

  DO k = 1,ncell(3)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
    DO j = 1,ncell(2)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
      DO i = 1,ncell(1)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        IF( px .GE. -0.5_MK .AND. px .LE. 0.5_MK .AND. &
          & py .GE. -0.5_MK .AND. py .LE. 0.5_MK .AND. &
          & pz .GE. -0.5_MK .AND. pz .LE. 0.5_MK )THEN

          CALL RANDOM_NUMBER(u)

          mesh%vort(1,i,j,k) = 0.5_MK * SQRT(12.0_MK) * (2.0_MK*u(1) - 1.0_MK)
          mesh%vort(2,i,j,k) = 0.5_MK * SQRT(12.0_MK) * (2.0_MK*u(2) - 1.0_MK)
          mesh%vort(3,i,j,k) = 0.5_MK * SQRT(12.0_MK) * (2.0_MK*u(3) - 1.0_MK)

        END IF

      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! De-allocate auxilary vector
!---------------------------------------------------------------------------------!
  DEALLOCATE(seed)

!---------------------------------------------------------------------------------!
! Convolve the random field with covariance tensor function
!---------------------------------------------------------------------------------!
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

  ncell_ptc = (2 - bc_ptc) * (patch%ncell + 2 * ( 1 - bc_ptc ) * patch%nghost)

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

!  eps = 2.0_MK * SQRT(visc * dtime)

!  IF(rank .EQ. 0)THEN
!    IF( eps .LT. pmlib_regularisation_radius * MAXVAL(dx) )THEN
!      WRITE(*,*)'WARNING: ', eps, pmlib_regularisation_radius * MAXVAL(dx)
!    END IF
!  END IF

  eps = 3.0_MK * pmlib_regularisation_radius * MAXVAL(dx)

!---------------------------------------------------------------------------------!
! Setup coefficients for the regularisation filter
!---------------------------------------------------------------------------------!
  c1 = 1.0_MK
  c2 = 0.5_MK
  c3 = 0.125_MK
  c4 = 1.0_MK/48.0_MK
  c5 = 1.0_MK/384.0_MK

  c17_6 = 17.0_MK/6.0_MK

!---------------------------------------------------------------------------------!
! FFT the field in the third direction, convolve and IFFT
!---------------------------------------------------------------------------------!
  IF(fft_seq(3) .EQ. 1)THEN ! Convolve on x-pencils

    IF(bc_ptc(1) .EQ. 0)THEN
      nfft = 2*ncell(1)
    ELSE
      nfft = ncell(1)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )

    ALLOCATE( ifft4_in(nfft),  ifft5_in(nfft),  ifft6_in(nfft) )
    ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )

    DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
      IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
        ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
      ELSE
        ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
      END IF

      DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
        IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
          ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
        ELSE
          ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
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
        CALL fft( fft1_in, fft1_out, nfft )
        CALL fft( fft2_in, fft2_out, nfft )
        CALL fft( fft3_in, fft3_out, nfft )

        DO i = 1,nfft
! Calculate wave number in x-direction (FFT-shifted)
          IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
            ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
          ELSE
            ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
          END IF


! Calculate the energy spectrum
          kmag = 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )
          E = length**4 * kmag**4 / (1.0_MK + length**2 * kmag**2)**c17_6

! Calculate spectral covariance function
          kmag = SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )
          IF(kmag .GT. 0.0_MK)THEN
            DO jj = 1,pmlib_ndim
              DO ii = 1,pmlib_ndim

                IF(ii .EQ. jj)THEN
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( kmag**2 - ki(ii)*ki(jj) )
                ELSE
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( - ki(ii)*ki(jj) )
                END IF

              END DO
            END DO
          ELSE
            C = 0.0_MK
          END IF

! Do Cholesky decomposition
          H = 0.0_MK
          IF( i .EQ. 1 .OR. j .EQ. 1 .OR. k .EQ. 1) THEN 
            H(1,1) = SQRT( ABS( C(1,1) ) )
            H(2,2) = SQRT( ABS( C(2,2) ) )
            H(3,3) = SQRT( ABS( C(3,3) ) )
          ELSE
            H(1,1) = SQRT( MAX(0.0_MK, C(1,1) ) )
            H(2,1) = C(2,1) / H(1,1)
            H(3,1) = C(3,1) / H(1,1)
            H(2,2) = SQRT( MAX(0.0_MK, C(2,2) - H(2,1) * H(2,1) ) )
            H(3,2) = ( C(3,2) - H(3,1) * H(2,1) ) / H(2,2)
            H(3,3) = SQRT( MAX(0.0_MK, C(3,3) - ( H(3,1)*H(3,1) + H(3,2)*H(3,2) )))
          END IF

! Convolve random vorticity field with covariance function
          fft1_out(i) = H(1,1) * fft1_out(i) &
                    & + H(1,2) * fft2_out(i) &
                    & + H(1,3) * fft3_out(i)

          fft2_out(i) = H(2,1) * fft1_out(i) &
                    & + H(2,2) * fft2_out(i) &
                    & + H(2,3) * fft3_out(i)

          fft3_out(i) = H(3,1) * fft1_out(i) &
                    & + H(3,2) * fft2_out(i) &
                    & + H(3,3) * fft3_out(i)

! Regularise the vorticity field
          IF( PRESENT(reg_vort) .AND. reg_vort )THEN
            s = eps * 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )

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
            ifft1_in(i) = fft1_out(i) * zeta_fft
            ifft2_in(i) = fft2_out(i) * zeta_fft
            ifft3_in(i) = fft3_out(i) * zeta_fft
          ELSE
            ifft1_in(i) = fft1_out(i)
            ifft2_in(i) = fft2_out(i)
            ifft3_in(i) = fft3_out(i)

            zeta_fft = 1.0_MK
          END IF

! Get the Poisson integration kernel
          IF( patch%bound_cond(1) .NE. 1 .AND. &
            & patch%bound_cond(2) .NE. 1 .AND. &
            & patch%bound_cond(3) .NE. 1 )THEN

            G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
          ELSE
            kmag = ki(1)**2 + ki(2)**2 + ki(3)**2
            IF(kmag .GT. 0.0_MK)THEN
              G = zeta_fft/(4.0_MK * PI**2 * kmag )
            ELSE
              G = 0.0_MK
            END IF
          END IF

! Correct the vorticity by the divergence of the vector potential
          IF(PRESENT(reproj) .AND. reproj)THEN

            div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC)*fft1_out(i) &
                    &                   + CMPLX(0.0_MK,ki(2),MKC)*fft2_out(i) &
                    &                   + CMPLX(0.0_MK,ki(3),MKC)*fft3_out(i) )

            ifft1_in(i) = ifft1_in(i) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(1),MKC) * div_psi
            ifft2_in(i) = ifft2_in(i) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(2),MKC) * div_psi
            ifft3_in(i) = ifft3_in(i) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(3),MKC) * div_psi
          END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
          ifft4_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(2),MKC) * fft3_out(i) &
                      &                   - CMPLX(0.0_MK,ki(3),MKC) * fft2_out(i) )
          ifft5_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(3),MKC) * fft1_out(i) &
                      &                   - CMPLX(0.0_MK,ki(1),MKC) * fft3_out(i) )
          ifft6_in(i) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC) * fft2_out(i) &
                      &                   - CMPLX(0.0_MK,ki(2),MKC) * fft1_out(i) )

        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        CALL ifft( ifft2_in, ifft2_out, nfft )
        CALL ifft( ifft3_in, ifft3_out, nfft )

        CALL ifft( ifft4_in, ifft4_out, nfft )
        CALL ifft( ifft5_in, ifft5_out, nfft )
        CALL ifft( ifft6_in, ifft6_out, nfft )

! Store fourier coefficients discard zero-padding
        DO i = 1,ncell(1)
          vort_fft(1,i,j,k) = ifft1_out(i)/REAL(nfft,MK)
          vort_fft(2,i,j,k) = ifft2_out(i)/REAL(nfft,MK)
          vort_fft(3,i,j,k) = ifft3_out(i)/REAL(nfft,MK)

          vel_fft(1,i,j,k) = ifft4_out(i)/REAL(nfft,MK)
          vel_fft(2,i,j,k) = ifft5_out(i)/REAL(nfft,MK)
          vel_fft(3,i,j,k) = ifft6_out(i)/REAL(nfft,MK)
        END DO

      END DO !j
    END DO !k

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

    ALLOCATE( ifft4_in(nfft),  ifft5_in(nfft),  ifft6_in(nfft) )
    ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )

    DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
      IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
        ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
      ELSE
        ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
      END IF

      DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
        IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
        ELSE
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
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
        CALL fft( fft1_in, fft1_out, nfft )
        CALL fft( fft2_in, fft2_out, nfft )
        CALL fft( fft3_in, fft3_out, nfft )

        DO j = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
          IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
            ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
          ELSE
            ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
          END IF

! Calculate the energy spectrum
          kmag = 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )
          E = length**4 * kmag**4 / (1.0_MK + length**2 * kmag**2)**c17_6

! Calculate spectral covariance function
          kmag = SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )
          IF(kmag .GT. 0.0_MK)THEN
            DO jj = 1,pmlib_ndim
              DO ii = 1,pmlib_ndim

                IF(ii .EQ. jj)THEN
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( kmag**2 - ki(ii)*ki(jj) )
                ELSE
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( - ki(ii)*ki(jj) )
                END IF

              END DO
            END DO
          ELSE
            C = 0.0_MK
          END IF

! Do Cholesky decomposition
          H = 0.0_MK
          IF( i .EQ. 1 .OR. j .EQ. 1 .OR. k .EQ. 1) THEN 
            H(1,1) = SQRT( ABS( C(1,1) ) )
            H(2,2) = SQRT( ABS( C(2,2) ) )
            H(3,3) = SQRT( ABS( C(3,3) ) )
          ELSE
            H(1,1) = SQRT( MAX(0.0_MK, C(1,1) ) )
            H(2,1) = C(2,1) / H(1,1)
            H(3,1) = C(3,1) / H(1,1)
            H(2,2) = SQRT( MAX(0.0_MK, C(2,2) - H(2,1) * H(2,1) ) )
            H(3,2) = ( C(3,2) - H(3,1) * H(2,1) ) / H(2,2)
            H(3,3) = SQRT( MAX(0.0_MK, C(3,3) - ( H(3,1)*H(3,1) + H(3,2)*H(3,2) )))
          END IF

! Convolve random vorticity field with covariance function
          fft1_out(j) = H(1,1) * fft1_out(j) &
                    & + H(1,2) * fft2_out(j) &
                    & + H(1,3) * fft3_out(j)

          fft2_out(j) = H(2,1) * fft1_out(j) &
                    & + H(2,2) * fft2_out(j) &
                    & + H(2,3) * fft3_out(j)

          fft3_out(j) = H(3,1) * fft1_out(j) &
                    & + H(3,2) * fft2_out(j) &
                    & + H(3,3) * fft3_out(j)

! Regularise the output vorticity field
          IF( PRESENT(reg_vort) .AND. reg_vort )THEN
            s = eps * 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )

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
            ifft1_in(j) = fft1_out(j) * zeta_fft
            ifft2_in(j) = fft2_out(j) * zeta_fft
            ifft3_in(j) = fft3_out(j) * zeta_fft
          ELSE
            ifft1_in(j) = fft1_out(j)
            ifft2_in(j) = fft2_out(j)
            ifft3_in(j) = fft3_out(j)

            zeta_fft = 1.0_MK
          END IF

! Get the Poisson integration kernel
          IF( patch%bound_cond(1) .NE. 1 .AND. &
            & patch%bound_cond(2) .NE. 1 .AND. &
            & patch%bound_cond(3) .NE. 1 )THEN

            G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
          ELSE
            kmag = ki(1)**2 + ki(2)**2 + ki(3)**2
            IF(kmag .GT. 0.0_MK)THEN
              G = zeta_fft/(4.0_MK * PI**2 * kmag )
            ELSE
              G = 0.0_MK
            END IF
          END IF

! Correct the output vorticity by the divergence of the stream function
          IF(PRESENT(reproj) .AND. reproj)THEN
            div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC)*fft1_out(j) &
                    &                   + CMPLX(0.0_MK,ki(2),MKC)*fft2_out(j) &
                    &                   + CMPLX(0.0_MK,ki(3),MKC)*fft3_out(j) )

            ifft1_in(j) = ifft1_in(j) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(1),MKC) * div_psi
            ifft2_in(j) = ifft2_in(j) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(2),MKC) * div_psi
            ifft3_in(j) = ifft3_in(j) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(3),MKC) * div_psi
          END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
          ifft4_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(2),MKC) * fft3_out(j) &
                      &                   - CMPLX(0.0_MK,ki(3),MKC) * fft2_out(j) )
          ifft5_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(3),MKC) * fft1_out(j) &
                      &                   - CMPLX(0.0_MK,ki(1),MKC) * fft3_out(j) )
          ifft6_in(j) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC) * fft2_out(j) &
                      &                   - CMPLX(0.0_MK,ki(2),MKC) * fft1_out(j) )

        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        CALL ifft( ifft2_in, ifft2_out, nfft )
        CALL ifft( ifft3_in, ifft3_out, nfft )

        CALL ifft( ifft4_in, ifft4_out, nfft )
        CALL ifft( ifft5_in, ifft5_out, nfft )
        CALL ifft( ifft6_in, ifft6_out, nfft )

! Store fourier coefficients discard zero-padding
        DO j = 1,ncell(2)
          vort_fft(1,i,j,k) = ifft1_out(j)/REAL(nfft,MK)
          vort_fft(2,i,j,k) = ifft2_out(j)/REAL(nfft,MK)
          vort_fft(3,i,j,k) = ifft3_out(j)/REAL(nfft,MK)

          vel_fft(1,i,j,k) = ifft4_out(j)/REAL(nfft,MK)
          vel_fft(2,i,j,k) = ifft5_out(j)/REAL(nfft,MK)
          vel_fft(3,i,j,k) = ifft6_out(j)/REAL(nfft,MK)
        END DO

      END DO !i
    END DO !k


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

    ALLOCATE( ifft4_in(nfft),  ifft5_in(nfft),  ifft6_in(nfft) )
    ALLOCATE( ifft4_out(nfft), ifft5_out(nfft), ifft6_out(nfft) )

    DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
      IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
        ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) + ksmp(2)
      ELSE
        ki(2) = kmin(2) + REAL(j - 1,MK)*dk(2) - ksmp(2)
      END IF

      DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
        IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
        ELSE
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
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
        CALL fft( fft1_in, fft1_out, nfft )
        CALL fft( fft2_in, fft2_out, nfft )
        CALL fft( fft3_in, fft3_out, nfft )

        DO k = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
          IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
            ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) + ksmp(3)
          ELSE
            ki(3) = kmin(3) + REAL(k - 1,MK)*dk(3) - ksmp(3)
          END IF

! Calculate the energy spectrum
          kmag = 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )
          E = length**4 * kmag**4 / (1.0_MK + length**2 * kmag**2)**c17_6

! Calculate spectral covariance function
          kmag = SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )

          IF(kmag .GT. 0.0_MK)THEN

            DO jj = 1,pmlib_ndim
              DO ii = 1,pmlib_ndim

                IF(ii .EQ. jj)THEN
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( kmag**2 - ki(ii)*ki(jj) )
                ELSE
                  C(ii,jj) = PI * E / (rho * kmag**2) * ( - ki(ii)*ki(jj) )
                END IF

              END DO
            END DO
          ELSE
            C = 0.0_MK
          END IF

! Do Cholesky decomposition
          H = 0.0_MK
          IF( i .EQ. 1 .OR. j .EQ. 1 .OR. k .EQ. 1) THEN 
            H(1,1) = SQRT( ABS( C(1,1) ) )
            H(2,2) = SQRT( ABS( C(2,2) ) )
            H(3,3) = SQRT( ABS( C(3,3) ) )
          ELSE
            H(1,1) = SQRT( MAX(0.0_MK, C(1,1) ) )
            H(2,1) = C(2,1) / H(1,1)
            H(3,1) = C(3,1) / H(1,1)
            H(2,2) = SQRT( MAX(0.0_MK, C(2,2) - H(2,1) * H(2,1) ) )
            H(3,2) = ( C(3,2) - H(3,1) * H(2,1) ) / H(2,2)
            H(3,3) = SQRT( MAX(0.0_MK, C(3,3) - ( H(3,1)*H(3,1) + H(3,2)*H(3,2) )))
          END IF

! Convolve random vorticity field with covariance function
          fft1_out(k) = H(1,1) * fft1_out(k) &
                    & + H(1,2) * fft2_out(k) &
                    & + H(1,3) * fft3_out(k)

          fft2_out(k) = H(2,1) * fft1_out(k) &
                    & + H(2,2) * fft2_out(k) &
                    & + H(2,3) * fft3_out(k)

          fft3_out(k) = H(3,1) * fft1_out(k) &
                    & + H(3,2) * fft2_out(k) &
                    & + H(3,3) * fft3_out(k)

! Regularise the output vorticity field
          IF( PRESENT(reg_vort) .AND. reg_vort )THEN
            s = eps * 2.0_MK * PI * SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )

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
            ifft1_in(k) = fft1_out(k) * zeta_fft
            ifft2_in(k) = fft2_out(k) * zeta_fft
            ifft3_in(k) = fft3_out(k) * zeta_fft
          ELSE
            ifft1_in(k) = fft1_out(k)
            ifft2_in(k) = fft2_out(k)
            ifft3_in(k) = fft3_out(k)

            zeta_fft = 1.0_MK
          END IF

! Get the Poisson integration kernel
          IF( patch%bound_cond(1) .NE. 1 .AND. &
            & patch%bound_cond(2) .NE. 1 .AND. &
            & patch%bound_cond(3) .NE. 1 )THEN

            G = mesh%pkernel_fft(1,i,j,k) *dx(1)*dx(2)*dx(3)
          ELSE
            kmag = ki(1)**2 + ki(2)**2 + ki(3)**2
            IF(kmag .GT. 0.0_MK)THEN
              G = zeta_fft/(4.0_MK * PI**2 * kmag )
            ELSE
              G = 0.0_MK
            END IF
          END IF

! Correct the output vorticity by the divergence of the stream function
          IF(PRESENT(reproj) .AND. reproj)THEN
            div_psi = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC)*fft1_out(k) &
                    &                   + CMPLX(0.0_MK,ki(2),MKC)*fft2_out(k) &
                    &                   + CMPLX(0.0_MK,ki(3),MKC)*fft3_out(k) )

            ifft1_in(k) = ifft1_in(k) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(1),MKC) * div_psi
            ifft2_in(k) = ifft2_in(k) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(2),MKC) * div_psi
            ifft3_in(k) = ifft3_in(k) &
                      & + 2.0_MK * PI * CMPLX(0.0_MK,ki(3),MKC) * div_psi
          END IF

! Convolve and do spectral curl to get the FFT coefficients of the velocity
          ifft4_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(2),MKC) * fft3_out(k) &
                      &                   - CMPLX(0.0_MK,ki(3),MKC) * fft2_out(k) )
          ifft5_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(3),MKC) * fft1_out(k) &
                      &                   - CMPLX(0.0_MK,ki(1),MKC) * fft3_out(k) )
          ifft6_in(k) = G * 2.0_MK * PI * ( CMPLX(0.0_MK,ki(1),MKC) * fft2_out(k) &
                      &                   - CMPLX(0.0_MK,ki(2),MKC) * fft1_out(k) )

        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        CALL ifft( ifft2_in, ifft2_out, nfft )
        CALL ifft( ifft3_in, ifft3_out, nfft )

        CALL ifft( ifft4_in, ifft4_out, nfft )
        CALL ifft( ifft5_in, ifft5_out, nfft )
        CALL ifft( ifft6_in, ifft6_out, nfft )

! Store fourier coefficients discard zero-padding
        DO k = 1,ncell(3)
          vort_fft(1,i,j,k) = ifft1_out(k)/REAL(nfft,MK)
          vort_fft(2,i,j,k) = ifft2_out(k)/REAL(nfft,MK)
          vort_fft(3,i,j,k) = ifft3_out(k)/REAL(nfft,MK)

          vel_fft(1,i,j,k)  = ifft4_out(k)/REAL(nfft,MK)
          vel_fft(2,i,j,k)  = ifft5_out(k)/REAL(nfft,MK)
          vel_fft(3,i,j,k)  = ifft6_out(k)/REAL(nfft,MK)
        END DO

      END DO !i
    END DO !j

  END IF

!---------------------------------------------------------------------------------!
! IFFT vorticity field in the first two directions
!---------------------------------------------------------------------------------!
  ifft_seq(1) = patch%fft_seq(2)
  ifft_seq(2) = patch%fft_seq(1)
  ifft_seq(3) = 0

  CALL pmlib_ifft( topo_all,topo_third,vort_fft,topo_all%cuboid,mesh%vort, &
                 & ifft_seq,ierr,fft_shift = .FALSE.)

  CALL pmlib_ifft( topo_all,topo_third,vel_fft,topo_all%cuboid,mesh%vel, &
                 & ifft_seq,ierr,fft_shift = .FALSE.)

!---------------------------------------------------------------------------------!
! Rescale turbulence field
!---------------------------------------------------------------------------------!
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

  max_vort = 0.0_MK
  sum_std_vel  = 0.0_MK
  DO k = 1,ncell(3)
    DO j = 1,ncell(2)
      DO i = 1,ncell(1)

        sum_std_vel(1) = sum_std_vel(1) + mesh%vel(1,i,j,k)**2
        sum_std_vel(2) = sum_std_vel(2) + mesh%vel(2,i,j,k)**2
        sum_std_vel(3) = sum_std_vel(3) + mesh%vel(3,i,j,k)**2

      END DO !i
    END DO !j
  END DO !k

  CALL MPI_ALLREDUCE( sum_std_vel,std_vel,pmlib_ndim,mpi_prec_real,MPI_SUM,mpi_comm,ierr )

  std_vel = SQRT( std_vel / (patch%ncell(1)*patch%ncell(2)*patch%ncell(3)) )

  scaling = ( sigma**3/(std_vel(1)*std_vel(2)*std_vel(3)) )**(1.0_MK/3.0_MK)

  velocity = (std_vel(1)*std_vel(2)*std_vel(3))**(1.0_MK/3.0_MK)

!  IF(rank .eq. 0)THEN
!    write(*,*)'Scaling: ',scaling, 'std. vel: ',std
!  END IF

!  OPEN(20,FILE = 'init_vel.dat')

!  WRITE(20,'(3I20)') ncell(1),ncell(2),ncell(3)
!  WRITE(20,'(3E20.12)') dx(1), dx(2), dx(2)

  sum_std_vel  = 0.0_MK
  sum_std_vort = 0.0_MK
  std_vel   = 0.0_MK
  std_vort  = 0.0_MK
  DO k = 1-nghost(5),ncell(3)+nghost(6)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        mesh%vort(1,i,j,k) = scaling * mesh%vort(1,i,j,k)
        mesh%vort(2,i,j,k) = scaling * mesh%vort(2,i,j,k)
        mesh%vort(3,i,j,k) = scaling * mesh%vort(3,i,j,k)

        mesh%vel(1,i,j,k)  = scaling * mesh%vel(1,i,j,k)
        mesh%vel(2,i,j,k)  = scaling * mesh%vel(2,i,j,k)
        mesh%vel(3,i,j,k)  = scaling * mesh%vel(3,i,j,k)

        IF( SQRT( mesh%vort(1,i,j,k)**2 &
              & + mesh%vort(2,i,j,k)**2 &
              & + mesh%vort(3,i,j,k)**2) .GT. max_vort ) THEN

          max_vort = SQRT( mesh%vort(1,i,j,k)**2 &
                       & + mesh%vort(2,i,j,k)**2 &
                       & + mesh%vort(3,i,j,k)**2 )
        END IF

        IF( i .GE. 1 .AND. i .LE. ncell(1) .AND. & 
          & j .GE. 1 .AND. j .LE. ncell(2) .AND. & 
          & k .GE. 1 .AND. k .LE. ncell(3) )THEN

!          WRITE(20,'(3E20.12)') mesh%vel(1,i,j,k),mesh%vel(2,i,j,k),mesh%vel(3,i,j,k)

          sum_std_vel(1) = sum_std_vel(1) + mesh%vel(1,i,j,k)**2
          sum_std_vel(2) = sum_std_vel(2) + mesh%vel(2,i,j,k)**2
          sum_std_vel(3) = sum_std_vel(3) + mesh%vel(3,i,j,k)**2

          sum_std_vort(1) = sum_std_vort(1) + mesh%vort(1,i,j,k)**2
          sum_std_vort(2) = sum_std_vort(2) + mesh%vort(2,i,j,k)**2
          sum_std_vort(3) = sum_std_vort(3) + mesh%vort(3,i,j,k)**2
        END IF

      END DO !i
    END DO !j
  END DO !k

!  CLOSE(20)

  CALL MPI_ALLREDUCE( max_vort,vort_max,1,mpi_prec_real,MPI_MAX,comm,ierr  )

  CALL MPI_ALLREDUCE(  sum_std_vel, std_vel,pmlib_ndim,mpi_prec_real,MPI_SUM,mpi_comm,ierr )
  CALL MPI_ALLREDUCE( sum_std_vort,std_vort,pmlib_ndim,mpi_prec_real,MPI_SUM,mpi_comm,ierr )

  std_vel  = SQRT(  std_vel / (patch%ncell(1)*patch%ncell(2)*patch%ncell(3)) )
  std_vort = SQRT( std_vort / (patch%ncell(1)*patch%ncell(2)*patch%ncell(3)) )

!---------------------------------------------------------------------------------!
! Regularise to Kolmogorov length
!---------------------------------------------------------------------------------!
  visc = (0.5_MK*eps)**2 * SQRT(std_vort(1)**2 + std_vort(2)**2 + std_vort(3)**2)

  velocity = sigma

  kolm_len = ( visc**2 / (std_vort(1)**2 + std_vort(2)**2 + std_vort(3)**2) )**0.25_MK

!  IF( kolm_len .GT. eps )THEN
!    CALL pmlib_regularise(patch,topo_all,mesh%vort,pmlib_poisson_order, kolm_len ,ierr)
!  ELSE
!    IF(rank .eq. 0)THEN
!      WRITE(*,*)'WARNING: Kolmogorov scale is not resolved'
!    END IF
!  END IF


!---------------------------------------------------------------------------------!
! Output
!---------------------------------------------------------------------------------!
  IF(rank .eq. 0)THEN
    WRITE(*,'(A,E12.3)')' Reynolds number:         ', length*sigma/visc
    WRITE(*,'(A,E12.3)')' Smallest resolved scale: ', eps
    WRITE(*,'(A,2E12.3)')' Kolmogorov scale:        ', kolm_len, ( visc**2/vort_max**2 )**0.25_MK
    WRITE(*,*)''
    WRITE(*,'(A,E12.3)') ' Scaling:        ', scaling
    WRITE(*,'(A,3E12.3)')' std. velocity:  ', std_vel
    WRITE(*,'(A,3E12.3)')' std. vorticity: ', std_vort
  END IF

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
  DEALLOCATE(topo_third,stat=ierr)

  DEALLOCATE( fft1_in, fft1_out, ifft1_in, ifft1_out )
  DEALLOCATE( fft2_in, fft2_out, ifft2_in, ifft2_out )
  DEALLOCATE( fft3_in, fft3_out, ifft3_in, ifft3_out )

  DEALLOCATE( ifft4_in,  ifft5_in,  ifft6_in )
  DEALLOCATE( ifft4_out, ifft5_out, ifft6_out )

  IF( ASSOCIATED(vort_fft) ) THEN
    DEALLOCATE(vort_fft)
  END IF
  IF( ASSOCIATED(vel_fft) ) THEN
    DEALLOCATE(vel_fft)
  END IF  

!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .FALSE.

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE flowcases_turbulence

