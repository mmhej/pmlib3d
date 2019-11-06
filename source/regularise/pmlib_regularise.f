!---------------------------------------------------------------------------------!
! pmlib_regularise.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_regularise(patch,topo_all,mesh,reg_order,filter_width,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_communication
USE pmlib_mod_fourier

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                                   :: patch
  TYPE(class_topology_all)                            :: topo_all
  REAL(MK), DIMENSION(:,:,:,:),POINTER,INTENT(INOUT)  :: mesh
  INTEGER, INTENT(IN)                                 :: reg_order
  REAL(MK), INTENT(IN)                                :: filter_width
  INTEGER, INTENT(OUT)                                :: ierr


!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=16)                 :: caller = 'pmlib_regularise'
  REAL(MK)                          :: PI
  INTEGER                           :: i,j,k,l

  REAL(MK)                          :: fact

  REAL(MK)                          :: px, py, pz
  REAL(MK)                          :: rho

  INTEGER                           :: nvar
  INTEGER, DIMENSION(pmlib_ndim)    :: ncell
  INTEGER, DIMENSION(2*pmlib_ndim)  :: nghost
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin
  REAL(MK), DIMENSION(pmlib_ndim)   :: dx 

  INTEGER                           :: i_ptc, j_ptc, k_ptc
  REAL(MK), DIMENSION(pmlib_ndim)   :: l_ptc, ks, dk, kmin
  INTEGER, DIMENSION(pmlib_ndim)    :: ncell_ptc
  REAL(MK)                          :: ki, kj, kk, s
  REAL(MK)                          :: eps
  REAL(MK), DIMENSION(pmlib_ndim)   :: xmin_ptc, xmax_ptc
  INTEGER, DIMENSION(pmlib_ndim)    :: bc_ptc

  REAL(MK)                          :: c1, c2, c3, c4, c5

  INTEGER, DIMENSION(pmlib_ndim)    :: fft_seq, ifft_seq

  COMPLEX(MKC)                      :: zeta_fft

  INTEGER                           :: nfft

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

  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER   :: mesh_fft => NULL()

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_third => NULL()

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

  nvar = SIZE(mesh,1)
  IF(nvar .GT. 3)THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, 'Cant regularise variables of ndim > 3.')
    GOTO 9999
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
! Allocate Fourier coefficients for mesh values
!---------------------------------------------------------------------------------!
  ncell   = topo_third(mpi_rank)%ncell
  ALLOCATE( mesh_fft(pmlib_ndim,ncell(1),ncell(2),ncell(3)) )

!---------------------------------------------------------------------------------!
! FFT the field in the first two directions and map to the third direction
!---------------------------------------------------------------------------------!
  fft_seq(3) = 0
  CALL pmlib_fft( topo_all,topo_all%cuboid,mesh,topo_third,mesh_fft, &
                & fft_seq,ierr,fft_shift = .FALSE. )

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
  dk(1) = 2.0_MK*PI/l_ptc(1)
  dk(2) = 2.0_MK*PI/l_ptc(2)
  dk(3) = 2.0_MK*PI/l_ptc(3)

! patch index
  i_ptc = NINT((xmin(1) - xmin_ptc(1))/dx(1))
  j_ptc = NINT((xmin(2) - xmin_ptc(2))/dx(2))
  k_ptc = NINT((xmin(3) - xmin_ptc(3))/dx(3))

  kmin(1) = - PI/dx(1) + REAL(i_ptc,MK) * dk(1)
  kmin(2) = - PI/dx(2) + REAL(j_ptc,MK) * dk(2)
  kmin(3) = - PI/dx(3) + REAL(k_ptc,MK) * dk(3)

  eps = filter_width

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
      nfft = 2*ncell(1)
    ELSE
      nfft = ncell(1)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    IF(nvar .GT. 1)THEN
      ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    END IF
    IF(nvar .GT. 2)THEN
      ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )
    END IF

    DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
      IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
        kk = kmin(3) + REAL(k - 1,MK)*dk(3) + PI/dx(3)
      ELSE
        kk = kmin(3) + REAL(k - 1,MK)*dk(3) - PI/dx(3)
      END IF

      DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
        IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
          kj = kmin(2) + REAL(j - 1,MK)*dk(2) + PI/dx(2)
        ELSE
          kj = kmin(2) + REAL(j - 1,MK)*dk(2) - PI/dx(2)
        END IF

! Store fft pencil
        DO i = 1, nfft
          IF( i .LE. ncell(1) )THEN
            fft1_in(i) = mesh_fft(1,i,j,k)
            IF(nvar .GT. 1)THEN
              fft2_in(i) = mesh_fft(2,i,j,k)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(i) = mesh_fft(3,i,j,k)
            END IF
          ELSE ! zero padding
            fft1_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
            IF(nvar .GT. 1)THEN
              fft2_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(i) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
          END IF
        END DO

! Fourier transform
        CALL fft( fft1_in, fft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL fft( fft2_in, fft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL fft( fft3_in, fft3_out, nfft )
        END IF

        DO i = 1,nfft
! Calculate wave number in x-direction (FFT-shifted)
          IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
            ki = kmin(1) + REAL(i - 1,MK)*dk(1) + PI/dx(1)
          ELSE
            ki = kmin(1) + REAL(i - 1,MK)*dk(1) - PI/dx(1)
          END IF 

! Get basis function
          s = eps*SQRT( ki**2 + kj**2 + kk**2 )

          IF(reg_order .EQ. 0)THEN             
            zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 2)THEN
            zeta_fft = CMPLX( c1  &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 4)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 6)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 8)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 10)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                     & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          END IF

! Regularise
          ifft1_in(i) = fft1_out(i) * zeta_fft
          IF(nvar .GT. 1)THEN
            ifft2_in(i) = fft2_out(i) * zeta_fft
          END IF
          IF(nvar .GT. 2)THEN
            ifft3_in(i) = fft3_out(i) * zeta_fft
          END IF
        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL ifft( ifft2_in, ifft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL ifft( ifft3_in, ifft3_out, nfft )
        END IF

! Store fourier coefficients discard zero-padding
        DO i = 1,ncell(1)
          mesh_fft(1,i,j,k) = ifft1_out(i)/REAL(nfft,MK)
          IF(nvar .GT. 1)THEN
            mesh_fft(2,i,j,k) = ifft2_out(i)/REAL(nfft,MK)
          END IF
          IF(nvar .GT. 2)THEN
            mesh_fft(3,i,j,k) = ifft3_out(i)/REAL(nfft,MK)
          END IF 
        END DO

      END DO !j
    END DO !k

  ELSEIF(fft_seq(3) .EQ. 2)THEN ! Convolve on y-pencils

    IF(bc_ptc(2) .EQ. 0)THEN
      nfft = 2*ncell(2)
    ELSE
      nfft = ncell(2)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    IF(nvar .GT. 1)THEN
      ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    END IF
    IF(nvar .GT. 2)THEN
      ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )
    END IF

    DO k = 1, ncell(3)
! Calculate wave number (FFT-shifted)
      IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
        kk = kmin(3) + REAL(k - 1,MK)*dk(3) + PI/dx(3)
      ELSE
        kk = kmin(3) + REAL(k - 1,MK)*dk(3) - PI/dx(3)
      END IF

      DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
        IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
          ki = kmin(1) + REAL(i - 1,MK)*dk(1) + PI/dx(1)
        ELSE
          ki = kmin(1) + REAL(i - 1,MK)*dk(1) - PI/dx(1)
        END IF 

! Store fft pencil
        DO j = 1, nfft
          IF( j .LE. ncell(2) )THEN
            fft1_in(j) = mesh_fft(1,i,j,k)
            IF(nvar .GT. 1)THEN
              fft2_in(j) = mesh_fft(2,i,j,k)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(j) = mesh_fft(3,i,j,k)
            END IF
          ELSE ! zero padding
            fft1_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
            IF(nvar .GT. 1)THEN
              fft2_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(j) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
          END IF
        END DO

! Fourier transform
        CALL fft( fft1_in, fft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL fft( fft2_in, fft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL fft( fft3_in, fft3_out, nfft )
        END IF

        DO j = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
          IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
            kj = kmin(2) + REAL(j - 1,MK)*dk(2) + PI/dx(2)
          ELSE
            kj = kmin(2) + REAL(j - 1,MK)*dk(2) - PI/dx(2)
          END IF

! Get basis function
          s = eps*SQRT( ki**2 + kj**2 + kk**2 )

          IF(reg_order .EQ. 0)THEN             
            zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 2)THEN
            zeta_fft = CMPLX( c1  &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 4)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 6)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 8)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 10)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                     & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          END IF

! Regularise
          ifft1_in(j) = fft1_out(j) * zeta_fft
          IF(nvar .GT. 1)THEN
            ifft2_in(j) = fft2_out(j) * zeta_fft
          END IF
          IF(nvar .GT. 2)THEN
            ifft3_in(j) = fft3_out(j) * zeta_fft
          END IF
        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL ifft( ifft2_in, ifft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL ifft( ifft3_in, ifft3_out, nfft )
        END IF

! Store fourier coefficients discard zero-padding
        DO j = 1,ncell(2)
          mesh_fft(1,i,j,k) = ifft1_out(j)/REAL(nfft,MK)
          IF(nvar .GT. 1)THEN
            mesh_fft(2,i,j,k) = ifft2_out(j)/REAL(nfft,MK)
          END IF
          IF(nvar .GT. 2)THEN
            mesh_fft(3,i,j,k) = ifft3_out(j)/REAL(nfft,MK)
          END IF
        END DO

      END DO !i
    END DO !k

  ELSEIF(fft_seq(3) .EQ. 3)THEN ! Convolve on z-pencils

    IF(bc_ptc(3) .EQ. 0)THEN
      nfft = 2*ncell(3)
    ELSE
      nfft = ncell(3)
    END IF

! Allocate fft pencils
    ALLOCATE( fft1_in(nfft), fft1_out(nfft), ifft1_in(nfft), ifft1_out(nfft) )
    IF(nvar .GT. 1)THEN
      ALLOCATE( fft2_in(nfft), fft2_out(nfft), ifft2_in(nfft), ifft2_out(nfft) )
    END IF
    IF(nvar .GT. 2)THEN
      ALLOCATE( fft3_in(nfft), fft3_out(nfft), ifft3_in(nfft), ifft3_out(nfft) )
    END IF

    DO j = 1, ncell(2)
! Calculate wave number (FFT-shifted)
      IF( j + j_ptc .LE. ncell_ptc(2)/2 )THEN
        kj = kmin(2) + REAL(j - 1,MK)*dk(2) + PI/dx(2)
      ELSE
        kj = kmin(2) + REAL(j - 1,MK)*dk(2) - PI/dx(2)
      END IF

      DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
        IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
          ki = kmin(1) + REAL(i - 1,MK)*dk(1) + PI/dx(1)
        ELSE
          ki = kmin(1) + REAL(i - 1,MK)*dk(1) - PI/dx(1)
        END IF 

! Store fft pencil
        DO k = 1, nfft
          IF( k .LE. ncell(3) )THEN
            fft1_in(k) = mesh_fft(1,i,j,k)
            IF(nvar .GT. 1)THEN
              fft2_in(k) = mesh_fft(2,i,j,k)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(k) = mesh_fft(3,i,j,k)
            END IF
          ELSE ! zero padding
            fft1_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
            IF(nvar .GT. 1)THEN
              fft2_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
            IF(nvar .GT. 2)THEN
              fft3_in(k) = CMPLX(0.0_MK,0.0_MK,MKC)
            END IF
          END IF
        END DO

! Fourier transform
        CALL fft( fft1_in, fft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL fft( fft2_in, fft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL fft( fft3_in, fft3_out, nfft )
        END IF

        DO k = 1, nfft
! Calculate wave number in y-direction (FFT-shifted)
          IF( k + k_ptc .LE. ncell_ptc(3)/2 )THEN
            kk = kmin(3) + REAL(k - 1,MK)*dk(3) + PI/dx(3)
          ELSE
            kk = kmin(3) + REAL(k - 1,MK)*dk(3) - PI/dx(3)
          END IF

! Get basis function
          s = eps*SQRT( ki**2 + kj**2 + kk**2 )

          IF(reg_order .EQ. 0)THEN             
            zeta_fft = CMPLX(1.0_MK,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 2)THEN
            zeta_fft = CMPLX( c1  &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 4)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 6)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 8)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 ) &
                     & * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          ELSEIF(reg_order .EQ. 10)THEN
            zeta_fft = CMPLX(( c1 + c2*s**2 + c3*s**4 + c4*s**6 &
                     & + c5*s**8 ) * EXP(-0.5_MK*s**2) ,0.0_MK,MKC)
          END IF

! Regularise
          ifft1_in(k) = fft1_out(k) * zeta_fft
          IF(nvar .GT. 1)THEN
            ifft2_in(k) = fft2_out(k) * zeta_fft
          END IF
          IF(nvar .GT. 2)THEN
            ifft3_in(k) = fft3_out(k) * zeta_fft
          END IF

        END DO

! Inverse Fourier transform
        CALL ifft( ifft1_in, ifft1_out, nfft )
        IF(nvar .GT. 1)THEN
          CALL ifft( ifft2_in, ifft2_out, nfft )
        END IF
        IF(nvar .GT. 2)THEN
          CALL ifft( ifft3_in, ifft3_out, nfft )
        END IF

! Store fourier coefficients discard zero-padding
        DO k = 1,ncell(3)
          mesh_fft(1,i,j,k) = ifft1_out(k)/REAL(nfft,MK)
          IF(nvar .GT. 1)THEN
            mesh_fft(2,i,j,k) = ifft2_out(k)/REAL(nfft,MK)
          END IF
          IF(nvar .GT. 2)THEN
            mesh_fft(3,i,j,k) = ifft3_out(k)/REAL(nfft,MK)
          END IF
        END DO

      END DO !i
    END DO !j

  END IF

!---------------------------------------------------------------------------------!
! De-allocate fft arrays
!---------------------------------------------------------------------------------!
  DEALLOCATE( fft1_in, fft1_out, ifft1_in, ifft1_out )
  IF(nvar .GT. 1)THEN
    DEALLOCATE( fft2_in, fft2_out, ifft2_in, ifft2_out )
  END IF
  IF(nvar .GT. 2)THEN
    DEALLOCATE( fft3_in, fft3_out, ifft3_in, ifft3_out )
  END IF

!---------------------------------------------------------------------------------!
! IFFT velocity field in the first direction
!---------------------------------------------------------------------------------!
  ifft_seq(1) = patch%fft_seq(2)
  ifft_seq(2) = patch%fft_seq(1)
  ifft_seq(3) = 0
  CALL pmlib_ifft( topo_all,topo_third,mesh_fft,topo_all%cuboid,mesh, &
                 & ifft_seq,ierr,fft_shift = .FALSE.)

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh_fft) ) THEN
    DEALLOCATE(mesh_fft)
  END IF

!---------------------------------------------------------------------------------!
! Deallocate auxilirary topology
!---------------------------------------------------------------------------------!
  DEALLOCATE(topo_third)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_regularise
