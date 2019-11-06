!---------------------------------------------------------------------------------!
! diagnostics_spectral.f
!---------------------------------------------------------------------------------!
! Calculates and outputs the diagnostics of the flowcase
!---------------------------------------------------------------------------------!
SUBROUTINE diagnostics_spectral(patch,topo_all,mesh,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_fourier

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch),INTENT(IN)                 :: patch
  TYPE(class_topology_all),INTENT(IN)          :: topo_all
  TYPE(class_mesh),INTENT(IN)                  :: mesh
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                     :: PI
  INTEGER                                      :: i,j,k, ii
  CHARACTER(LEN=256)                           :: outfile
  INTEGER                                      :: nf
  REAL(MK)                                     :: fmax, df

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  INTEGER                                      :: i_ptc, j_ptc, k_ptc
  REAL(MK), DIMENSION(pmlib_ndim)              :: l_ptc, ks, dk, kmin, ksmp, ki
  REAL(MK), DIMENSION(pmlib_ndim)              :: xmin_ptc, xmax_ptc
  INTEGER, DIMENSION(pmlib_ndim)               :: bc_ptc, ncell_ptc
  REAL(MK)                                     :: kmag, vol
  REAL(MK)                                     :: energy, enstrophy

  REAL(MK)                                     :: S11, S22, S33
  REAL(MK)                                     :: O11, O22, O33

  COMPLEX(MKC)                                 :: vel_fft1,vel_fft2,vel_fft3
  COMPLEX(MKC)                                 :: vort_fft1,vort_fft2,vort_fft3

  COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER      :: vort_fft, vel_fft

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

  IF (ioutput_diag .NE. 0)THEN

    IF( MOD(itime,ioutput_diag) .EQ. 0 &
      & .OR. abort .OR. itime .EQ. ntime )THEN


!---------------------------------------------------------------------------------!
! Fourier transform vorticity and velocity fields
!---------------------------------------------------------------------------------!
  CALL pmlib_fft( topo_all,topo_all%cuboid,mesh%vort,topo_all%cuboid,vort_fft, &
                & patch%fft_seq,ierr,fft_shift = .FALSE.)

  CALL pmlib_fft( topo_all,topo_all%cuboid,mesh%vel,topo_all%cuboid,vel_fft, &
                & patch%fft_seq,ierr,fft_shift = .FALSE.)

!---------------------------------------------------------------------------------!
! Get topology
!---------------------------------------------------------------------------------!
  ncell    = topo_all%cuboid(mpi_rank)%ncell
  dx       = topo_all%cuboid(mpi_rank)%dx
  xmin     = topo_all%cuboid(mpi_rank)%xmin

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

  ksmp(1) = 1.0_MK/( 2.0_MK * dx(1) )
  ksmp(2) = 1.0_MK/( 2.0_MK * dx(2) )
  ksmp(3) = 1.0_MK/( 2.0_MK * dx(3) )

  kmin(1) = - ksmp(1) + REAL(i_ptc,MK) * dk(1)
  kmin(2) = - ksmp(2) + REAL(j_ptc,MK) * dk(2)
  kmin(3) = - ksmp(3) + REAL(k_ptc,MK) * dk(3)

  vol = patch%ncell(1)*patch%ncell(2)*patch%ncell(3) * dx(1)*dx(2)*dx(3)


  WRITE(outfile,'(A,I5.5,A,I5.5,A)') 'spectrum_I',itime,'_P',rank,'.dat' 
  OPEN(20,FILE = outfile)

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

      DO i = 1,ncell(1)
! Calculate wave number in x-direction (FFT-shifted)
        IF( i + i_ptc .LE. ncell_ptc(1)/2 )THEN
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) + ksmp(1)
        ELSE
          ki(1) = kmin(1) + REAL(i - 1,MK)*dk(1) - ksmp(1)
        END IF

        kmag = SQRT( ki(1)**2 + ki(2)**2 + ki(3)**2 )

        vel_fft1 = vel_fft(1,i,j,k) * dx(1)*dx(2)*dx(3)
        vel_fft2 = vel_fft(2,i,j,k) * dx(1)*dx(2)*dx(3)
        vel_fft3 = vel_fft(3,i,j,k) * dx(1)*dx(2)*dx(3)   

        vort_fft1 = vort_fft(1,i,j,k) * dx(1)*dx(2)*dx(3)
        vort_fft2 = vort_fft(2,i,j,k) * dx(1)*dx(2)*dx(3)
        vort_fft3 = vort_fft(3,i,j,k) * dx(1)*dx(2)*dx(3)          

        S11 = REAL( vel_fft1 * vel_fft1 ,MK)/vol
        S22 = REAL( vel_fft2 * vel_fft2 ,MK)/vol
        S33 = REAL( vel_fft3 * vel_fft3 ,MK)/vol

        O11 = REAL( vort_fft1 * vort_fft1 ,MK)/vol
        O22 = REAL( vort_fft2 * vort_fft2 ,MK)/vol
        O33 = REAL( vort_fft3 * vort_fft3 ,MK)/vol

        energy    = 2.0_MK * PI * kmag**2 * (S11 + S22 + S33)
        enstrophy = 2.0_MK * PI * kmag**2 * (O11 + O22 + O33)

        WRITE(20,'(7E20.12)') kmag, energy, enstrophy

      END DO
    END DO
  END DO

  CLOSE(20)

!---------------------------------------------------------------------------------!
! Output gnuplot
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0)THEN

    WRITE(outfile,'(A)') 'spectrum.gnu' 
    OPEN(17,FILE = outfile)


    WRITE(17,'(A)') '# Spectrum PLOT'    

    WRITE(17,'(A)')'reset'
    WRITE(17,'(A)')'set terminal postscript eps enhanced'
    WRITE(17,'(A)')'set output "energy_spectrum.eps"'
    WRITE(17,'(A)')'set size 0.65,0.65'
    WRITE(17,'(A)')'set pointsize 1'
    WRITE(17,'(A)')'set log x'
    WRITE(17,'(A)')'set log y'
    WRITE(17,'(A)')'set format x "10^{%T}"'
    WRITE(17,'(A)')'set format y "10^{%T}"'
    WRITE(17,'(A)')'set yl "S"'
    WRITE(17,'(A)')'set xl "k"'
    WRITE(17,'(A)')'set xran [ 1e0 : 1e2 ]'
    WRITE(17,'(A)')'set yran [ 1e-10 : 1e-1 ]'
    WRITE(17,'(A)')''
    WRITE(17,'(A,E12.4)')'L = ', length
    WRITE(17,'(A,E12.4)')'std = ', velocity
    WRITE(17,'(A)')'scale = 3. * std**2 * L /0.32866'
    WRITE(17,'(A)')''
    WRITE(17,'(A)')'plot \'     !'(cancel compiler warning)
    WRITE(17,'(A)')'x**(-5/3) with line linetype 2 lc rgb "black" lw 1 notitle,\' !'
    DO i = 0,nproc-1
      WRITE(17,'(A,I5.5,A,I5.5,2A)') '"spectrum_I',itime,'_P',i,'.dat"', &
           & ' u 1:2 with points pointtype 1 lc rgb "blue" lw 1 notitle ,\'  !'
    END DO
    WRITE(17,'(2A)')'scale * L**4 * (2.*pi*x)**4 /(1. + L**2 * (2.*pi*x)**2)**(17./6.)', &
                  & ' with line linetype 1 lc rgb "red" lw 1 notitle'

    CLOSE(17)

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



    END IF
  END IF


!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE diagnostics_spectral
