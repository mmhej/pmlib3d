!---------------------------------------------------------------------------------!
! pmlib_fourier.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_fourier

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_communication

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Module variables
!---------------------------------------------------------------------------------!
  COMPLEX(MKC),DIMENSION(:),POINTER      :: fft_tmp

CONTAINS

#ifdef __mkl

#include "fourier/pmlib_fourier_fft_mkl.f"
#include "fourier/pmlib_fourier_ifft_mkl.f" 

#elif __fftw

#include "fourier/pmlib_fourier_fft_fftw.f"
#include "fourier/pmlib_fourier_ifft_fftw.f"

#else

#include "fourier/pmlib_fourier_fft.f"
#include "fourier/pmlib_fourier_ifft.f" 

#endif

END MODULE pmlib_mod_fourier
