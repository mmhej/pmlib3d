!---------------------------------------------------------------------------------!
! pmlib_poisson.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_poisson

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "poisson/pmlib_poisson_setup.f"

#ifdef __mkl
#include "poisson/pmlib_poisson_solve_mkl.f"
#elif __fftw
#include "poisson/pmlib_poisson_solve_fftw.f"
#else
#include "poisson/pmlib_poisson_solve.f"
#endif


END MODULE pmlib_mod_poisson
