!---------------------------------------------------------------------------------!
! pmlib_regularise.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_regularise

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!
  INTERFACE pmlib_regularise
    MODULE PROCEDURE pmlib_regularise
  END INTERFACE

CONTAINS

#ifdef __mkl

#include "regularise/pmlib_regularise_mkl.f"

#elif __fftw

#include "regularise/pmlib_regularise_fftw.f"

#else

#include "regularise/pmlib_regularise.f"

#endif

END MODULE pmlib_mod_regularise
