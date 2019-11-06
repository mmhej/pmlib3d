!---------------------------------------------------------------------------------!
! pmlib_poisson.f
!   version: 2D single-core
!---------------------------------------------------------------------------------!
! <<HEADER>>
!---------------------------------------------------------------------------------!
MODULE mod_poisson

USE pmlib_mod_poisson

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!
!  INTERFACE pmlib_poisson_setup
!    MODULE PROCEDURE pmlib_poisson_setup_sp
!    MODULE PROCEDURE pmlib_poisson_setup_mp
!  END INTERFACE pmlib_poisson_setup

!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "poisson/poisson_setup_ppf.f"

END MODULE mod_poisson
