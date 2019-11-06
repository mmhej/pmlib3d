!---------------------------------------------------------------------------------!
! pmlib_repatch.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_repatch


USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!
  INTERFACE pmlib_repatch
    MODULE PROCEDURE pmlib_repatch
  END INTERFACE

CONTAINS

#include "repatch/pmlib_repatch.f"

END MODULE pmlib_mod_repatch
