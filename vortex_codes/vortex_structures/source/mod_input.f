!---------------------------------------------------------------------------------!
! mod_input.f
!---------------------------------------------------------------------------------!
! This module contains the input routines
!---------------------------------------------------------------------------------!
MODULE mod_input

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

CONTAINS

#include "input/input_setup.f"
#include "input/input_particles.f"

END MODULE mod_input

