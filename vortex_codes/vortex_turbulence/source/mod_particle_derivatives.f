!---------------------------------------------------------------------------------!
! mod_particle_derivatives.f
!---------------------------------------------------------------------------------!
! This module contains the particle derivaties of the particles.
!---------------------------------------------------------------------------------!
MODULE mod_particle_derivatives

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

CONTAINS

#include "particle_derivatives/part_deriv_vorticity.f"

END MODULE mod_particle_derivatives
