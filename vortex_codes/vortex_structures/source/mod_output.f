!---------------------------------------------------------------------------------!
! mod_output.f
!---------------------------------------------------------------------------------!
! This module contains the output routines
!---------------------------------------------------------------------------------!
MODULE mod_output

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

CONTAINS

#include "output/output_particles.f"
#include "output/output_mesh.f"

END MODULE mod_output

