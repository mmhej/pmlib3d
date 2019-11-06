!---------------------------------------------------------------------------------!
! pmlib_interpolation.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_interpolation

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

#include "interpolation/pmlib_interp_particle_mesh.f"
#include "interpolation/pmlib_interp_mesh_particle.f"
#include "interpolation/pmlib_interp_mesh_mesh.f"

END MODULE pmlib_mod_interpolation
