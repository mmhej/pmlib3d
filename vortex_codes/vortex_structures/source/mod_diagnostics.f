!---------------------------------------------------------------------------------!
! mod_diagnostics.f
!---------------------------------------------------------------------------------!
! This routine initialises the vorticity field
!---------------------------------------------------------------------------------!
MODULE mod_diagnostics

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

  INTERFACE diagnostics_energy
    MODULE PROCEDURE diagnostics_energy_part
    MODULE PROCEDURE diagnostics_energy_mesh
  END INTERFACE

CONTAINS

#include "diagnostics/diagnostics_energy_part.f"
#include "diagnostics/diagnostics_energy_mesh.f"

#include "diagnostics/diagnostics_ring.f"

#include "diagnostics/diagnostics_bump.f"

#include "diagnostics/diagnostics_collision.f"

#include "diagnostics/diagnostics_helix.f"

#include "diagnostics/diagnostics_knot.f"

#include "diagnostics/diagnostics_sphere.f"

#include "diagnostics/diagnostics_spectral.f"

END MODULE mod_diagnostics

