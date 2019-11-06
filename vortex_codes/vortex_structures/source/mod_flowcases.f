!---------------------------------------------------------------------------------!
! mod_flowcases.f
!---------------------------------------------------------------------------------!
! This routine initialises the vorticity field
!---------------------------------------------------------------------------------!
MODULE mod_flowcases

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

CONTAINS

#include "flowcases/flowcases_ring.f"

#include "flowcases/flowcases_gauss.f"

#include "flowcases/flowcases_ringfit.f"

#include "flowcases/flowcases_bump.f"

#include "flowcases/flowcases_collision.f"

#include "flowcases/flowcases_helix.f"

#include "flowcases/flowcases_knot.f"

END MODULE mod_flowcases

