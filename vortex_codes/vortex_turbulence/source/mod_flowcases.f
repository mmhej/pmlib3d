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

#include "flowcases/flowcases_turbulence.f"

#include "flowcases/flowcases_taylor_green.f"

END MODULE mod_flowcases
