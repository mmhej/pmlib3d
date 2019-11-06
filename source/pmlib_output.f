!---------------------------------------------------------------------------------!
! pmlib_output.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_output

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!
  INTERFACE pmlib_output_topology
    MODULE PROCEDURE pmlib_output_topology
  END INTERFACE

CONTAINS

#include "output/pmlib_output_topology.f"
#include "output/pmlib_output_mesh.f"

END MODULE pmlib_mod_output
