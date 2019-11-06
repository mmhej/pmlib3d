!---------------------------------------------------------------------------------!
! mod_penalisation.f
!---------------------------------------------------------------------------------!
! This routine enforces a solid body by iterative Brinkman penalisation
!---------------------------------------------------------------------------------!
MODULE mod_penalisation

USE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Set penalisation classes
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                         :: patch_penal
  TYPE(class_topology_all)                  :: topo_penal
  TYPE(class_mesh),SAVE                     :: mesh_penal

CONTAINS

#include "penalisation/penalisation_solve.f"

END MODULE mod_penalisation
