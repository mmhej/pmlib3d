!---------------------------------------------------------------------------------!
! pmlib_visualise.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_visualise

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE


CONTAINS

#include "visualise/pmlib_visualise_particles.f"
#include "visualise/pmlib_visualise_sort_particles.f"
#include "visualise/pmlib_visualise_postscript.f"

END MODULE pmlib_mod_visualise
