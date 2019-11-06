!---------------------------------------------------------------------------------!
! pmlib_particles.f
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_particles

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Define particle class
!---------------------------------------------------------------------------------!
  TYPE class_particles
! position (dim,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: pos => NULL()
! velocity (dim,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: vel => NULL()
! vorticity (ncom,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: vort => NULL()
! material derivative of vorticity (ncom,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: dvort => NULL()

! substep velocity (dim,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: vel_rk => NULL()
! substep material derivative of vorticity (ncom,ipart)
    REAL(MK),DIMENSION(:,:),POINTER :: dvort_rk => NULL()
  END TYPE class_particles

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "particles/pmlib_particles_map.f"
#include "particles/pmlib_particles_advance.f"

END MODULE pmlib_mod_particles
