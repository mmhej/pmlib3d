!---------------------------------------------------------------------------------!
! pmlib_patch.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_patch

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Define patch class
!---------------------------------------------------------------------------------!
  TYPE class_patch
    INTEGER                          :: ptype      !patch type
    INTEGER                          :: level      !resolution level
    INTEGER                          :: nsubs      !total number of subdomains
    INTEGER, DIMENSION(2)            :: parent     !patch parent
    REAL(MK),DIMENSION(pmlib_ndim)   :: xmin       !minimum physical extent
    REAL(MK),DIMENSION(pmlib_ndim)   :: xmax       !maximum physical extent
    INTEGER, DIMENSION(pmlib_ndim)   :: ncell      !number of mesh cells
    REAL(MK),DIMENSION(pmlib_ndim)   :: dx         !resolution
    INTEGER, DIMENSION(pmlib_ndim)   :: bound_cond !boundary conditions
    INTEGER                          :: nghost     !number of ghost points
    INTEGER, DIMENSION(pmlib_ndim)   :: fft_seq    !sequence in which to do ffts
  END TYPE class_patch

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "patch/pmlib_patch_adjust.f"

END MODULE pmlib_mod_patch
