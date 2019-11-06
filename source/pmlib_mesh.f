!---------------------------------------------------------------------------------!
! pmlib_mesh.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_mesh

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Define field class
!---------------------------------------------------------------------------------!
  TYPE class_mesh
! velocity (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: vel => NULL()
! vorticity (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: vort => NULL()
! vector potential (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: vec_pot => NULL()
! material derivative of vorticity (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: dvort => NULL()
! poisson kernel fft
    COMPLEX(MKC),DIMENSION(:,:,:,:),POINTER :: pkernel_fft  => NULL()
! patch mask
    LOGICAL,DIMENSION(:,:,:),POINTER        :: patch_mask => NULL()

! Penalisation variables
! penalisation velocity (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER      :: vel_pen => NULL()
! characteristic function (1,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER      :: solid => NULL()
! velocity of solid (dim,x,y,z)
    REAL(MK),DIMENSION(:,:,:,:),POINTER      :: solid_vel => NULL()
  END TYPE class_mesh

!---------------------------------------------------------------------------------!
! Module variables
!---------------------------------------------------------------------------------!
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE      :: pmlib_mesh_mean => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE      :: pmlib_mesh_parent => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE      :: pmlib_mesh_aux => NULL()

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "mesh/pmlib_mesh_allocate.f"
#include "mesh/pmlib_mesh_map.f"
#include "mesh/pmlib_mesh_map_ghost.f"

END MODULE pmlib_mod_mesh
