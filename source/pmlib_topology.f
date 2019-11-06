!---------------------------------------------------------------------------------!
! pmlib_topology.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_topology

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Define topology class
!---------------------------------------------------------------------------------!
  TYPE class_topology
    INTEGER                          :: ptype      !patch type
    INTEGER                          :: level      !resolution level
    INTEGER, DIMENSION(2)            :: parent     !patch parent

    REAL(MK),DIMENSION(pmlib_ndim)   :: xmin     !minimum physical extent
    REAL(MK),DIMENSION(pmlib_ndim)   :: xmax     !maximum physical extent
    INTEGER, DIMENSION(pmlib_ndim)   :: ncell    !number of mesh cells
    INTEGER, DIMENSION(pmlib_ndim)   :: icell    !global index of first mesh cell
    REAL(MK),DIMENSION(pmlib_ndim)   :: dx       !resolution
    INTEGER, DIMENSION(2*pmlib_ndim) :: bound_cond !boundary conditions
    INTEGER, DIMENSION(2*pmlib_ndim) :: nghost   !number of ghost points
    INTEGER                          :: npart    !number of particles (local?)

    INTEGER, DIMENSION(2*pmlib_ndim) :: ineigh ! proc(west,east,south,...)
  END TYPE class_topology

  TYPE class_topology_all
    TYPE(class_topology),DIMENSION(:),POINTER     :: cuboid  => NULL()
    TYPE(class_topology),DIMENSION(:),POINTER     :: xpencil => NULL()
    TYPE(class_topology),DIMENSION(:),POINTER     :: ypencil => NULL()
    TYPE(class_topology),DIMENSION(:),POINTER     :: zpencil => NULL()
    TYPE(class_topology),DIMENSION(:),POINTER     :: parent  => NULL()
  END TYPE class_topology_all

!---------------------------------------------------------------------------------!
! Interfaces
!---------------------------------------------------------------------------------!
  INTERFACE pmlib_topology_cuboid
    MODULE PROCEDURE pmlib_topology_cuboid
  END INTERFACE

  INTERFACE pmlib_topology_pencil
    MODULE PROCEDURE pmlib_topology_pencil
  END INTERFACE

  INTERFACE pmlib_topology_kernel
    MODULE PROCEDURE pmlib_topology_kernel
  END INTERFACE

  INTERFACE pmlib_topology_parent
    MODULE PROCEDURE pmlib_topology_parent
  END INTERFACE

!---------------------------------------------------------------------------------!
! Include subroutines
!---------------------------------------------------------------------------------!
CONTAINS

#include "topology/pmlib_topology_cuboid.f"

#include "topology/pmlib_topology_pencil.f"
#include "topology/pmlib_topology_kernel.f"
#include "topology/pmlib_topology_parent.f"


END MODULE pmlib_mod_topology
