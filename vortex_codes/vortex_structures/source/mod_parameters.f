!---------------------------------------------------------------------------------!
! mod_parameters.f
!---------------------------------------------------------------------------------!
! This routines declares all data that are global to the client.
!---------------------------------------------------------------------------------!
MODULE mod_parameters

USE pmlib_mod_parameters
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_particles
USE pmlib_mod_mesh


!---------------------------------------------------------------------------------!
! Flow case
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=256)               :: flowcase
  LOGICAL                          :: restart
  CHARACTER(LEN=256)               :: restart_file
  LOGICAL                          :: abort

!---------------------------------------------------------------------------------!
! Dimensions
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER               :: ndim  = pmlib_ndim
  INTEGER, PARAMETER               :: nvort = pmlib_nvort

!---------------------------------------------------------------------------------!
! MPI parameters
!---------------------------------------------------------------------------------! 
  INTEGER                          :: comm
  INTEGER                          :: rank
  INTEGER                          :: nproc

!---------------------------------------------------------------------------------!
! Physical parameters
!---------------------------------------------------------------------------------!
  REAL(MK)                         :: length
  REAL(MK)                         :: velocity
  REAL(MK)                         :: visc

  REAL(MK), DIMENSION(ndim)        :: vel_inf

  REAL(MK)                         :: vort_max
  REAL(MK),DIMENSION(ndim)         :: vort_centroid

!---------------------------------------------------------------------------------!
! Time parameters
!---------------------------------------------------------------------------------!
  INTEGER                          :: ntime, itime
  REAL(MK)                         :: dtime, time
  REAL(MK)                         :: ctime

!---------------------------------------------------------------------------------!
! Simulation parameters
!---------------------------------------------------------------------------------!
  INTEGER                          :: resolution
  LOGICAL                          :: reproject
  LOGICAL                          :: regularise_vort
  LOGICAL                          :: penalisation

  INTEGER                          :: iplot_mesh, iplot_part
  INTEGER                          :: ioutput_part, ioutput_diag
  INTEGER                          :: iremesh, irepatch
  LOGICAL                          :: remesh, repatch
  REAL(MK)                         :: repatch_trunc, remesh_trunc

!---------------------------------------------------------------------------------!
! Penalisation parameters
!---------------------------------------------------------------------------------!
  REAL(MK)                         :: penal_tolerance
  REAL(MK)                         :: penal_relaxation
  INTEGER                          :: penal_max

!---------------------------------------------------------------------------------!
! Diagnostics files
!---------------------------------------------------------------------------------!
  INTEGER                               :: idiag, ndiag
  REAL(MK), DIMENSION(:,:), ALLOCATABLE :: diag


END MODULE mod_parameters
