!---------------------------------------------------------------------------------!
! pmlib_parameters.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_parameters

USE pmlib_mod_write

IMPLICIT NONE

include 'mpif.h'
!---------------------------------------------------------------------------------!
! Arithmetic precision
!---------------------------------------------------------------------------------!
#ifdef __single
! Single precision
  INTEGER, PARAMETER           :: MK               = KIND(1E0)
  INTEGER, PARAMETER           :: MKC              = KIND(CMPLX(1E0,1E0,MK))
  INTEGER, PARAMETER           :: mpi_prec_real    = MPI_REAL
  INTEGER, PARAMETER           :: mpi_prec_complex = MPI_COMPLEX
#elif __double
! Double precision
  INTEGER, PARAMETER           :: MK               = KIND(1D0)
  INTEGER, PARAMETER           :: MKC              = KIND(CMPLX(1D0,1D0,MK))
  INTEGER, PARAMETER           :: mpi_prec_real    = MPI_DOUBLE
  INTEGER, PARAMETER           :: mpi_prec_complex = MPI_DOUBLE_COMPLEX
#else

#error "Arithmetic presicion is not defined"

#endif

  INTEGER, PARAMETER               :: mpi_prec_int     = MPI_INT
  INTEGER, PARAMETER               :: mpi_prec_log     = MPI_LOGICAL
  INTEGER, PARAMETER               :: MAXCHAR          = 256

!---------------------------------------------------------------------------------!
! Number of dimensions
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER               :: pmlib_ndim  = 3
  INTEGER, PARAMETER               :: pmlib_nvort = 3

!---------------------------------------------------------------------------------!
! Patch structure
!---------------------------------------------------------------------------------!
  INTEGER                          :: pmlib_nlevel = 1
  INTEGER                          :: pmlib_max_npatch = 1
  INTEGER,DIMENSION(:),ALLOCATABLE :: pmlib_npatch
  INTEGER                          :: pmlib_patch_ext = 32

!---------------------------------------------------------------------------------!
! Numerical parameters
!---------------------------------------------------------------------------------!
  REAL(MK)                         :: pmlib_regularisation_radius = 0.0_MK
  INTEGER                          :: pmlib_poisson_kernel = -99
  INTEGER                          :: pmlib_poisson_order = -99
  INTEGER                          :: pmlib_fd_order = -99
  INTEGER                          :: pmlib_interpolation_order = -99
  INTEGER                          :: pmlib_time_integration_order = -99

!---------------------------------------------------------------------------------!
! MPI parameters
!---------------------------------------------------------------------------------!
  INTEGER                          :: mpi_comm  = -99
  INTEGER                          :: mpi_rank  = -99
  INTEGER                          :: mpi_nproc = -99

CONTAINS

#include "parameters/pmlib_parameters_check.f"

END MODULE pmlib_mod_parameters
