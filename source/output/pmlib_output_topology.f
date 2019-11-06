!------------------------------------------------------------------------------!
! pmlib_output_topology.f
!   version: 3D multi-core
!------------------------------------------------------------------------------!
SUBROUTINE pmlib_output_topology(patch,topo,tag,ierr)

USE pmlib_mod_topology
USE pmlib_mod_patch

IMPLICIT NONE
!-------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------
  TYPE(class_patch)                         :: patch
  TYPE(class_topology),DIMENSION(:),POINTER :: topo
  CHARACTER(LEN=*), INTENT(IN)              :: tag
  INTEGER, INTENT(OUT)                      :: ierr

!-------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------
  INTEGER                            :: i
  CHARACTER(LEN=64)                  :: topofile

!----------------------------------------------------------------------------
! Initialise routine
!----------------------------------------------------------------------------
  ierr = 0

  WRITE(topofile,'(A,A)') TRIM(tag),'.dat'
  OPEN(17, FILE=TRIM(topofile))
!----------------------------------------------------------------------------
! Write file
!----------------------------------------------------------------------------
  WRITE(17,'(A)')'#==========================================================================='
  WRITE(17,'(A)')'## Topology'
  WRITE(17,'(A)')'#==========================================================================='
  WRITE(17,'(A)')''
  WRITE(17,'(A,I2)')   '# Patch type:             ', patch%ptype
  WRITE(17,'(A,3F8.2)')'# Minimum extend:         ', patch%xmin
  WRITE(17,'(A,3F8.2)')'# Maximum extend:         ', patch%xmax
  WRITE(17,'(A,I4)')   '# Resolution:                1/', NINT(1.0_MK/patch%dx(1))
  WRITE(17,'(A,3I8)')  '# Number of mesh cells:   ', patch%ncell
  WRITE(17,'(A,3I2)')  '# FFT-seq:                ', patch%fft_seq

  WRITE(17,'(A)')''
  DO i = 0,mpi_nproc-1

    WRITE(17,'(A)')   '#---------------------------------------------------------------------------'
    WRITE(17,'(A,I4)')'# Topology for processor: ', i
    WRITE(17,'(A)')   '#---------------------------------------------------------------------------'
    WRITE(17,'(A)')''
    WRITE(17,'(A,I8)')    '#  Number of particles:    ', topo(i)%npart
    WRITE(17,'(A,3F8.2)') '#  Minimum extend:         ', topo(i)%xmin
    WRITE(17,'(A,3F8.2)') '#  Maximum extend:         ', topo(i)%xmax
    WRITE(17,'(A,3I8)')   '#  Number of mesh cells:   ', topo(i)%ncell
    WRITE(17,'(A,6I4)')   '#  Number of ghost cells:  ', topo(i)%nghost
    WRITE(17,'(A,6I4)')   '#  Boundary conditions:    ', topo(i)%bound_cond
    WRITE(17,'(A,6I4)')   '#  Neighbours:             ', topo(i)%ineigh
    WRITE(17,'(A)')''

  END DO

  CLOSE(17)

!----------------------------------------------------------------------------
! Return
!----------------------------------------------------------------------------
! 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_output_topology
