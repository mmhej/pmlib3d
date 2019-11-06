!---------------------------------------------------------------------------------!
! pmlib_topology_parent.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_topology_parent(topo,topo_parent,ierr)

USE pmlib_mod_patch

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN)     :: topo
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(OUT)    :: topo_parent
  INTEGER, INTENT(OUT)                                      :: ierr

!---------------------------------------------------------------------------------!
!  Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=21)                   :: caller = 'pmlib_topology_parent'
  INTEGER                             :: i, iproc
  INTEGER                             :: nghost_ext, nghost_int
  REAL(MK)                            :: reldx

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Allocate 
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(topo_parent) ) THEN
    DEALLOCATE(topo_parent,stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to deallocate topology.')
      GOTO 9999
    ENDIF
  END IF
  ALLOCATE( topo_parent(0:mpi_nproc-1) )

!---------------------------------------------------------------------------------!
! Pass patch related data to topology
!---------------------------------------------------------------------------------!
  DO iproc = 0,mpi_nproc-1
    topo_parent(iproc)%ptype      = topo(iproc)%ptype
    topo_parent(iproc)%level      = topo(iproc)%level
    topo_parent(iproc)%parent     = topo(iproc)%parent

    IF( topo(iproc)%parent(1) .EQ. 0 )THEN
      ierr = - 1
      CALL pmlib_write(mpi_rank,caller,'Cannot create parent topology on level 0' )
      GO TO 9999
    END IF

    topo_parent(iproc)%xmin       = topo(iproc)%xmin
    topo_parent(iproc)%xmax       = topo(iproc)%xmax

    reldx = 2.0_MK**(topo(iproc)%level - topo(iproc)%parent(1))

    DO i = 1,pmlib_ndim
      IF( MOD(topo(iproc)%ncell(i), NINT(reldx) ) .EQ. 0 )THEN
        topo_parent(iproc)%ncell(i) = NINT( REAL(topo(iproc)%ncell(i), MK)/reldx )
      ELSE
        ierr = - 1
        CALL pmlib_write(mpi_rank,caller,'Number mesh cells does not match' )
        GO TO 9999
      END IF
    END DO
    topo_parent(iproc)%dx         = reldx * topo(iproc)%dx

!---------------------------------------------------------------------------------!
! Number of ghost cells
!---------------------------------------------------------------------------------!
    IF(pmlib_interpolation_order .EQ. 3)THEN
      topo_parent(iproc)%nghost = 2
    ELSEIF(pmlib_interpolation_order .EQ. 4)THEN
      topo_parent(iproc)%nghost = 3
    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Can not define number of ghost cells.')
      GOTO 9999    
    END IF

!---------------------------------------------------------------------------------!
! Handle boundaries
!---------------------------------------------------------------------------------!
    topo_parent(iproc)%bound_cond = topo(iproc)%bound_cond
    topo_parent(iproc)%ineigh     = topo(iproc)%ineigh

  END DO ! iproc

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_topology_parent
