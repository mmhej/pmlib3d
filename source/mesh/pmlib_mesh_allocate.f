!---------------------------------------------------------------------------------!
! pmlib_mesh_allocate.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_mesh_allocate(topo,mesh,ierr,vort,dvort,vel,mask)

USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER :: topo
  TYPE(class_mesh)                          :: mesh
  INTEGER, INTENT(OUT)                      :: ierr
  LOGICAL, INTENT(IN), OPTIONAL             :: vort
  LOGICAL, INTENT(IN), OPTIONAL             :: dvort
  LOGICAL, INTENT(IN), OPTIONAL             :: vel
  LOGICAL, INTENT(IN), OPTIONAL             :: mask

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=19)                         :: caller = 'pmlib_mesh_allocate'
  INTEGER,DIMENSION(pmlib_ndim)             :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)           :: nghost

!---------------------------------------------------------------------------------!
! Allocate the multi level, multi patch type (not the actual mesh)
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh%vort) ) DEALLOCATE( mesh%vort, stat=ierr)
  IF( ASSOCIATED(mesh%vel) )  DEALLOCATE( mesh%vel, stat=ierr)
  IF( ASSOCIATED(mesh%dvort)) DEALLOCATE( mesh%dvort, stat=ierr)

  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Allocate the mesh
!---------------------------------------------------------------------------------!
  ncell  = topo(mpi_rank)%ncell
  nghost = topo(mpi_rank)%nghost

  IF( PRESENT(vort) .AND. vort )THEN
    ALLOCATE( mesh%vort(pmlib_nvort, &
            & 1-nghost(1):ncell(1)+nghost(2),  &
            & 1-nghost(3):ncell(2)+nghost(4),  &
            & 1-nghost(5):ncell(3)+nghost(6)), &
            & stat=ierr)
  END IF

  IF( PRESENT(dvort) .AND. dvort )THEN
    ALLOCATE( mesh%dvort(pmlib_nvort, &
            & 1-nghost(1):ncell(1)+nghost(2),  &
            & 1-nghost(3):ncell(2)+nghost(4),  &
            & 1-nghost(5):ncell(3)+nghost(6)), &
            & stat=ierr)
  END IF

  IF( PRESENT(vel) .AND. vel )THEN
    ALLOCATE( mesh%vel(pmlib_ndim, &
            & 1-nghost(1):ncell(1)+nghost(2),  &
            & 1-nghost(3):ncell(2)+nghost(4),  &
            & 1-nghost(5):ncell(3)+nghost(6)), &
            & stat=ierr)    
  END IF

  IF( PRESENT(mask) .AND. mask )THEN
    ALLOCATE( mesh%patch_mask( &
            & 1-nghost(1):ncell(1)+nghost(2),  &
            & 1-nghost(3):ncell(2)+nghost(4),  &
            & 1-nghost(5):ncell(3)+nghost(6)), &
            & stat=ierr)    
  END IF

  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
    GOTO 9999
  ENDIF


!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_mesh_allocate
