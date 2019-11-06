!---------------------------------------------------------------------------------!
! pmlib_finalise.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_finalise

USE pmlib_mod_parameters
USE pmlib_mod_write

IMPLICIT NONE

  INTERFACE pmlib_finalise
    MODULE PROCEDURE pmlib_finalise
  END INTERFACE

CONTAINS


SUBROUTINE pmlib_finalise(part,mesh,ierr)

USE pmlib_mod_mesh
USE pmlib_mod_particles
USE pmlib_mod_topology
USE pmlib_mod_poisson

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_particles)                        :: part
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! De-allocate particles
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(part%pos) ) DEALLOCATE( part%pos )
  IF( ASSOCIATED(part%vel) ) DEALLOCATE( part%vel )
  IF( ASSOCIATED(part%vort) ) DEALLOCATE( part%vort )
  IF( ASSOCIATED(part%dvort)) DEALLOCATE( part%dvort )
  IF( ASSOCIATED(part%vel_rk) ) DEALLOCATE( part%vel_rk )
  IF( ASSOCIATED(part%dvort_rk) ) DEALLOCATE( part%dvort_rk )
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,'pmlib_allocate_mesh','Failed to deallocate array.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! De-allocate mesh
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh%vel) ) DEALLOCATE( mesh%vel )
  IF( ASSOCIATED(mesh%vort) ) DEALLOCATE( mesh%vort )
  IF( ASSOCIATED(mesh%dvort) ) DEALLOCATE( mesh%dvort )
  IF( ASSOCIATED(mesh%pkernel_fft) ) DEALLOCATE(mesh%pkernel_fft)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,'pmlib_allocate_mesh','Failed to deallocate array.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_finalise

END MODULE pmlib_mod_finalise
