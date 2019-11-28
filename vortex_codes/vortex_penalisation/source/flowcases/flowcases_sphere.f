!---------------------------------------------------------------------------------!
! flowcases_sphere.f
!---------------------------------------------------------------------------------!
! This routine sets up the penalisation routine
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_sphere(patch,ierr)

USE mod_penalisation

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_poisson

USE poisson_solver_module

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                            :: patch
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=16)                     :: caller = 'penal_setup'
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  INTEGER                                      :: nbuff, nvar

  REAL(MK),DIMENSION(ndim)                     :: xmin, xmax, dx
  INTEGER,DIMENSION(ndim)                      :: ncell, offset
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  REAL(MK)                                     :: r,r0

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr  = 0
  nvar  = 3
  nbuff = 5
  r0    = 0.5_MK

!---------------------------------------------------------------------------------!
! Create penalisation patch
!---------------------------------------------------------------------------------!
  patch_penal%dx         = patch%dx
  patch_penal%ptype      = 1
  patch_penal%level      = 1
  patch_penal%parent     = (/ 0, 0 /)

  patch_penal%nghost     = 4
  patch_penal%bound_cond = (/ 0, 0, 0 /)

!---------------------------------------------------------------------------------!
! Find bounding box of the solid budy
!---------------------------------------------------------------------------------!
  patch_penal%xmin(1) =   0.0_MK - REAL(nbuff,MK) * patch%dx(1)
  patch_penal%xmin(2) = - 0.5_MK - REAL(nbuff,MK) * patch%dx(2)
  patch_penal%xmin(3) = - 0.5_MK - REAL(nbuff,MK) * patch%dx(3)
  patch_penal%xmax(1) =   1.0_MK + REAL(nbuff,MK) * patch%dx(1)
  patch_penal%xmax(2) =   0.5_MK + REAL(nbuff,MK) * patch%dx(2)
  patch_penal%xmax(3) =   0.5_MK + REAL(nbuff,MK) * patch%dx(3)

!---------------------------------------------------------------------------------!
! Adjust penalisation patch
!---------------------------------------------------------------------------------!
  CALL pmlib_patch_adjust(patch_penal,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to adjust patches.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Correct penalisation patch to the main mesh
!---------------------------------------------------------------------------------!
  DO i = 1,ndim
    patch_penal%xmin(i) = REAL( NINT( ( patch_penal%xmin(i) - patch%xmin(i) ) &
                        & / patch%dx(i) ),MK) * patch%dx(i) + patch%xmin(i)
    patch_penal%xmax(i) = patch_penal%xmin(i) &
                      & + REAL(patch_penal%ncell(i),MK) * patch_penal%dx(i)
  END DO

!---------------------------------------------------------------------------------!
! Create penalisation topology
!---------------------------------------------------------------------------------!
  CALL pmlib_topology_cuboid(patch_penal,topo_penal%cuboid,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to create topology.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Allocate penalisation mesh
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh_penal%vort) ) DEALLOCATE( mesh_penal%vort)
  IF( ASSOCIATED(mesh_penal%vel) ) DEALLOCATE( mesh_penal%vel)
  IF( ASSOCIATED(mesh_penal%vel_pen) ) DEALLOCATE( mesh_penal%vel_pen)
  IF( ASSOCIATED(mesh_penal%solid)) DEALLOCATE( mesh_penal%solid)
  IF( ASSOCIATED(mesh_penal%solid_vel)) DEALLOCATE( mesh_penal%solid_vel)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to deallocate array.')
    GOTO 9999
  ENDIF

  dx     = topo_penal%cuboid(rank)%dx
  xmin   = topo_penal%cuboid(rank)%xmin
  ncell  = topo_penal%cuboid(rank)%ncell
  nghost = topo_penal%cuboid(rank)%nghost

  ALLOCATE( mesh_penal%vort(nvar, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)

  ALLOCATE( mesh_penal%vel(ndim, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)

  ALLOCATE( mesh_penal%vel_pen(ndim, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)

  ALLOCATE( mesh_penal%solid(1, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)

  ALLOCATE( mesh_penal%solid_vel(ndim, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)

  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to allocate array.')
    GOTO 9999
  ENDIF

  mesh_penal%vort      = 0.0_MK
  mesh_penal%vel       = 0.0_MK
  mesh_penal%vel_pen   = 0.0_MK
  mesh_penal%solid     = 0.0_MK
  mesh_penal%solid_vel = 0.0_MK

!---------------------------------------------------------------------------------!
! Create characteristic function
!---------------------------------------------------------------------------------!
  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)
        
        r = SQRT( (px - 0.5_MK)**2 + py**2 + pz**2)

        IF(r .LE. r0)THEN
          mesh_penal%solid(1,i,j,k) = 1.0_MK
        ELSE
          mesh_penal%solid(1,i,j,k) = 0.0_MK
        END IF

        mesh_penal%solid_vel(1,i,j,k) = 0.0_MK
        mesh_penal%solid_vel(2,i,j,k) = 0.0_MK
        mesh_penal%solid_vel(3,i,j,k) = 0.0_MK
      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! Setup Poisson solver on the penalisation mesh
!---------------------------------------------------------------------------------!
! Store cuboid partition to solver 2
    ALLOCATE( poisson_solver(2)%partition( 0:nproc-1 ) ) 
    offset = (/ 1, 1, 1 /)
    DO i = 0,nproc-1
      poisson_solver(2)%partition(i)%ncell = topo_penal%cuboid(i)%ncell
      poisson_solver(2)%partition(i)%icell = topo_penal%cuboid(i)%icell-offset
      poisson_solver(2)%partition(i)%dx    = topo_penal%cuboid(i)%dx
    END DO

! Setup solver 2
    CALL poisson_solver_setup3d(2,patch%ncell,patch%bound_cond,patch%dx)
    CALL poisson_solver_set_return_curl(2,.TRUE.) ! specify lhs operator

!  CALL pmlib_poisson_setup(patch_penal,topo_penal,mesh_penal,ierr)
!  IF (ierr .NE. 0) THEN
!    CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
!    GOTO 9999
!  ENDIF

!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .TRUE.

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE flowcases_sphere
