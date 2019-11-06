!---------------------------------------------------------------------------------!
! pmlib_remesh.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_remesh( patch,topo,mesh,part,min_vort,min_dvort,ierr, &
                       & dvort,vel,mask,incl_ghost)

USE pmlib_mod_mesh
USE pmlib_mod_particles
USE pmlib_mod_patch
USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch), INTENT(IN)                              :: patch
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(INOUT)   :: topo
  TYPE(class_mesh), INTENT(IN)                               :: mesh
  TYPE(class_particles), INTENT(INOUT)                       :: part
  REAL(MK), INTENT(IN)                                       :: min_vort
  REAL(MK), INTENT(IN)                                       :: min_dvort
  INTEGER, INTENT(OUT)                                       :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                              :: dvort
  LOGICAL, INTENT(IN), OPTIONAL                              :: vel
  LOGICAL, INTENT(IN), OPTIONAL                              :: mask
  LOGICAL, INTENT(IN), OPTIONAL                              :: incl_ghost

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=22)                       :: caller = 'pmlib_remesh_particles'
  REAL(MK)                                :: px, py, pz
  REAL(MK)                                :: strength, dstrength
  INTEGER                                 :: i, j, k, l
  INTEGER                                 :: inp
  REAL(MK),DIMENSION(pmlib_ndim)          :: xmin, dx
  INTEGER,DIMENSION(pmlib_ndim)           :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)         :: nghost,bound_cond
  LOGICAL                                 :: place, ghost

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

  dx         = topo(mpi_rank)%dx
  ncell      = topo(mpi_rank)%ncell
  bound_cond = topo(mpi_rank)%bound_cond
  xmin       = topo(mpi_rank)%xmin
  nghost     = topo(mpi_rank)%nghost

  IF(PRESENT(incl_ghost) .AND. incl_ghost)THEN
    ghost = .TRUE.
  ELSE
    ghost = .FALSE.
  END IF

!---------------------------------------------------------------------------------!
! Deallocate the particles if they are already allocated
!---------------------------------------------------------------------------------!
  IF(ASSOCIATED( part%pos )) DEALLOCATE(part%pos,stat=ierr)

  IF(ASSOCIATED( part%vel )) DEALLOCATE(part%vel,stat=ierr)

  IF(ASSOCIATED( part%vort )) DEALLOCATE(part%vort,stat=ierr)

  IF(ASSOCIATED( part%dvort )) DEALLOCATE(part%dvort,stat=ierr)

  IF(ASSOCIATED( part%vel_rk )) DEALLOCATE(part%vel_rk,stat=ierr)

  IF(ASSOCIATED( part%dvort_rk )) DEALLOCATE(part%dvort_rk,stat=ierr)

  IF(ierr .NE. 0) THEN
    CALL pmlib_write( mpi_rank,'pmlib_remesh_particles', &
                    & 'Failed to deallocate existing particles.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Find the number of particles
!---------------------------------------------------------------------------------!
  inp = 0

  DO k = 1-nghost(5), ncell(3)+nghost(6)
    DO j = 1-nghost(3), ncell(2)+nghost(4)
      DO i = 1-nghost(1), ncell(1)+nghost(2)

        IF( ( ( i .GE. 1 .AND. i .LE. ncell(1) ) .OR. &
          &   ( ghost .AND. i .LT. 1 .AND. bound_cond(1) .EQ. 0 ) .OR. &
          &   ( ghost .AND. i .GT. ncell(1) .AND. bound_cond(2) .EQ. 0 ) ) .AND. &
          & ( ( j .GE. 1 .AND. j .LE. ncell(2) ) .OR. &
          &   ( ghost .AND. j .LT. 1 .AND. bound_cond(3) .EQ. 0 ) .OR. &
          &   ( ghost .AND. j .GT. ncell(2) .AND. bound_cond(4) .EQ. 0 ) ) .AND. &
          & ( ( k .GE. 1 .AND. k .LE. ncell(3) ) .OR. &
          &   ( ghost .AND. k .LT. 1 .AND. bound_cond(5) .EQ. 0 ) .OR. &
          &   ( ghost .AND. k .GT. ncell(3) .AND. bound_cond(6) .EQ. 0 ) ) )THEN
          place = .TRUE.
        ELSE
          place = .FALSE.
        END IF

        strength = SQRT( &
                 & mesh%vort(1,i,j,k)**2 + &
                 & mesh%vort(2,i,j,k)**2 + &
                 & mesh%vort(3,i,j,k)**2 &
                 & )

        dstrength = SQRT( &
                 & mesh%dvort(1,i,j,k)**2 + &
                 & mesh%dvort(2,i,j,k)**2 + &
                 & mesh%dvort(3,i,j,k)**2 )

        IF( place .AND. &
          & ( strength .GT. min_vort .OR. &
          &   dstrength .GT. min_dvort ) ) THEN
          IF(PRESENT(mask) .AND. mask)THEN
            place = mesh%patch_mask(i,j,k)
          ELSE
            place = .TRUE.
          END IF
        ELSE
          place = .FALSE.
        END IF

        IF(place) THEN
          inp = inp + 1
        END IF

      END DO !i
    END DO !j
  END DO !k

  topo(mpi_rank)%npart = inp

!---------------------------------------------------------------------------------!
! Allocate particles
!---------------------------------------------------------------------------------!
  ALLOCATE( part%pos(pmlib_ndim,inp), part%vort(pmlib_nvort,inp), &
          & part%vel(pmlib_ndim,inp), part%dvort(pmlib_nvort,inp) )

  ALLOCATE( part%vel_rk(pmlib_ndim,inp), &
          & part%dvort_rk(pmlib_nvort,inp) )

!---------------------------------------------------------------------------------!
! Create particles
!---------------------------------------------------------------------------------!
  inp = 0

  DO k = 1-nghost(5), ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)

    DO j = 1-nghost(3), ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)

      DO i = 1-nghost(1), ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        IF( ( ( i .GE. 1 .AND. i .LE. ncell(1) ) .OR. &
          &   ( ghost .AND. i .LT. 1 .AND. bound_cond(1) .EQ. 0 ) .OR. &
          &   ( ghost .AND. i .GT. ncell(1) .AND. bound_cond(2) .EQ. 0 ) ) .AND. &
          & ( ( j .GE. 1 .AND. j .LE. ncell(2) ) .OR. &
          &   ( ghost .AND. j .LT. 1 .AND. bound_cond(3) .EQ. 0 ) .OR. &
          &   ( ghost .AND. j .GT. ncell(2) .AND. bound_cond(4) .EQ. 0 ) ) .AND. &
          & ( ( k .GE. 1 .AND. k .LE. ncell(3) ) .OR. &
          &   ( ghost .AND. k .LT. 1 .AND. bound_cond(5) .EQ. 0 ) .OR. &
          &   ( ghost .AND. k .GT. ncell(3) .AND. bound_cond(6) .EQ. 0 ) ) )THEN
          place = .TRUE.
        ELSE
          place = .FALSE.
        END IF

        strength = SQRT( &
                 & mesh%vort(1,i,j,k)**2 + &
                 & mesh%vort(2,i,j,k)**2 + &
                 & mesh%vort(3,i,j,k)**2 &
                 & )

        dstrength = SQRT( &
                 & mesh%dvort(1,i,j,k)**2 + &
                 & mesh%dvort(2,i,j,k)**2 + &
                 & mesh%dvort(3,i,j,k)**2 )

        IF( place .AND. &
          & ( strength .GT. min_vort .OR. &
          &   dstrength .GT. min_dvort ) ) THEN
          IF(PRESENT(mask) .AND. mask)THEN
            place = mesh%patch_mask(i,j,k)
          ELSE
            place = .TRUE.
          END IF
        ELSE
          place = .FALSE.
        END IF

        IF(place) THEN

          inp = inp + 1

          part%pos(1,inp) = px
          part%pos(2,inp) = py
          part%pos(3,inp) = pz

          IF( PRESENT(vel) .AND. vel )THEN
            part%vel(1:pmlib_ndim,inp) = mesh%vel(1:pmlib_ndim,i,j,k)
          ELSE
            part%vel(1:pmlib_ndim,inp) = 0.0_MK
          END IF

          part%vel_rk(1:pmlib_ndim,inp) = 0.0_MK

          part%vort(1:pmlib_nvort,inp) = mesh%vort(1:pmlib_nvort,i,j,k)

          IF( PRESENT(dvort) .AND. dvort )THEN
            part%dvort(1:pmlib_nvort,inp) = mesh%dvort(1:pmlib_nvort,i,j,k)
          ELSE
            part%dvort(1:pmlib_nvort,inp) = 0.0_MK
          END IF
          part%dvort_rk(1:pmlib_nvort,inp) = 0.0_MK

        END IF

      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_remesh
