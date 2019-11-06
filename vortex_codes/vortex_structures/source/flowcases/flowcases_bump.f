!---------------------------------------------------------------------------------!
! flowcases_bump.f
!---------------------------------------------------------------------------------!
! This routine initialises the vorticity field
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_bump(topo,mesh,ierr)

USE pmlib_mod_topology
USE pmlib_mod_mesh

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER    :: topo
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                     :: PI
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  REAL(MK)                                     :: r,r0
  REAL(MK)                                     :: a,b,c

  REAL(MK)                                     :: rho, phi, theta
  REAL(MK)                                     :: vort_mag, sigma

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  r0       = 0.5_MK
  vort_max = 1.0_MK
  c        = 20.0_MK

  dx     = topo(rank)%dx
  xmin   = topo(rank)%xmin
  ncell  = topo(rank)%ncell
  nghost = topo(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel = 0.0_MK

  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)

    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)

      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        rho   = SQRT(pz*pz + py*py)
        phi   = SQRT((rho - r0)**2 + px**2)
        theta = ATAN2(pz,py)

        IF ( phi .LT. r0) THEN

          vort_mag = -EXP(- c*r0**2/(2.0_MK*r0*rho - rho**2 - px**2)) * &
                   & ( 4.0_MK*c**2*r0**4*px**2*rho**2 &
                   & - 16.0_MK*r0**4*rho**4 &
                   & + 32.0_MK*r0**3*rho**5 &
                   & - 24.0_MK*r0**2*rho**6 &
                   & + 8.0_MK*r0*rho**7 & 
                   & - 4.0_MK*rho**6*px**2 &
                   & - 6.0_MK*rho**4*px**4 & 
                   & - 4.0_MK*rho**2*px**6 &
                   & - 8.0_MK*c*r0**5*rho**3 &
                   & + 8.0_MK*c*r0**4*rho**4 &
                   & - 6.0_MK*c*r0**3*rho**5 & 
                   & + 4.0_MK*c**2*r0**6*rho**2 &
                   & - 8.0_MK*c**2*r0**5*rho**3 & 
                   & + 4.0_MK*c**2*r0**4*rho**4 & 
                   & + 2.0_MK*c*r0**2*rho**6 &
                   & + 32.0_MK*r0**3*rho**3*px**2 &
                   & - 48.0_MK*r0**2*rho**4*px**2 &
                   & - 24.0_MK*r0**2*rho**2*px**4 & 
                   & + 24.0_MK*r0*rho**5*px**2 &
                   & + 24.0_MK*r0*rho**3*px**4 & 
                   & + 8.0_MK*r0*rho*px**6 & 
                   & + 2.0_MK*c*r0**3*rho*px**4 &
                   & + 2.0_MK*c*r0**2*rho**2*px**4 &
                   & - 4.0_MK*c*r0**3*rho**3*px**2 & 
                   & + 4.0_MK*c*r0**2*rho**4*px**2 & 
                   & - rho**8 - px**8) & 
                   & * (2.0_MK*r0*rho - rho**2 - px**2)**(-4)*rho**(-2) &
                   & / ( EXP(-c)*(4.0_MK*c + 1.0_MK)/r0**2 )

            mesh%vort(1,i,j,k) =   0.0_MK 
            mesh%vort(2,i,j,k) = - SIN(theta)*vort_mag
            mesh%vort(3,i,j,k) =   COS(theta)*vort_mag

          ELSE

            mesh%vort(1,i,j,k) = 0.0_MK
            mesh%vort(2,i,j,k) = 0.0_MK
            mesh%vort(3,i,j,k) = 0.0_MK

          ENDIF

      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .FALSE.

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE flowcases_bump

