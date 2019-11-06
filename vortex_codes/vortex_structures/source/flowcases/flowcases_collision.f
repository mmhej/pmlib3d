!---------------------------------------------------------------------------------!
! flowcases_collision.f
!---------------------------------------------------------------------------------!
! This routine initialises the vorticity field
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_collision(topo,mesh,ierr,perturbe)

USE pmlib_mod_topology
USE pmlib_mod_mesh

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER    :: topo
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr
  REAL(MK), OPTIONAL, INTENT(IN)               :: perturbe

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                     :: PI
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  REAL(MK)                                     :: r,r0,a

  REAL(MK)                                     :: rho, phi1, phi2, theta
  REAL(MK)                                     :: sigma0, sigma1, sigma2
  REAL(MK)                                     :: vort_mag
  REAL(MK)                                     :: perturbe1, perturbe2

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  r0   = 0.5_MK
  sigma0 = 0.2_MK*r0
  vort_max = 1.0_MK/(2.0_MK*PI*sigma0**2)

  dx     = topo(rank)%dx
  xmin   = topo(rank)%xmin
  ncell  = topo(rank)%ncell
  nghost = topo(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel  = 0.0_MK

  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)

    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)

      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        theta = ATAN2(pz,py)

        IF( PRESENT(perturbe) )THEN
          sigma1 = (1.0_MK + perturbe*SIN(10.0_MK*theta + 26.954_MK))*sigma0
          sigma2 = (1.0_MK + perturbe*SIN(14.0_MK*theta + 14.132_MK))*sigma0
        ELSE
          sigma1 = sigma0
          sigma2 = sigma0
        END IF

        rho     = SQRT(pz*pz + py*py)
        phi1    = SQRT((rho - r0)**2 + (px - r0)**2)
        phi2    = SQRT((rho - r0)**2 + (px + r0)**2)

        IF ( phi1 .LT. r0) THEN
          vort_mag = 1.0/(2.0_MK*PI*sigma1**2) * EXP(- phi1**2/sigma1**2)

          mesh%vort(1,i,j,k) =  0.0_MK 
          mesh%vort(2,i,j,k) =  SIN(theta)*vort_mag
          mesh%vort(3,i,j,k) = -COS(theta)*vort_mag
        END IF

        IF ( phi2 .LT. r0) THEN
          vort_mag = 1.0/(2.0_MK*PI*sigma2**2) * EXP(- phi2**2/sigma2**2)

          mesh%vort(1,i,j,k) =  0.0_MK 
          mesh%vort(2,i,j,k) = -SIN(theta)*vort_mag
          mesh%vort(3,i,j,k) =  COS(theta)*vort_mag
        END IF

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

END SUBROUTINE flowcases_collision

