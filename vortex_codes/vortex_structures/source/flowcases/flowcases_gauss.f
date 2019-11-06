!---------------------------------------------------------------------------------!
! flowcases_gauss.f
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_gauss(topo,mesh,ierr,perturbe)

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

  REAL(MK)                                     :: r,r0
  REAL(MK)                                     :: a,b,c

  REAL(MK)                                     :: rho, phi, theta, psi
  REAL(MK)                                     :: sigma
  REAL(MK)                                     :: circ, vort_mag

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
!  r0   = 0.535_MK
!  sigma0 = 0.0673547
!  vort_max = 33.228986031609693

  r0       = 0.5_MK
  sigma    = 0.2_MK*r0
  circ     = 1.0_MK
  vort_max = circ/(2.0_MK*PI*sigma**2)

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

        theta = ATAN2(pz,py)

        rho = SQRT(pz*pz + py*py)
        phi = SQRT((rho - r0)**2 + px**2)

        vort_mag = circ/(2.0_MK*PI*sigma**2) * EXP(- 0.5_MK * phi**2/sigma**2)


        mesh%vort(1,i,j,k) =  0.0_MK
        mesh%vort(2,i,j,k) = -SIN(theta)*vort_mag
        mesh%vort(3,i,j,k) =  COS(theta)*vort_mag

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

END SUBROUTINE flowcases_gauss

