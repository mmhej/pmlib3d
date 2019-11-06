!---------------------------------------------------------------------------------!
! flowcases_taylor_green.f
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_taylor_green(topo,mesh,ierr)

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
  vort_max = 1.0_MK

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

        mesh%vort(1,i,j,k) = - vort_max * SIN(2.0_MK * PI * px) &
                           & * COS(2.0_MK * PI * py) * SIN(2.0_MK * PI * pz)
        mesh%vort(2,i,j,k) = - vort_max * COS(2.0_MK * PI * px) & 
                           & * SIN(2.0_MK * PI * py) * SIN(2.0_MK * PI * pz)
        mesh%vort(3,i,j,k) = - vort_max * COS(2.0_MK * PI * px) &
                           & * COS(2.0_MK * PI * py) * COS(2.0_MK * PI * pz)

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

END SUBROUTINE flowcases_taylor_green

