!---------------------------------------------------------------------------------!
! flowcases_ringfit.f
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_ringfit(topo,mesh,ierr)

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
  REAL(MK),DIMENSION(0:4,1:5)                  :: c
  REAL(MK),DIMENSION(0:4,1:5,1:3)              :: b 

  REAL(MK)                                     :: rho, phi, theta, psi
  REAL(MK)                                     :: sigma0, sigma
  REAL(MK)                                     :: circ, vort_mag, norm
  REAL(MK)                                     :: sum_vort, vort_sum, max_vort

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  r0       = 0.5_MK
  sigma0   = 0.2_MK*r0
  circ     = 1.0_MK

  b = 0.0_MK

  b(0,1,1) =    9.2_MK; b(0,1,2) =   -164_MK; b(0,1,3) =   338.5_MK;
  b(0,2,1) =  -10.0_MK; b(0,2,2) =  586.8_MK; b(0,2,3) = -1342.8_MK;
  b(0,3,1) =  -83.9_MK; b(0,3,2) = -276.5_MK; b(0,3,3) =  1040.6_MK;

  b(1,1,1) =    6.9_MK; b(1,1,2) =   17.2_MK; b(1,1,3) =   -55.6_MK;
  b(1,2,1) =   14.9_MK; b(1,2,2) = -584.2_MK; b(1,2,3) =  1328.7_MK;
  b(1,3,1) =   98.0_MK; b(1,3,2) =  307.9_MK; b(1,3,3) = -1123.8_MK;
  
  b(2,1,1) =  -11.3_MK; b(2,1,2) =   27.1_MK; b(2,1,3) =   -19.3_MK;
  b(2,2,1) =  125.5_MK; b(2,2,2) = -741.1_MK; b(2,2,3) =  1438.5_MK;
  b(2,3,1) = -372.6_MK; b(2,3,2) = 2632.5_MK; b(2,3,3) = -5493.9_MK;

  b(3,3,1) =   27.5_MK; b(3,3,2) = -114.8_MK; b(3,3,3) =   182.6_MK;
  b(3,2,1) =   -4.1_MK; b(3,2,2) =   -7.8_MK; b(3,2,3) =    43.5_MK;

  b(4,1,1) =   -0.1_MK; b(4,1,2) =    3.4_MK; b(4,1,3) =    -7.6_MK;
  b(4,2,1) =   -0.5_MK; b(4,2,2) =   -4.1_MK; b(4,2,3) =    13.6_MK;

  c(0,1) = b(0,1,1) * (sigma0/r0)**1 + b(0,1,2) * (sigma0/r0)**2 + b(0,1,3) * (sigma0/r0)**3
  c(0,2) = b(0,2,1) * (sigma0/r0)**1 + b(0,2,2) * (sigma0/r0)**2 + b(0,2,3) * (sigma0/r0)**3
  c(0,3) = b(0,3,1) * (sigma0/r0)**1 + b(0,3,2) * (sigma0/r0)**2 + b(0,3,3) * (sigma0/r0)**3
  c(0,4) = b(0,4,1) * (sigma0/r0)**1 + b(0,4,2) * (sigma0/r0)**2 + b(0,4,3) * (sigma0/r0)**3
  c(0,5) = b(0,5,1) * (sigma0/r0)**1 + b(0,5,2) * (sigma0/r0)**2 + b(0,5,3) * (sigma0/r0)**3

  c(1,1) = b(1,1,1) * (sigma0/r0)**1 + b(1,1,2) * (sigma0/r0)**2 + b(1,1,3) * (sigma0/r0)**3
  c(1,2) = b(1,2,1) * (sigma0/r0)**1 + b(1,2,2) * (sigma0/r0)**2 + b(1,2,3) * (sigma0/r0)**3
  c(1,3) = b(1,3,1) * (sigma0/r0)**1 + b(1,3,2) * (sigma0/r0)**2 + b(1,3,3) * (sigma0/r0)**3
  c(1,4) = b(1,4,1) * (sigma0/r0)**1 + b(1,4,2) * (sigma0/r0)**2 + b(1,4,3) * (sigma0/r0)**3
  c(1,5) = b(1,5,1) * (sigma0/r0)**1 + b(1,5,2) * (sigma0/r0)**2 + b(1,5,3) * (sigma0/r0)**3

  c(2,1) = b(2,1,1) * (sigma0/r0)**1 + b(2,1,2) * (sigma0/r0)**2 + b(2,1,3) * (sigma0/r0)**3
  c(2,2) = b(2,2,1) * (sigma0/r0)**1 + b(2,2,2) * (sigma0/r0)**2 + b(2,2,3) * (sigma0/r0)**3
  c(2,3) = b(2,3,1) * (sigma0/r0)**1 + b(2,3,2) * (sigma0/r0)**2 + b(2,3,3) * (sigma0/r0)**3
  c(2,4) = b(2,4,1) * (sigma0/r0)**1 + b(2,4,2) * (sigma0/r0)**2 + b(2,4,3) * (sigma0/r0)**3
  c(2,5) = b(2,5,1) * (sigma0/r0)**1 + b(2,5,2) * (sigma0/r0)**2 + b(2,5,3) * (sigma0/r0)**3

  c(3,1) = b(3,1,1) * (sigma0/r0)**1 + b(3,1,2) * (sigma0/r0)**2 + b(3,1,3) * (sigma0/r0)**3
  c(3,2) = b(3,2,1) * (sigma0/r0)**1 + b(3,2,2) * (sigma0/r0)**2 + b(3,2,3) * (sigma0/r0)**3
  c(3,3) = b(3,3,1) * (sigma0/r0)**1 + b(3,3,2) * (sigma0/r0)**2 + b(3,3,3) * (sigma0/r0)**3
  c(3,4) = b(3,4,1) * (sigma0/r0)**1 + b(3,4,2) * (sigma0/r0)**2 + b(3,4,3) * (sigma0/r0)**3
  c(3,5) = b(3,5,1) * (sigma0/r0)**1 + b(3,5,2) * (sigma0/r0)**2 + b(3,5,3) * (sigma0/r0)**3

  c(4,1) = b(4,1,1) * (sigma0/r0)**1 + b(4,1,2) * (sigma0/r0)**2 + b(4,1,3) * (sigma0/r0)**3
  c(4,2) = b(4,2,1) * (sigma0/r0)**1 + b(4,2,2) * (sigma0/r0)**2 + b(4,2,3) * (sigma0/r0)**3
  c(4,3) = b(4,3,1) * (sigma0/r0)**1 + b(4,3,2) * (sigma0/r0)**2 + b(4,3,3) * (sigma0/r0)**3
  c(4,4) = b(4,4,1) * (sigma0/r0)**1 + b(4,4,2) * (sigma0/r0)**2 + b(4,4,3) * (sigma0/r0)**3
  c(4,5) = b(4,5,1) * (sigma0/r0)**1 + b(4,5,2) * (sigma0/r0)**2 + b(4,5,3) * (sigma0/r0)**3

  dx     = topo(rank)%dx
  xmin   = topo(rank)%xmin
  ncell  = topo(rank)%ncell
  nghost = topo(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel  = 0.0_MK
  sum_vort  = 0.0_MK
  max_vort  = 0.0_MK
  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        theta = ATAN2(pz,py)

        rho = SQRT(pz*pz + py*py)

        phi = SQRT((rho - r0)**2 + px**2)
        psi = ATAN2(px,-(rho - r0))

        sigma = sigma0*( 1.0_MK  &
                       & + c(0,1) * (phi/r0)**1                      &
                       & + c(0,2) * (phi/r0)**2                      &
                       & + c(0,3) * (phi/r0)**3                      &
                       & + c(0,4) * (phi/r0)**4                      &
                       & + c(0,5) * (phi/r0)**5                      &
                       & + c(1,1) * COS(1.0_MK*psi) * (phi/r0)**1    &
                       & + c(1,2) * COS(1.0_MK*psi) * (phi/r0)**2    &
                       & + c(1,3) * COS(1.0_MK*psi) * (phi/r0)**3    &
                       & + c(1,4) * COS(1.0_MK*psi) * (phi/r0)**4    &
                       & + c(1,5) * COS(1.0_MK*psi) * (phi/r0)**5    &
                       & + c(2,1) * COS(2.0_MK*psi) * (phi/r0)**1    &
                       & + c(2,2) * COS(2.0_MK*psi) * (phi/r0)**2    &
                       & + c(2,3) * COS(2.0_MK*psi) * (phi/r0)**3    &
                       & + c(2,4) * COS(2.0_MK*psi) * (phi/r0)**4    &
                       & + c(2,5) * COS(2.0_MK*psi) * (phi/r0)**5    &
                       & + c(3,1) * COS(3.0_MK*psi) * (phi/r0)**1    &
                       & + c(3,2) * COS(3.0_MK*psi) * (phi/r0)**2    &
                       & + c(3,3) * COS(3.0_MK*psi) * (phi/r0)**3    &
                       & + c(3,4) * COS(3.0_MK*psi) * (phi/r0)**4    &
                       & + c(3,5) * COS(3.0_MK*psi) * (phi/r0)**5    &
                       & + c(4,1) * COS(4.0_MK*psi) * (phi/r0)**1    &
                       & + c(4,2) * COS(4.0_MK*psi) * (phi/r0)**2    &
                       & + c(4,3) * COS(4.0_MK*psi) * (phi/r0)**3    &
                       & + c(4,4) * COS(4.0_MK*psi) * (phi/r0)**4    &
                       & + c(4,5) * COS(4.0_MK*psi) * (phi/r0)**5 )

        IF( sigma .GT. 0.0_MK )THEN
          vort_mag = EXP(- 0.5_MK * phi**2/sigma**2)
        ELSE
          vort_mag = 0.0_MK
        END IF

!---------------------------------------------------------------------------------!
! Project the vorticity onto the Cartesian components
!---------------------------------------------------------------------------------!
        mesh%vort(1,i,j,k) =  0.0_MK
        mesh%vort(2,i,j,k) = -SIN(theta)*vort_mag
        mesh%vort(3,i,j,k) =  COS(theta)*vort_mag

!---------------------------------------------------------------------------------!
! Sum the vorticity
!---------------------------------------------------------------------------------!
        IF( k .GE. 1 .AND. k .LE. ncell(3) .AND. &
          & j .GE. 1 .AND. j .LE. ncell(2) .AND. &
          & i .GE. 1 .AND. i .LE. ncell(1) )THEN
          sum_vort = sum_vort + vort_mag/rho
          IF(vort_mag .GT. max_vort)THEN
            max_vort = vort_mag
          END IF
        END IF

      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! Normalise to specified circulation
!---------------------------------------------------------------------------------!
  CALL MPI_ALLREDUCE( max_vort,vort_max,1,mpi_prec_real, MPI_MAX,comm,ierr )
  CALL MPI_ALLREDUCE( sum_vort,vort_sum,1,mpi_prec_real, MPI_SUM,comm,ierr )

  norm = 2.0_MK * PI * circ / ( vort_sum * dx(1)*dx(2)*dx(3) )
  vort_max = norm * vort_max

  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

        mesh%vort(1,i,j,k) = norm * mesh%vort(1,i,j,k)
        mesh%vort(2,i,j,k) = norm * mesh%vort(2,i,j,k)
        mesh%vort(3,i,j,k) = norm * mesh%vort(3,i,j,k)

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

END SUBROUTINE flowcases_ringfit

