!---------------------------------------------------------------------------------!
! flowcases_knot.f
!---------------------------------------------------------------------------------!
! Initiates the vorticity field of a vortex knot by defining a vortex 
! filament as a parametric curve which is regularised with the given 
! vortex core radius
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_knot(patch,topo_all,mesh,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_interpolation

USE poisson_solver_module

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                            :: patch
  TYPE(class_topology_all)                     :: topo_all
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                 :: PI
  INTEGER                                  :: i,j,k
  REAL(MK)                                 :: x1,x2,x3
  REAL(MK)                                 :: y1,y2,y3
  REAL(MK)                                 :: z1,z2,z3

  REAL(MK)                                 :: theta, px, py, pz
  REAL(MK)                                 :: theta1, px1, py1, pz1
  REAL(MK)                                 :: theta2, px2, py2, pz2
  REAL(MK)                                 :: len_segment

  INTEGER                                  :: npart,inp

  REAL(MK),DIMENSION(ndim)                 :: xmin, xmax, dx
  INTEGER,DIMENSION(ndim)                  :: ncell, offset
  INTEGER,DIMENSION(2*ndim)                :: nghost

  REAL(MK)                                 :: r,r0
  REAL(MK)                                 :: circ
  REAL(MK)                                 :: norm

  REAL(MK)                                 :: sigma

  REAL(MK),DIMENSION(:,:),POINTER          :: part_pos, part_circ

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  r0       = 0.5_MK
  sigma    = 0.15_MK*r0
  circ     = 1.0_MK
  vort_max = circ/(2.0_MK*PI*sigma**2)

  dx     = topo_all%cuboid(rank)%dx
  xmin   = topo_all%cuboid(rank)%xmin
  xmax   = topo_all%cuboid(rank)%xmax
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel = 0.0_MK

!---------------------------------------------------------------------------------!
! Vortex filament by a parametric curve
!---------------------------------------------------------------------------------!
! Initiate and allocate circulation particles on processors
  npart = 3000
  inp   = 0
  DO i = 1,npart
    theta = ( REAL(i,MK) - 0.5_MK )/REAL(npart,MK) *(2.0_MK*PI)

    px = - SIN(3.0_MK * theta)/3.0_MK
    py = ( COS(theta) - 2.0_MK * COS(2.0_MK * theta) )/3.0_MK
    pz = ( SIN(theta) + 2.0_MK * SIN(2.0_MK * theta) )/3.0_MK

    IF( px .GE. xmin(1) .AND. px .LT. xmax(1) .AND. & 
      & py .GE. xmin(2) .AND. py .LT. xmax(2) .AND. & 
      & pz .GE. xmin(3) .AND. pz .LT. xmax(3)  )THEN
      
      inp = inp + 1
    END IF
  END DO

  ALLOCATE( part_pos(ndim,inp), part_circ(nvort,inp) )

! Initiate circulation particles
  inp   = 0
  DO i = 1,npart
    theta  = ( REAL(i,MK) - 0.5_MK )/REAL(npart,MK) *(2.0_MK*PI)
    theta1 = ( REAL(i,MK) - 1.0_MK )/REAL(npart,MK) *(2.0_MK*PI)
    theta2 = ( REAL(i,MK) )/REAL(npart,MK) *(2.0_MK*PI)

    px = - SIN(3.0_MK * theta)/3.0_MK
    py = ( COS(theta) - 2.0_MK * COS(2.0_MK * theta) )/3.0_MK
    pz = ( SIN(theta) + 2.0_MK * SIN(2.0_MK * theta) )/3.0_MK

    px1 = - SIN(3.0_MK * theta1)/3.0_MK
    py1 = ( COS(theta1) - 2.0_MK * COS(2.0_MK * theta1) )/3.0_MK
    pz1 = ( SIN(theta1) + 2.0_MK * SIN(2.0_MK * theta1) )/3.0_MK

    px2 = - SIN(3.0_MK * theta2)/3.0_MK
    py2 = ( COS(theta2) - 2.0_MK * COS(2.0_MK * theta2) )/3.0_MK
    pz2 = ( SIN(theta2) + 2.0_MK * SIN(2.0_MK * theta2) )/3.0_MK

    IF( px .GE. xmin(1) .AND. px .LT. xmax(1) .AND. & 
      & py .GE. xmin(2) .AND. py .LT. xmax(2) .AND. & 
      & pz .GE. xmin(3) .AND. pz .LT. xmax(3)  )THEN
      inp = inp + 1

      part_pos(1,inp) = px
      part_pos(2,inp) = py
      part_pos(3,inp) = pz

      len_segment = SQRT( (px2 - px1)**2 + (py2 - py1)**2 + (pz2 - pz1)**2 ) 

      part_circ(1,inp) = - COS(3.0_MK * theta)
      part_circ(2,inp) = -( SIN(theta) - 4.0_MK*SIN(2.0_MK*theta) )/3.0_MK
      part_circ(3,inp) =  ( COS(theta) + 4.0_MK*COS(2.0_MK*theta) )/3.0_MK

      norm = SQRT( part_circ(1,inp)**2 + &
                 & part_circ(2,inp)**2 + &
                 & part_circ(3,inp)**2) * dx(1)*dx(2)*dx(3)

      part_circ(1,inp) = len_segment * circ * part_circ(1,inp)/norm
      part_circ(2,inp) = len_segment * circ * part_circ(2,inp)/norm
      part_circ(3,inp) = len_segment * circ * part_circ(3,inp)/norm

    END IF
  END DO

  topo_all%cuboid(rank)%npart = inp

! Interpolate to mesh
  CALL pmlib_interp_particle_mesh( topo_all%cuboid, part_pos, part_circ, &
                                 & mesh%vort, ierr, clear = .TRUE.)

! Regularise mesh to get Gauss distribution of radius sigma
  offset = (/ 1, 1, 1 /)
  CALL poisson_solver_push(1,offset,mesh%vort )
  CALL poisson_solver_smooth3d(1,sigma)
  CALL poisson_solver_pull(1,offset,mesh%vel,mesh%vort)

!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .FALSE.

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN


END SUBROUTINE flowcases_knot

