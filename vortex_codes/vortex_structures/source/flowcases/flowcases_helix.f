!---------------------------------------------------------------------------------!
! flowcases_helix.f
!---------------------------------------------------------------------------------!
! Initiates the vorticity field of a helix vortex by defining a vortex 
! filament as a parametric curve which is regularised with the given 
! vortex core radius
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_helix(patch,topo_all,mesh,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_interpolation
USE pmlib_mod_regularise

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
  REAL(MK)                                 :: len_segment, len_domain

  INTEGER                                  :: npart,inp

  REAL(MK),DIMENSION(ndim)                 :: xmin, xmax, dx
  INTEGER,DIMENSION(ndim)                  :: ncell
  INTEGER,DIMENSION(2*ndim)                :: nghost

  REAL(MK)                                 :: r,r0
  REAL(MK)                                 :: circ
  REAL(MK)                                 :: nrev
  REAL(MK)                                 :: norm
  REAL(MK)                                 :: x_pert_amp, x_pert_freq
  REAL(MK)                                 :: r_pert_amp, r_pert_freq

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
  nrev     = 4.0_MK
  r0       = 0.5_MK
  sigma    = 0.15_MK*r0
  circ     = 1.0_MK
  vort_max = circ/(2.0_MK*PI*sigma**2)

  x_pert_amp  = 0.01_MK*r0
  x_pert_freq = 4.0_MK

  r_pert_amp  = 0.00_MK*r0
  r_pert_freq = 0.0_MK

  len_domain = patch%xmax(1) - patch%xmin(1)

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

    px = len_domain*theta/(2.0_MK*PI) + patch%xmin(1) + & 
       & x_pert_amp * SIN(nrev * x_pert_freq * theta)
    py = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta) ) &
       & * r0 * COS(nrev*theta)
    pz = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta) ) & 
       & * r0 * SIN(nrev*theta)

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

    px = len_domain*theta/(2.0_MK*PI) + patch%xmin(1) + & 
       & x_pert_amp * SIN(nrev * x_pert_freq * theta)
    py = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta) ) &
       & * r0 * COS(nrev*theta)
    pz = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta) ) & 
       & * r0 * SIN(nrev*theta)

    px1 = len_domain*theta1/(2.0_MK*PI) + patch%xmin(1) + & 
       & x_pert_amp * SIN(nrev * x_pert_freq * theta1)
    py1 = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta1) ) &
       & * r0 * COS(nrev*theta1)
    pz1 = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta1) ) & 
       & * r0 * SIN(nrev*theta1)

    px2 = len_domain*theta2/(2.0_MK*PI) + patch%xmin(1) + & 
       & x_pert_amp * SIN(nrev * x_pert_freq * theta2)
    py2 = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta2) ) &
       & * r0 * COS(nrev*theta2)
    pz2 = ( 1.0_MK +  r_pert_amp * SIN(nrev*r_pert_freq*theta2) ) & 
       & * r0 * SIN(nrev*theta2)

    IF( px .GE. xmin(1) .AND. px .LT. xmax(1) .AND. & 
      & py .GE. xmin(2) .AND. py .LT. xmax(2) .AND. & 
      & pz .GE. xmin(3) .AND. pz .LT. xmax(3)  )THEN
      inp = inp + 1

      part_pos(1,inp) = px
      part_pos(2,inp) = py
      part_pos(3,inp) = pz

      len_segment = SQRT( (px2 - px1)**2 + (py2 - py1)**2 + (pz2 - pz1)**2 ) 

      part_circ(1,inp) = len_domain/(2.0_MK*PI) + x_pert_amp * nrev & 
                       & * x_pert_freq * COS(nrev*x_pert_freq *theta)
      part_circ(2,inp) = r_pert_amp * r_pert_freq * nrev * r0 &
                       & * COS(nrev*r_pert_freq*theta) &
                       & * COS(nrev * theta) - r0 * nrev & 
                       & * (1.0_MK + r_pert_amp & 
                       & * SIN(nrev*r_pert_freq*theta) ) &
                       & * SIN(nrev * theta)
      part_circ(3,inp) = r_pert_amp * r_pert_freq * nrev * r0 &
                       & * COS(nrev*r_pert_freq*theta) &
                       & * SIN(nrev * theta) + r0 * nrev &
                       & * (1.0_MK + r_pert_amp &
                       & * SIN(nrev*r_pert_freq*theta) ) &
                       & * COS(nrev * theta);

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
  CALL pmlib_regularise(patch,topo_all,mesh%vort,2,sigma,ierr)

!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .FALSE.

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE flowcases_helix

