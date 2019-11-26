!---------------------------------------------------------------------------------!
! flowcases_ring.f
!---------------------------------------------------------------------------------!
! Initiates the vorticity field of a vortex ring by defining a vortex 
! filament as a parametric curve which is regularised with the given 
! vortex core radius
!---------------------------------------------------------------------------------!
SUBROUTINE flowcases_ring(patch,topo_all,mesh,part,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_particles
USE pmlib_mod_interpolation
!USE pmlib_mod_regularise

USE poisson_solver_module

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                            :: patch
  TYPE(class_topology_all)                     :: topo_all
  TYPE(class_mesh)                             :: mesh
  TYPE(class_particles)                        :: part
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                 :: PI
  INTEGER                                  :: i,j,k
  REAL(MK)                                 :: x1,x2,x3
  REAL(MK)                                 :: y1,y2,y3
  REAL(MK)                                 :: z1,z2,z3

  REAL(MK)                                 ::  px,  py,  pz
  REAL(MK)                                 :: px1, py1, pz1
  REAL(MK)                                 :: px2, py2, pz2
  REAL(MK)                                 :: len_segment

  INTEGER                                  :: npart,inp

  REAL(MK),DIMENSION(ndim)                 :: xmin, xmax, dx
  INTEGER,DIMENSION(ndim)                  :: ncell, offset
  INTEGER,DIMENSION(2*ndim)                :: nghost

  REAL(MK)                                 :: r,r0
  REAL(MK)                                 :: sigma
  REAL(MK)                                 :: circ
  REAL(MK)                                 :: norm
  REAL(MK)                                 :: pert_x, pert_r

  REAL(MK), DIMENSION(4)                   :: u
  REAL(MK)                                 :: a,b,c,d
  REAL(MK)                                 :: theta, theta1, theta2
  REAL(MK),DIMENSION(:),POINTER            ::  gr,  dgr,  gx,  dgx 
  REAL(MK),DIMENSION(:),POINTER            :: gr1, gx1, gr2, gx2
!  REAL(MK),DIMENSION(:,:),POINTER          :: part_pos, part_circ

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  r0       = 0.5_MK
  sigma    = 0.2_MK*r0
  circ     = 1.0_MK
  vort_max = circ/(2.0_MK*PI*sigma**2)

  pert_x  = 0.0000_MK
  pert_r  = 0.0000_MK

  dx     = topo_all%cuboid(rank)%dx
  xmin   = topo_all%cuboid(rank)%xmin
  xmax   = topo_all%cuboid(rank)%xmax
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

  mesh%vort = 0.0_MK
  mesh%vel  = 0.0_MK

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
    CALL pmlib_write( rank,'flowcases_ring', &
                    & 'Failed to deallocate existing particles.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Create vortex filament by a parametric curve
!---------------------------------------------------------------------------------!
  npart = 1000

  ALLOCATE(  gr(npart), dgr(npart) )
  ALLOCATE( gr1(npart), gr2(npart) )

  ALLOCATE(  gx(npart), dgx(npart) )
  ALLOCATE( gx1(npart), gx2(npart) )

  gr  = 0.0_MK
  dgr = 0.0_MK
  gx  = 0.0_MK
  dgx = 0.0_MK

! Create perturbation functions
!  CALL RANDOM_SEED

  DO j = 1,32
    CALL RANDOM_NUMBER(u)
    u(1) = u(1) - 0.5_MK
    u(2) = u(2) - 0.5_MK
    u(3) = u(3) - 0.5_MK
    u(4) = u(4) - 0.5_MK

    a = u(1)/SQRT(u(1)**2 + u(2)**2)
    b = u(2)/SQRT(u(1)**2 + u(2)**2)
    c = u(3)/SQRT(u(3)**2 + u(4)**2)
    d = u(4)/SQRT(u(3)**2 + u(4)**2)

    DO i = 1,npart
      theta  = 2.0_MK * PI * ( REAL(i,MK) - 0.5_MK )/REAL(npart,MK)
      theta1 = 2.0_MK * PI * ( REAL(i,MK) - 1.0_MK )/REAL(npart,MK)
      theta2 = 2.0_MK * PI * ( REAL(i,MK) )/REAL(npart,MK)


      gr(i)   = gr(i) + a * SIN( REAL(j,MK) * theta ) &
                    & + b * COS( REAL(j,MK) * theta )
      dgr(i)  = dgr(i) + REAL(j,MK) * a * COS( REAL(j,MK) * theta ) &
                       - REAL(j,MK) * b * SIN( REAL(j,MK) * theta )

      gr1(i)  = gr1(i) + a * SIN( REAL(j,MK) * theta1 ) &
                     & + b * COS( REAL(j,MK) * theta1 )
      gr2(i)  = gr2(i) + a * SIN( REAL(j,MK) * theta2 ) &
                     & + b * COS( REAL(j,MK) * theta2 )


      gx(i)   = gx(i) + c * SIN( REAL(j,MK) * theta ) &
                    & + d * COS( REAL(j,MK) * theta )

      dgx(i)  = dgx(i) + REAL(j,MK) * c * COS( REAL(j,MK) * theta ) &
                       - REAL(j,MK) * d * SIN( REAL(j,MK) * theta )

      gx1(i)  = gx1(i) + c * SIN( REAL(j,MK) * theta1 ) &
                     & + d * COS( REAL(j,MK) * theta1 )
      gx2(i)  = gx2(i) + c * SIN( REAL(j,MK) * theta2 ) &
                     & + d * COS( REAL(j,MK) * theta2 )
    END DO
  END DO

! Initiate and allocate circulation particles on processors
  inp   = 0
  DO i = 1,npart

    theta  = 2.0_MK * PI * ( REAL(i,MK) - 0.5_MK )/REAL(npart,MK)

    px = pert_x * gx(i)
    py = ( 1.0_MK + pert_r * gr(i) ) * r0 * COS(theta)
    pz = ( 1.0_MK + pert_r * gr(i) ) * r0 * SIN(theta)

    IF( px .GE. xmin(1) .AND. px .LT. xmax(1) .AND. & 
      & py .GE. xmin(2) .AND. py .LT. xmax(2) .AND. & 
      & pz .GE. xmin(3) .AND. pz .LT. xmax(3)  )THEN
      
      inp = inp + 1
    END IF
  END DO


!---------------------------------------------------------------------------------!
! Allocate particles
!---------------------------------------------------------------------------------!
  ALLOCATE( part%pos(ndim,inp), part%vort(nvort,inp), &
          & part%vel(ndim,inp), part%dvort(nvort,inp) )

  ALLOCATE( part%vel_rk(ndim,inp), &
          & part%dvort_rk(nvort,inp) )


! Initiate circulation particles
  inp   = 0
  DO i = 1,npart

    theta  = 2.0_MK * PI * ( REAL(i,MK) - 0.5_MK )/REAL(npart,MK)
    theta1 = 2.0_MK * PI * ( REAL(i,MK) - 1.0_MK )/REAL(npart,MK)
    theta2 = 2.0_MK * PI * ( REAL(i,MK) )/REAL(npart,MK)

    px = pert_x * gx(i)
    py = ( 1.0_MK +  pert_r * gr(i) ) * r0 * COS(theta)
    pz = ( 1.0_MK +  pert_r * gr(i) ) * r0 * SIN(theta)

    px1 = pert_x * gx1(i)
    py1 = ( 1.0_MK +  pert_r * gr1(i) ) * r0 * COS(theta1)
    pz1 = ( 1.0_MK +  pert_r * gr1(i) ) * r0 * SIN(theta1)

    px2 = pert_x * gx2(i)
    py2 = ( 1.0_MK +  pert_r * gr2(i) ) * r0 * COS(theta2)
    pz2 = ( 1.0_MK +  pert_r * gr2(i) ) * r0 * SIN(theta2)

    IF( px .GE. xmin(1) .AND. px .LT. xmax(1) .AND. & 
      & py .GE. xmin(2) .AND. py .LT. xmax(2) .AND. & 
      & pz .GE. xmin(3) .AND. pz .LT. xmax(3)  )THEN
      inp = inp + 1

      part%pos(1,inp) = px
      part%pos(2,inp) = py
      part%pos(3,inp) = pz

      len_segment = SQRT( (px2 - px1)**2 + (py2 - py1)**2 + (pz2 - pz1)**2 ) 

      part%vort(1,inp) = pert_x * dgx(i)

      part%vort(2,inp) = pert_r * r0 * dgr(i) * COS(theta) & 
                     & - r0 * (1 + pert_r * gr(i)) * SIN(theta)

      part%vort(3,inp) = pert_r * r0 * dgr(i) * SIN(theta) &
                     & + r0 * (1 + pert_r * gr(i)) * COS(theta)

      norm = SQRT( part%vort(1,inp)**2 + &
                 & part%vort(2,inp)**2 + &
                 & part%vort(3,inp)**2) * dx(1)*dx(2)*dx(3)

      part%vort(1,inp) = len_segment * circ * part%vort(1,inp)/norm
      part%vort(2,inp) = len_segment * circ * part%vort(2,inp)/norm
      part%vort(3,inp) = len_segment * circ * part%vort(3,inp)/norm

      part%dvort_rk(1,inp) = 0.0_MK
      part%dvort_rk(2,inp) = 0.0_MK
      part%dvort_rk(3,inp) = 0.0_MK

      part%vel_rk(1,inp) = 0.0_MK
      part%vel_rk(2,inp) = 0.0_MK
      part%vel_rk(3,inp) = 0.0_MK

    END IF
  END DO

! open(10,file='./plot/filament.dat')
! do i = 1,inp
!   write(10,'(6E20.12)')part_pos(:,i),part_circ(:,i)
! end do
! close(10)

  topo_all%cuboid(rank)%npart = inp

! Interpolate to mesh
  CALL pmlib_interp_particle_mesh( topo_all%cuboid, part%pos, part%vort, &
                                 & mesh%vort, ierr, clear = .TRUE.)


! Regularise mesh to get Gauss distribution of radius sigma
! CALL pmlib_regularise(patch,topo_all,mesh%vort,2,sigma,ierr)
  offset = (/ 1, 1, 1 /)
  CALL poisson_solver_push( mesh%vort, offset )
  CALL poisson_solver_smooth3d( sigma )
  CALL poisson_solver_pull( mesh%vort , mesh%vel, offset )


!---------------------------------------------------------------------------------!
! Toggle penalisation
!---------------------------------------------------------------------------------!
  penalisation = .FALSE.

!---------------------------------------------------------------------------------!
! De-allocate local pointers
!---------------------------------------------------------------------------------!
  DEALLOCATE(  gr, dgr )
  DEALLOCATE( gr1, gr2 )

  DEALLOCATE(  gx, dgx )
  DEALLOCATE( gx1, gx2 )

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE flowcases_ring

