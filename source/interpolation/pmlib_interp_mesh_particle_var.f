!---------------------------------------------------------------------------------!
! pmlib_interp_mesh_particle_var.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_interp_mesh_particle_var(topo,xpart,mesh,part,ierr,clear)

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
!  Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN) :: topo
  REAL(MK) , DIMENSION(:,:), POINTER, INTENT(IN)        :: xpart 
  REAL(MK), DIMENSION(:,:,:,:), POINTER, INTENT(IN)     :: mesh
  REAL(MK), DIMENSION(:,:), POINTER, INTENT(OUT)        :: part
  INTEGER, INTENT(OUT)                                  :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                         :: clear

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=26)                      :: caller = 'pmlib_interp_mesh_particle'
  INTEGER                                :: i, nvar
  REAL(MK),DIMENSION(pmlib_ndim)         :: xmin, dx
  INTEGER,DIMENSION(pmlib_ndim)          :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)        :: nghost

!---------------------------------------------------------------------------------!
! Variables for unrolled versions
!---------------------------------------------------------------------------------!
  REAL(MK) :: c_1_dx, c_1_dy, c_1_dz, c1, c2, c3
  REAL(MK) :: dx1,dx2,dx3,dx4,dx5,dx6
  REAL(MK) :: dy1,dy2,dy3,dy4,dy5,dy6
  REAL(MK) :: dz1,dz2,dz3,dz4,dz5,dz6
  REAL(MK) :: ax1,ax2,ax3,ax4,ax5,ax6
  REAL(MK) :: ay1,ay2,ay3,ay4,ay5,ay6
  REAL(MK) :: az1,az2,az3,az4,az5,az6
  INTEGER  :: i1,i2,i3,i4,i5,i6
  INTEGER  :: j1,j2,j3,j4,j5,j6
  INTEGER  :: k1,k2,k3,k4,k5,k6

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  nvar = SIZE(mesh,1)
  IF(nvar .NE. SIZE(part,1))THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
         & 'Dimension of interpolation variables does not correspond.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Obtain topology
!---------------------------------------------------------------------------------!
  dx     = topo(mpi_rank)%dx
  xmin   = topo(mpi_rank)%xmin
  nghost = topo(mpi_rank)%nghost
  ncell  = topo(mpi_rank)%ncell

!---------------------------------------------------------------------------------!
! Clear particles
!---------------------------------------------------------------------------------!
  IF(PRESENT(clear) .AND. clear)THEN
    part = 0.0_MK
  END IF

  IF(pmlib_interpolation_order .EQ. 3)THEN
!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
  c_1_dx = 1.0_MK/dx(1)
  c_1_dy = 1.0_MK/dx(2)
  c_1_dz = 1.0_MK/dx(3)

!---------------------------------------------------------------------------------!
! Interpolate mesh to particles
!---------------------------------------------------------------------------------!
  DO i = 1,topo(mpi_rank)%npart

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
    i2 = NINT( (xpart(1,i)-xmin(1)) * c_1_dx )
    j2 = NINT( (xpart(2,i)-xmin(2)) * c_1_dy )
    k2 = NINT( (xpart(3,i)-xmin(3)) * c_1_dz )

!---------------------------------------------------------------------------------!
! Find index of other cells
!---------------------------------------------------------------------------------!
    i1 = i2 - 1
    i3 = i2 + 1
    i4 = i2 + 2

    j1 = j2 - 1
    j3 = j2 + 1
    j4 = j2 + 2

    k1 = k2 - 1
    k3 = k2 + 1
    k4 = k2 + 2

!---------------------------------------------------------------------------------!
! Find normalised distance to other cells
!---------------------------------------------------------------------------------!
    dx2 = ( xpart(1,i) - xmin(1) - dx(1)*(REAL(i2,MK) - 0.5_MK) ) * c_1_dx
    dy2 = ( xpart(2,i) - xmin(2) - dx(2)*(REAL(j2,MK) - 0.5_MK) ) * c_1_dy
    dz2 = ( xpart(3,i) - xmin(3) - dx(3)*(REAL(k2,MK) - 0.5_MK) ) * c_1_dz

    dx1 = dx2 + 1.0_MK
    dx3 = 1.0_MK - dx2 
    dx4 = 2.0_MK - dx2 

    dy1 = dy2 + 1.0_MK
    dy3 = 1.0_MK - dy2 
    dy4 = 2.0_MK - dy2 

    dz1 = dz2 + 1.0_MK
    dz3 = 1.0_MK - dz2 
    dz4 = 2.0_MK - dz2 

!---------------------------------------------------------------------------------!
! The M-prime-4 kernel
!---------------------------------------------------------------------------------!
    ax1 = 0.5_MK * (2.0_MK - dx1)**2 * (1.0_MK - dx1)
    ax2 = 1.0_MK - 2.5_MK*dx2**2 + 1.5_MK*dx2**3
    ax3 = 1.0_MK - 2.5_MK*dx3**2 + 1.5_MK*dx3**3
    ax4 = 0.5_MK * (2.0_MK - dx4)**2 * (1.0_MK - dx4)

    ay1 = 0.5_MK * (2.0_MK - dy1)**2 * (1.0_MK - dy1)
    ay2 = 1.0_MK - 2.5_MK*dy2**2 + 1.5_MK*dy2**3
    ay3 = 1.0_MK - 2.5_MK*dy3**2 + 1.5_MK*dy3**3
    ay4 = 0.5_MK * (2.0_MK - dy4)**2 * (1.0_MK - dy4)

    az1 = 0.5_MK * (2.0_MK - dz1)**2 * (1.0_MK - dz1)
    az2 = 1.0_MK - 2.5_MK*dz2**2 + 1.5_MK*dz2**3
    az3 = 1.0_MK - 2.5_MK*dz3**2 + 1.5_MK*dz3**3
    az4 = 0.5_MK * (2.0_MK - dz4)**2 * (1.0_MK - dz4)

!---------------------------------------------------------------------------------!
! Cancel stencils that are not supported by the mesh
!---------------------------------------------------------------------------------!
    IF( i1 .LT. 1-nghost(1) ) ax1 = 0.0_MK
    IF( i1 .LT. 1-nghost(1) ) i1  = 1-nghost(1)
    IF( i2 .LT. 1-nghost(1) ) ax2 = 0.0_MK
    IF( i2 .LT. 1-nghost(1) ) i2  = 1-nghost(1)
    IF( i3 .LT. 1-nghost(1) ) ax3 = 0.0_MK
    IF( i3 .LT. 1-nghost(1) ) i3  = 1-nghost(1)
    IF( i4 .LT. 1-nghost(1) ) ax4 = 0.0_MK
    IF( i4 .LT. 1-nghost(1) ) i4  = 1-nghost(1)

    IF( i1 .GT. ncell(1)+nghost(2) ) ax1 = 0.0_MK
    IF( i1 .GT. ncell(1)+nghost(2) ) i1  = ncell(1)+nghost(2)
    IF( i2 .GT. ncell(1)+nghost(2) ) ax2 = 0.0_MK
    IF( i2 .GT. ncell(1)+nghost(2) ) i2  = ncell(1)+nghost(2)
    IF( i3 .GT. ncell(1)+nghost(2) ) ax3 = 0.0_MK
    IF( i3 .GT. ncell(1)+nghost(2) ) i3  = ncell(1)+nghost(2)
    IF( i4 .GT. ncell(1)+nghost(2) ) ax4 = 0.0_MK
    IF( i4 .GT. ncell(1)+nghost(2) ) i4  = ncell(1)+nghost(2)

    IF( j1 .LT. 1-nghost(3) ) ay1 = 0.0_MK
    IF( j1 .LT. 1-nghost(3) ) j1  = 1-nghost(3)
    IF( j2 .LT. 1-nghost(3) ) ay2 = 0.0_MK
    IF( j2 .LT. 1-nghost(3) ) j2  = 1-nghost(3)
    IF( j3 .LT. 1-nghost(3) ) ay3 = 0.0_MK
    IF( j3 .LT. 1-nghost(3) ) j3  = 1-nghost(3)
    IF( j4 .LT. 1-nghost(3) ) ay4 = 0.0_MK
    IF( j4 .LT. 1-nghost(3) ) j4  = 1-nghost(3)

    IF( j1 .GT. ncell(2)+nghost(4) ) ay1 = 0.0_MK
    IF( j1 .GT. ncell(2)+nghost(4) ) j1  = ncell(2)+nghost(4)
    IF( j2 .GT. ncell(2)+nghost(4) ) ay2 = 0.0_MK
    IF( j2 .GT. ncell(2)+nghost(4) ) j2  = ncell(2)+nghost(4)
    IF( j3 .GT. ncell(2)+nghost(4) ) ay3 = 0.0_MK
    IF( j3 .GT. ncell(2)+nghost(4) ) j3  = ncell(2)+nghost(4)
    IF( j4 .GT. ncell(2)+nghost(4) ) ay4 = 0.0_MK
    IF( j4 .GT. ncell(2)+nghost(4) ) j4  = ncell(2)+nghost(4)

    IF( k1 .LT. 1-nghost(5) ) az1 = 0.0_MK
    IF( k1 .LT. 1-nghost(5) ) k1  = 1-nghost(5)
    IF( k2 .LT. 1-nghost(5) ) az2 = 0.0_MK
    IF( k2 .LT. 1-nghost(5) ) k2  = 1-nghost(5)
    IF( k3 .LT. 1-nghost(5) ) az3 = 0.0_MK
    IF( k3 .LT. 1-nghost(5) ) k3  = 1-nghost(5)
    IF( k4 .LT. 1-nghost(5) ) az4 = 0.0_MK
    IF( k4 .LT. 1-nghost(5) ) k4  = 1-nghost(5)

    IF( k1 .GT. ncell(3)+nghost(6) ) az1 = 0.0_MK
    IF( k1 .GT. ncell(3)+nghost(6) ) k1  = ncell(3)+nghost(6)
    IF( k2 .GT. ncell(3)+nghost(6) ) az2 = 0.0_MK
    IF( k2 .GT. ncell(3)+nghost(6) ) k2  = ncell(3)+nghost(6)
    IF( k3 .GT. ncell(3)+nghost(6) ) az3 = 0.0_MK
    IF( k3 .GT. ncell(3)+nghost(6) ) k3  = ncell(3)+nghost(6)
    IF( k4 .GT. ncell(3)+nghost(6) ) az4 = 0.0_MK
    IF( k4 .GT. ncell(3)+nghost(6) ) k4  = ncell(3)+nghost(6)

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 3D kernel and apply to the field
!---------------------------------------------------------------------------------!
    part(1:nvar,i) = ax1*ay1*az1*mesh(1:nvar,i1,j1,k1) + &
                   & ax2*ay1*az1*mesh(1:nvar,i2,j1,k1) + &
                   & ax3*ay1*az1*mesh(1:nvar,i3,j1,k1) + &
                   & ax4*ay1*az1*mesh(1:nvar,i4,j1,k1) + &
                   & ax1*ay2*az1*mesh(1:nvar,i1,j2,k1) + &
                   & ax2*ay2*az1*mesh(1:nvar,i2,j2,k1) + &
                   & ax3*ay2*az1*mesh(1:nvar,i3,j2,k1) + &
                   & ax4*ay2*az1*mesh(1:nvar,i4,j2,k1) + &
                   & ax1*ay3*az1*mesh(1:nvar,i1,j3,k1) + &
                   & ax2*ay3*az1*mesh(1:nvar,i2,j3,k1) + &
                   & ax3*ay3*az1*mesh(1:nvar,i3,j3,k1) + &
                   & ax4*ay3*az1*mesh(1:nvar,i4,j3,k1) + &
                   & ax1*ay4*az1*mesh(1:nvar,i1,j4,k1) + &
                   & ax2*ay4*az1*mesh(1:nvar,i2,j4,k1) + &
                   & ax3*ay4*az1*mesh(1:nvar,i3,j4,k1) + &
                   & ax4*ay4*az1*mesh(1:nvar,i4,j4,k1) + &
                   & ax1*ay1*az2*mesh(1:nvar,i1,j1,k2) + &
                   & ax2*ay1*az2*mesh(1:nvar,i2,j1,k2) + &
                   & ax3*ay1*az2*mesh(1:nvar,i3,j1,k2) + &
                   & ax4*ay1*az2*mesh(1:nvar,i4,j1,k2) + &
                   & ax1*ay2*az2*mesh(1:nvar,i1,j2,k2) + &
                   & ax2*ay2*az2*mesh(1:nvar,i2,j2,k2) + &
                   & ax3*ay2*az2*mesh(1:nvar,i3,j2,k2) + &
                   & ax4*ay2*az2*mesh(1:nvar,i4,j2,k2) + &
                   & ax1*ay3*az2*mesh(1:nvar,i1,j3,k2) + &
                   & ax2*ay3*az2*mesh(1:nvar,i2,j3,k2) + &
                   & ax3*ay3*az2*mesh(1:nvar,i3,j3,k2) + &
                   & ax4*ay3*az2*mesh(1:nvar,i4,j3,k2) + &
                   & ax1*ay4*az2*mesh(1:nvar,i1,j4,k2) + &
                   & ax2*ay4*az2*mesh(1:nvar,i2,j4,k2) + &
                   & ax3*ay4*az2*mesh(1:nvar,i3,j4,k2) + &
                   & ax4*ay4*az2*mesh(1:nvar,i4,j4,k2) + &
                   & ax1*ay1*az3*mesh(1:nvar,i1,j1,k3) + &
                   & ax2*ay1*az3*mesh(1:nvar,i2,j1,k3) + &
                   & ax3*ay1*az3*mesh(1:nvar,i3,j1,k3) + &
                   & ax4*ay1*az3*mesh(1:nvar,i4,j1,k3) + &
                   & ax1*ay2*az3*mesh(1:nvar,i1,j2,k3) + &
                   & ax2*ay2*az3*mesh(1:nvar,i2,j2,k3) + &
                   & ax3*ay2*az3*mesh(1:nvar,i3,j2,k3) + &
                   & ax4*ay2*az3*mesh(1:nvar,i4,j2,k3) + &
                   & ax1*ay3*az3*mesh(1:nvar,i1,j3,k3) + &
                   & ax2*ay3*az3*mesh(1:nvar,i2,j3,k3) + &
                   & ax3*ay3*az3*mesh(1:nvar,i3,j3,k3) + &
                   & ax4*ay3*az3*mesh(1:nvar,i4,j3,k3) + &
                   & ax1*ay4*az3*mesh(1:nvar,i1,j4,k3) + &
                   & ax2*ay4*az3*mesh(1:nvar,i2,j4,k3) + &
                   & ax3*ay4*az3*mesh(1:nvar,i3,j4,k3) + &
                   & ax4*ay4*az3*mesh(1:nvar,i4,j4,k3) + &
                   & ax1*ay1*az4*mesh(1:nvar,i1,j1,k4) + &
                   & ax2*ay1*az4*mesh(1:nvar,i2,j1,k4) + &
                   & ax3*ay1*az4*mesh(1:nvar,i3,j1,k4) + &
                   & ax4*ay1*az4*mesh(1:nvar,i4,j1,k4) + &
                   & ax1*ay2*az4*mesh(1:nvar,i1,j2,k4) + &
                   & ax2*ay2*az4*mesh(1:nvar,i2,j2,k4) + &
                   & ax3*ay2*az4*mesh(1:nvar,i3,j2,k4) + &
                   & ax4*ay2*az4*mesh(1:nvar,i4,j2,k4) + &
                   & ax1*ay3*az4*mesh(1:nvar,i1,j3,k4) + &
                   & ax2*ay3*az4*mesh(1:nvar,i2,j3,k4) + &
                   & ax3*ay3*az4*mesh(1:nvar,i3,j3,k4) + &
                   & ax4*ay3*az4*mesh(1:nvar,i4,j3,k4) + &
                   & ax1*ay4*az4*mesh(1:nvar,i1,j4,k4) + &
                   & ax2*ay4*az4*mesh(1:nvar,i2,j4,k4) + &
                   & ax3*ay4*az4*mesh(1:nvar,i3,j4,k4) + &
                   & ax4*ay4*az4*mesh(1:nvar,i4,j4,k4)

  END DO ! particle

  ELSEIF(pmlib_interpolation_order .EQ. 4)THEN
!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
  c_1_dx = 1.0_MK/dx(1)
  c_1_dy = 1.0_MK/dx(2)
  c_1_dz = 1.0_MK/dx(3)

  c1 = -1.0_MK/24.0_MK
  c2 =  1.0_MK/24.0_MK 
  c3 = -1.0_MK/12.0_MK

!---------------------------------------------------------------------------------!
! Interpolate mesh to particles
!---------------------------------------------------------------------------------!
  DO i = 1,topo(mpi_rank)%npart

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
    i3 = NINT( (xpart(1,i)-xmin(1)) * c_1_dx )
    j3 = NINT( (xpart(2,i)-xmin(2)) * c_1_dy )
    k3 = NINT( (xpart(3,i)-xmin(3)) * c_1_dz )

!---------------------------------------------------------------------------------!
! Find index of other cells
!---------------------------------------------------------------------------------!
    i1 = i3 - 2
    i2 = i3 - 1
    i4 = i3 + 1
    i5 = i3 + 2
    i6 = i3 + 3

    j1 = j3 - 2
    j2 = j3 - 1
    j4 = j3 + 1
    j5 = j3 + 2
    j6 = j3 + 3

    k1 = k3 - 2
    k2 = k3 - 1
    k4 = k3 + 1
    k5 = k3 + 2
    k6 = k3 + 3

!---------------------------------------------------------------------------------!
! Find normalised distance to other cells
!---------------------------------------------------------------------------------!
    dx3 = ( xpart(1,i) - xmin(1) - dx(1)*(REAL(i3,MK) - 0.5_MK) ) * c_1_dx
    dy3 = ( xpart(2,i) - xmin(2) - dx(2)*(REAL(j3,MK) - 0.5_MK) ) * c_1_dy
    dz3 = ( xpart(3,i) - xmin(3) - dx(3)*(REAL(k3,MK) - 0.5_MK) ) * c_1_dz

    dx1 = dx3 + 2.0_MK
    dx2 = dx3 + 1.0_MK
    dx4 = 1.0_MK - dx3 
    dx5 = 2.0_MK - dx3 
    dx6 = 3.0_MK - dx3

    dy1 = dy3 + 2.0_MK
    dy2 = dy3 + 1.0_MK
    dy4 = 1.0_MK - dy3 
    dy5 = 2.0_MK - dy3 
    dy6 = 3.0_MK - dy3

    dz1 = dz3 + 2.0_MK
    dz2 = dz3 + 1.0_MK
    dz4 = 1.0_MK - dz3 
    dz5 = 2.0_MK - dz3 
    dz6 = 3.0_MK - dz3

!---------------------------------------------------------------------------------!
! The M-ast-6 kernel (vanRees:2011 for the correct version)
!---------------------------------------------------------------------------------!
    ax1 = c1 * (dx1 - 2.0_MK)* (dx1 - 3.0_MK)**3 & 
      & * (5.0_MK * dx1 - 8.0_MK)
    ax2 = c2 * (dx2 - 1.0_MK) * (dx2 - 2.0_MK) &
      & * (25.0_MK*dx2**3 - 114.0_MK*dx2**2 + 153.0_MK*dx2 - 48.0_MK)
    ax3 = c3 * (dx3 - 1.0_MK) * (25.0_MK*dx3**4 & 
      & - 38.0_MK*dx3**3 - 3.0_MK*dx3**2 + 12.0_MK*dx3 + 12.0_MK)
    ax4 = c3 * (dx4 - 1.0_MK) * (25.0_MK*dx4**4 & 
      & - 38.0_MK*dx4**3 - 3.0_MK*dx4**2 + 12.0_MK*dx4 + 12.0_MK)
    ax5 = c2 * (dx5 - 1.0_MK) * (dx5 - 2.0_MK) &
      & * (25.0_MK*dx5**3 - 114.0_MK*dx5**2 + 153.0_MK*dx5 - 48.0_MK)
    ax6 = c1 * (dx6 - 2.0_MK)* (dx6 - 3.0_MK)**3 & 
      & * (5.0_MK * dx6 - 8.0_MK)

    ay1 = c1 * (dy1 - 2.0_MK)* (dy1 - 3.0_MK)**3 & 
      & * (5.0_MK * dy1 - 8.0_MK)
    ay2 = c2 * (dy2 - 1.0_MK) * (dy2 - 2.0_MK) &
      & * (25.0_MK*dy2**3 - 114.0_MK*dy2**2 + 153.0_MK*dy2 - 48.0_MK)
    ay3 = c3 * (dy3 - 1.0_MK) * (25.0_MK*dy3**4 & 
      & - 38.0_MK*dy3**3 - 3.0_MK*dy3**2 + 12.0_MK*dy3 + 12.0_MK)
    ay4 = c3 * (dy4 - 1.0_MK) * (25.0_MK*dy4**4 & 
      & - 38.0_MK*dy4**3 - 3.0_MK*dy4**2 + 12.0_MK*dy4 + 12.0_MK)
    ay5 = c2 * (dy5 - 1.0_MK) * (dy5 - 2.0_MK) &
      & * (25.0_MK*dy5**3 - 114.0_MK*dy5**2 + 153.0_MK*dy5 - 48.0_MK)
    ay6 = c1 * (dy6 - 2.0_MK)* (dy6 - 3.0_MK)**3 & 
      & * (5.0_MK * dy6 - 8.0_MK)

    az1 = c1 * (dz1 - 2.0_MK)* (dz1 - 3.0_MK)**3 & 
      & * (5.0_MK * dz1 - 8.0_MK)
    az2 = c2 * (dz2 - 1.0_MK) * (dz2 - 2.0_MK) &
      & * (25.0_MK*dz2**3 - 114.0_MK*dz2**2 + 153.0_MK*dz2 - 48.0_MK)
    az3 = c3 * (dz3 - 1.0_MK) * (25.0_MK*dz3**4 & 
      & - 38.0_MK*dz3**3 - 3.0_MK*dz3**2 + 12.0_MK*dz3 + 12.0_MK)
    az4 = c3 * (dz4 - 1.0_MK) * (25.0_MK*dz4**4 & 
      & - 38.0_MK*dz4**3 - 3.0_MK*dz4**2 + 12.0_MK*dz4 + 12.0_MK)
    az5 = c2 * (dz5 - 1.0_MK) * (dz5 - 2.0_MK) &
      & * (25.0_MK*dz5**3 - 114.0_MK*dz5**2 + 153.0_MK*dz5 - 48.0_MK)
    az6 = c1 * (dz6 - 2.0_MK)* (dz6 - 3.0_MK)**3 & 
      & * (5.0_MK * dz6 - 8.0_MK)

!---------------------------------------------------------------------------------!
! Cancel stencils that are not supported by the mesh
!---------------------------------------------------------------------------------!
    IF( i1 .LT. 1-nghost(1) ) ax1 = 0.0_MK
    IF( i1 .LT. 1-nghost(1) ) i1  = 1-nghost(1)
    IF( i2 .LT. 1-nghost(1) ) ax2 = 0.0_MK
    IF( i2 .LT. 1-nghost(1) ) i2  = 1-nghost(1)
    IF( i3 .LT. 1-nghost(1) ) ax3 = 0.0_MK
    IF( i3 .LT. 1-nghost(1) ) i3  = 1-nghost(1)
    IF( i4 .LT. 1-nghost(1) ) ax4 = 0.0_MK
    IF( i4 .LT. 1-nghost(1) ) i4  = 1-nghost(1)
    IF( i5 .LT. 1-nghost(1) ) ax5 = 0.0_MK
    IF( i5 .LT. 1-nghost(1) ) i5  = 1-nghost(1)
    IF( i6 .LT. 1-nghost(1) ) ax6 = 0.0_MK
    IF( i6 .LT. 1-nghost(1) ) i6  = 1-nghost(1)

    IF( i1 .GT. ncell(1)+nghost(2) ) ax1 = 0.0_MK
    IF( i1 .GT. ncell(1)+nghost(2) ) i1  = ncell(1)+nghost(2)
    IF( i2 .GT. ncell(1)+nghost(2) ) ax2 = 0.0_MK
    IF( i2 .GT. ncell(1)+nghost(2) ) i2  = ncell(1)+nghost(2)
    IF( i3 .GT. ncell(1)+nghost(2) ) ax3 = 0.0_MK
    IF( i3 .GT. ncell(1)+nghost(2) ) i3  = ncell(1)+nghost(2)
    IF( i4 .GT. ncell(1)+nghost(2) ) ax4 = 0.0_MK
    IF( i4 .GT. ncell(1)+nghost(2) ) i4  = ncell(1)+nghost(2)
    IF( i5 .GT. ncell(1)+nghost(2) ) ax5 = 0.0_MK
    IF( i5 .GT. ncell(1)+nghost(2) ) i5  = ncell(1)+nghost(2)
    IF( i6 .GT. ncell(1)+nghost(2) ) ax6 = 0.0_MK
    IF( i6 .GT. ncell(1)+nghost(2) ) i6  = ncell(1)+nghost(2)

    IF( j1 .LT. 1-nghost(3) ) ay1 = 0.0_MK
    IF( j1 .LT. 1-nghost(3) ) j1  = 1-nghost(3)
    IF( j2 .LT. 1-nghost(3) ) ay2 = 0.0_MK
    IF( j2 .LT. 1-nghost(3) ) j2  = 1-nghost(3)
    IF( j3 .LT. 1-nghost(3) ) ay3 = 0.0_MK
    IF( j3 .LT. 1-nghost(3) ) j3  = 1-nghost(3)
    IF( j4 .LT. 1-nghost(3) ) ay4 = 0.0_MK
    IF( j4 .LT. 1-nghost(3) ) j4  = 1-nghost(3)
    IF( j5 .LT. 1-nghost(3) ) ay5 = 0.0_MK
    IF( j5 .LT. 1-nghost(3) ) j5  = 1-nghost(3)
    IF( j6 .LT. 1-nghost(3) ) ay6 = 0.0_MK
    IF( j6 .LT. 1-nghost(3) ) j6  = 1-nghost(3)

    IF( j1 .GT. ncell(2)+nghost(4) ) ay1 = 0.0_MK
    IF( j1 .GT. ncell(2)+nghost(4) ) j1  = ncell(2)+nghost(4)
    IF( j2 .GT. ncell(2)+nghost(4) ) ay2 = 0.0_MK
    IF( j2 .GT. ncell(2)+nghost(4) ) j2  = ncell(2)+nghost(4)
    IF( j3 .GT. ncell(2)+nghost(4) ) ay3 = 0.0_MK
    IF( j3 .GT. ncell(2)+nghost(4) ) j3  = ncell(2)+nghost(4)
    IF( j4 .GT. ncell(2)+nghost(4) ) ay4 = 0.0_MK
    IF( j4 .GT. ncell(2)+nghost(4) ) j4  = ncell(2)+nghost(4)
    IF( j5 .GT. ncell(2)+nghost(4) ) ay5 = 0.0_MK
    IF( j5 .GT. ncell(2)+nghost(4) ) j5  = ncell(2)+nghost(4)
    IF( j6 .GT. ncell(2)+nghost(4) ) ay6 = 0.0_MK
    IF( j6 .GT. ncell(2)+nghost(4) ) j6  = ncell(2)+nghost(4)

    IF( k1 .LT. 1-nghost(5) ) az1 = 0.0_MK
    IF( k1 .LT. 1-nghost(5) ) k1  = 1-nghost(5)
    IF( k2 .LT. 1-nghost(5) ) az2 = 0.0_MK
    IF( k2 .LT. 1-nghost(5) ) k2  = 1-nghost(5)
    IF( k3 .LT. 1-nghost(5) ) az3 = 0.0_MK
    IF( k3 .LT. 1-nghost(5) ) k3  = 1-nghost(5)
    IF( k4 .LT. 1-nghost(5) ) az4 = 0.0_MK
    IF( k4 .LT. 1-nghost(5) ) k4  = 1-nghost(5)
    IF( k5 .LT. 1-nghost(5) ) ay5 = 0.0_MK
    IF( k5 .LT. 1-nghost(5) ) j5  = 1-nghost(5)
    IF( k6 .LT. 1-nghost(5) ) ay6 = 0.0_MK
    IF( k6 .LT. 1-nghost(5) ) j6  = 1-nghost(5)

    IF( k1 .GT. ncell(3)+nghost(6) ) az1 = 0.0_MK
    IF( k1 .GT. ncell(3)+nghost(6) ) k1  = ncell(3)+nghost(6)
    IF( k2 .GT. ncell(3)+nghost(6) ) az2 = 0.0_MK
    IF( k2 .GT. ncell(3)+nghost(6) ) k2  = ncell(3)+nghost(6)
    IF( k3 .GT. ncell(3)+nghost(6) ) az3 = 0.0_MK
    IF( k3 .GT. ncell(3)+nghost(6) ) k3  = ncell(3)+nghost(6)
    IF( k4 .GT. ncell(3)+nghost(6) ) az4 = 0.0_MK
    IF( k4 .GT. ncell(3)+nghost(6) ) k4  = ncell(3)+nghost(6)
    IF( k5 .GT. ncell(3)+nghost(6) ) az5 = 0.0_MK
    IF( k5 .GT. ncell(3)+nghost(6) ) k5  = ncell(3)+nghost(6)
    IF( k6 .GT. ncell(3)+nghost(6) ) az6 = 0.0_MK
    IF( k6 .GT. ncell(3)+nghost(6) ) k6  = ncell(3)+nghost(6)

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 3D kernel and apply to the field
!---------------------------------------------------------------------------------!
    part(1:nvar,i) = ax1*ay1*az1*mesh(1:nvar,i1,j1,k1) + &
                   & ax2*ay1*az1*mesh(1:nvar,i2,j1,k1) + &
                   & ax3*ay1*az1*mesh(1:nvar,i3,j1,k1) + &
                   & ax4*ay1*az1*mesh(1:nvar,i4,j1,k1) + &
                   & ax5*ay1*az1*mesh(1:nvar,i5,j1,k1) + &
                   & ax6*ay1*az1*mesh(1:nvar,i6,j1,k1) + &
                   & ax1*ay2*az1*mesh(1:nvar,i1,j2,k1) + &
                   & ax2*ay2*az1*mesh(1:nvar,i2,j2,k1) + &
                   & ax3*ay2*az1*mesh(1:nvar,i3,j2,k1) + &
                   & ax4*ay2*az1*mesh(1:nvar,i4,j2,k1) + &
                   & ax5*ay2*az1*mesh(1:nvar,i5,j2,k1) + &
                   & ax6*ay2*az1*mesh(1:nvar,i6,j2,k1) + &
                   & ax1*ay3*az1*mesh(1:nvar,i1,j3,k1) + &
                   & ax2*ay3*az1*mesh(1:nvar,i2,j3,k1) + &
                   & ax3*ay3*az1*mesh(1:nvar,i3,j3,k1) + &
                   & ax4*ay3*az1*mesh(1:nvar,i4,j3,k1) + &
                   & ax5*ay3*az1*mesh(1:nvar,i5,j3,k1) + &
                   & ax6*ay3*az1*mesh(1:nvar,i6,j3,k1) + &
                   & ax1*ay4*az1*mesh(1:nvar,i1,j4,k1) + &
                   & ax2*ay4*az1*mesh(1:nvar,i2,j4,k1) + &
                   & ax3*ay4*az1*mesh(1:nvar,i3,j4,k1) + &
                   & ax4*ay4*az1*mesh(1:nvar,i4,j4,k1) + &
                   & ax5*ay4*az1*mesh(1:nvar,i5,j4,k1) + &
                   & ax6*ay4*az1*mesh(1:nvar,i6,j4,k1) + &
                   & ax1*ay5*az1*mesh(1:nvar,i1,j5,k1) + &
                   & ax2*ay5*az1*mesh(1:nvar,i2,j5,k1) + &
                   & ax3*ay5*az1*mesh(1:nvar,i3,j5,k1) + &
                   & ax4*ay5*az1*mesh(1:nvar,i4,j5,k1) + &
                   & ax5*ay5*az1*mesh(1:nvar,i5,j5,k1) + &
                   & ax6*ay5*az1*mesh(1:nvar,i6,j5,k1) + &
                   & ax1*ay6*az1*mesh(1:nvar,i1,j6,k1) + &
                   & ax2*ay6*az1*mesh(1:nvar,i2,j6,k1) + &
                   & ax3*ay6*az1*mesh(1:nvar,i3,j6,k1) + &
                   & ax4*ay6*az1*mesh(1:nvar,i4,j6,k1) + &
                   & ax5*ay6*az1*mesh(1:nvar,i5,j6,k1) + &
                   & ax6*ay6*az1*mesh(1:nvar,i6,j6,k1) + &
                   & ax1*ay1*az2*mesh(1:nvar,i1,j1,k2) + &
                   & ax2*ay1*az2*mesh(1:nvar,i2,j1,k2) + &
                   & ax3*ay1*az2*mesh(1:nvar,i3,j1,k2) + &
                   & ax4*ay1*az2*mesh(1:nvar,i4,j1,k2) + &
                   & ax5*ay1*az2*mesh(1:nvar,i5,j1,k2) + &
                   & ax6*ay1*az2*mesh(1:nvar,i6,j1,k2) + &
                   & ax1*ay2*az2*mesh(1:nvar,i1,j2,k2) + &
                   & ax2*ay2*az2*mesh(1:nvar,i2,j2,k2) + &
                   & ax3*ay2*az2*mesh(1:nvar,i3,j2,k2) + &
                   & ax4*ay2*az2*mesh(1:nvar,i4,j2,k2) + &
                   & ax5*ay2*az2*mesh(1:nvar,i5,j2,k2) + &
                   & ax6*ay2*az2*mesh(1:nvar,i6,j2,k2) + &
                   & ax1*ay3*az2*mesh(1:nvar,i1,j3,k2) + &
                   & ax2*ay3*az2*mesh(1:nvar,i2,j3,k2) + &
                   & ax3*ay3*az2*mesh(1:nvar,i3,j3,k2) + &
                   & ax4*ay3*az2*mesh(1:nvar,i4,j3,k2) + &
                   & ax5*ay3*az2*mesh(1:nvar,i5,j3,k2) + &
                   & ax6*ay3*az2*mesh(1:nvar,i6,j3,k2) + &
                   & ax1*ay4*az2*mesh(1:nvar,i1,j4,k2) + &
                   & ax2*ay4*az2*mesh(1:nvar,i2,j4,k2) + &
                   & ax3*ay4*az2*mesh(1:nvar,i3,j4,k2) + &
                   & ax4*ay4*az2*mesh(1:nvar,i4,j4,k2) + &
                   & ax5*ay4*az2*mesh(1:nvar,i5,j4,k2) + &
                   & ax6*ay4*az2*mesh(1:nvar,i6,j4,k2) + &
                   & ax1*ay5*az2*mesh(1:nvar,i1,j5,k2) + &
                   & ax2*ay5*az2*mesh(1:nvar,i2,j5,k2) + &
                   & ax3*ay5*az2*mesh(1:nvar,i3,j5,k2) + &
                   & ax4*ay5*az2*mesh(1:nvar,i4,j5,k2) + &
                   & ax5*ay5*az2*mesh(1:nvar,i5,j5,k2) + &
                   & ax6*ay5*az2*mesh(1:nvar,i6,j5,k2) + &
                   & ax1*ay6*az2*mesh(1:nvar,i1,j6,k2) + &
                   & ax2*ay6*az2*mesh(1:nvar,i2,j6,k2) + &
                   & ax3*ay6*az2*mesh(1:nvar,i3,j6,k2) + &
                   & ax4*ay6*az2*mesh(1:nvar,i4,j6,k2) + &
                   & ax5*ay6*az2*mesh(1:nvar,i5,j6,k2) + &
                   & ax6*ay6*az2*mesh(1:nvar,i6,j6,k2) + &
                   & ax1*ay1*az3*mesh(1:nvar,i1,j1,k3) + &
                   & ax2*ay1*az3*mesh(1:nvar,i2,j1,k3) + &
                   & ax3*ay1*az3*mesh(1:nvar,i3,j1,k3) + &
                   & ax4*ay1*az3*mesh(1:nvar,i4,j1,k3) + &
                   & ax5*ay1*az3*mesh(1:nvar,i5,j1,k3) + &
                   & ax6*ay1*az3*mesh(1:nvar,i6,j1,k3) + &
                   & ax1*ay2*az3*mesh(1:nvar,i1,j2,k3) + &
                   & ax2*ay2*az3*mesh(1:nvar,i2,j2,k3) + &
                   & ax3*ay2*az3*mesh(1:nvar,i3,j2,k3) + &
                   & ax4*ay2*az3*mesh(1:nvar,i4,j2,k3) + &
                   & ax5*ay2*az3*mesh(1:nvar,i5,j2,k3) + &
                   & ax6*ay2*az3*mesh(1:nvar,i6,j2,k3) + &
                   & ax1*ay3*az3*mesh(1:nvar,i1,j3,k3) + &
                   & ax2*ay3*az3*mesh(1:nvar,i2,j3,k3) + &
                   & ax3*ay3*az3*mesh(1:nvar,i3,j3,k3) + &
                   & ax4*ay3*az3*mesh(1:nvar,i4,j3,k3) + &
                   & ax5*ay3*az3*mesh(1:nvar,i5,j3,k3) + &
                   & ax6*ay3*az3*mesh(1:nvar,i6,j3,k3) + &
                   & ax1*ay4*az3*mesh(1:nvar,i1,j4,k3) + &
                   & ax2*ay4*az3*mesh(1:nvar,i2,j4,k3) + &
                   & ax3*ay4*az3*mesh(1:nvar,i3,j4,k3) + &
                   & ax4*ay4*az3*mesh(1:nvar,i4,j4,k3) + &
                   & ax5*ay4*az3*mesh(1:nvar,i5,j4,k3) + &
                   & ax6*ay4*az3*mesh(1:nvar,i6,j4,k3) + &
                   & ax1*ay5*az3*mesh(1:nvar,i1,j5,k3) + &
                   & ax2*ay5*az3*mesh(1:nvar,i2,j5,k3) + &
                   & ax3*ay5*az3*mesh(1:nvar,i3,j5,k3) + &
                   & ax4*ay5*az3*mesh(1:nvar,i4,j5,k3) + &
                   & ax5*ay5*az3*mesh(1:nvar,i5,j5,k3) + &
                   & ax6*ay5*az3*mesh(1:nvar,i6,j5,k3) + &
                   & ax1*ay6*az3*mesh(1:nvar,i1,j6,k3) + &
                   & ax2*ay6*az3*mesh(1:nvar,i2,j6,k3) + &
                   & ax3*ay6*az3*mesh(1:nvar,i3,j6,k3) + &
                   & ax4*ay6*az3*mesh(1:nvar,i4,j6,k3) + &
                   & ax5*ay6*az3*mesh(1:nvar,i5,j6,k3) + &
                   & ax6*ay6*az3*mesh(1:nvar,i6,j6,k3) + &
                   & ax1*ay1*az4*mesh(1:nvar,i1,j1,k4) + &
                   & ax2*ay1*az4*mesh(1:nvar,i2,j1,k4) + &
                   & ax3*ay1*az4*mesh(1:nvar,i3,j1,k4) + &
                   & ax4*ay1*az4*mesh(1:nvar,i4,j1,k4) + &
                   & ax5*ay1*az4*mesh(1:nvar,i5,j1,k4) + &
                   & ax6*ay1*az4*mesh(1:nvar,i6,j1,k4) + &
                   & ax1*ay2*az4*mesh(1:nvar,i1,j2,k4) + &
                   & ax2*ay2*az4*mesh(1:nvar,i2,j2,k4) + &
                   & ax3*ay2*az4*mesh(1:nvar,i3,j2,k4) + &
                   & ax4*ay2*az4*mesh(1:nvar,i4,j2,k4) + &
                   & ax5*ay2*az4*mesh(1:nvar,i5,j2,k4) + &
                   & ax6*ay2*az4*mesh(1:nvar,i6,j2,k4) + &
                   & ax1*ay3*az4*mesh(1:nvar,i1,j3,k4) + &
                   & ax2*ay3*az4*mesh(1:nvar,i2,j3,k4) + &
                   & ax3*ay3*az4*mesh(1:nvar,i3,j3,k4) + &
                   & ax4*ay3*az4*mesh(1:nvar,i4,j3,k4) + &
                   & ax5*ay3*az4*mesh(1:nvar,i5,j3,k4) + &
                   & ax6*ay3*az4*mesh(1:nvar,i6,j3,k4) + &
                   & ax1*ay4*az4*mesh(1:nvar,i1,j4,k4) + &
                   & ax2*ay4*az4*mesh(1:nvar,i2,j4,k4) + &
                   & ax3*ay4*az4*mesh(1:nvar,i3,j4,k4) + &
                   & ax4*ay4*az4*mesh(1:nvar,i4,j4,k4) + &
                   & ax5*ay4*az4*mesh(1:nvar,i5,j4,k4) + &
                   & ax6*ay4*az4*mesh(1:nvar,i6,j4,k4) + &
                   & ax1*ay5*az4*mesh(1:nvar,i1,j5,k4) + &
                   & ax2*ay5*az4*mesh(1:nvar,i2,j5,k4) + &
                   & ax3*ay5*az4*mesh(1:nvar,i3,j5,k4) + &
                   & ax4*ay5*az4*mesh(1:nvar,i4,j5,k4) + &
                   & ax5*ay5*az4*mesh(1:nvar,i5,j5,k4) + &
                   & ax6*ay5*az4*mesh(1:nvar,i6,j5,k4) + &
                   & ax1*ay6*az4*mesh(1:nvar,i1,j6,k4) + &
                   & ax2*ay6*az4*mesh(1:nvar,i2,j6,k4) + &
                   & ax3*ay6*az4*mesh(1:nvar,i3,j6,k4) + &
                   & ax4*ay6*az4*mesh(1:nvar,i4,j6,k4) + &
                   & ax5*ay6*az4*mesh(1:nvar,i5,j6,k4) + &
                   & ax6*ay6*az4*mesh(1:nvar,i6,j6,k4) + &
                   & ax1*ay1*az5*mesh(1:nvar,i1,j1,k5) + &
                   & ax2*ay1*az5*mesh(1:nvar,i2,j1,k5) + &
                   & ax3*ay1*az5*mesh(1:nvar,i3,j1,k5) + &
                   & ax4*ay1*az5*mesh(1:nvar,i4,j1,k5) + &
                   & ax5*ay1*az5*mesh(1:nvar,i5,j1,k5) + &
                   & ax6*ay1*az5*mesh(1:nvar,i6,j1,k5) + &
                   & ax1*ay2*az5*mesh(1:nvar,i1,j2,k5) + &
                   & ax2*ay2*az5*mesh(1:nvar,i2,j2,k5) + &
                   & ax3*ay2*az5*mesh(1:nvar,i3,j2,k5) + &
                   & ax4*ay2*az5*mesh(1:nvar,i4,j2,k5) + &
                   & ax5*ay2*az5*mesh(1:nvar,i5,j2,k5) + &
                   & ax6*ay2*az5*mesh(1:nvar,i6,j2,k5) + &
                   & ax1*ay3*az5*mesh(1:nvar,i1,j3,k5) + &
                   & ax2*ay3*az5*mesh(1:nvar,i2,j3,k5) + &
                   & ax3*ay3*az5*mesh(1:nvar,i3,j3,k5) + &
                   & ax4*ay3*az5*mesh(1:nvar,i4,j3,k5) + &
                   & ax5*ay3*az5*mesh(1:nvar,i5,j3,k5) + &
                   & ax6*ay3*az5*mesh(1:nvar,i6,j3,k5) + &
                   & ax1*ay4*az5*mesh(1:nvar,i1,j4,k5) + &
                   & ax2*ay4*az5*mesh(1:nvar,i2,j4,k5) + &
                   & ax3*ay4*az5*mesh(1:nvar,i3,j4,k5) + &
                   & ax4*ay4*az5*mesh(1:nvar,i4,j4,k5) + &
                   & ax5*ay4*az5*mesh(1:nvar,i5,j4,k5) + &
                   & ax6*ay4*az5*mesh(1:nvar,i6,j4,k5) + &
                   & ax1*ay5*az5*mesh(1:nvar,i1,j5,k5) + &
                   & ax2*ay5*az5*mesh(1:nvar,i2,j5,k5) + &
                   & ax3*ay5*az5*mesh(1:nvar,i3,j5,k5) + &
                   & ax4*ay5*az5*mesh(1:nvar,i4,j5,k5) + &
                   & ax5*ay5*az5*mesh(1:nvar,i5,j5,k5) + &
                   & ax6*ay5*az5*mesh(1:nvar,i6,j5,k5) + &
                   & ax1*ay6*az5*mesh(1:nvar,i1,j6,k5) + &
                   & ax2*ay6*az5*mesh(1:nvar,i2,j6,k5) + &
                   & ax3*ay6*az5*mesh(1:nvar,i3,j6,k5) + &
                   & ax4*ay6*az5*mesh(1:nvar,i4,j6,k5) + &
                   & ax5*ay6*az5*mesh(1:nvar,i5,j6,k5) + &
                   & ax6*ay6*az5*mesh(1:nvar,i6,j6,k5) + &
                   & ax1*ay1*az6*mesh(1:nvar,i1,j1,k6) + &
                   & ax2*ay1*az6*mesh(1:nvar,i2,j1,k6) + &
                   & ax3*ay1*az6*mesh(1:nvar,i3,j1,k6) + &
                   & ax4*ay1*az6*mesh(1:nvar,i4,j1,k6) + &
                   & ax5*ay1*az6*mesh(1:nvar,i5,j1,k6) + &
                   & ax6*ay1*az6*mesh(1:nvar,i6,j1,k6) + &
                   & ax1*ay2*az6*mesh(1:nvar,i1,j2,k6) + &
                   & ax2*ay2*az6*mesh(1:nvar,i2,j2,k6) + &
                   & ax3*ay2*az6*mesh(1:nvar,i3,j2,k6) + &
                   & ax4*ay2*az6*mesh(1:nvar,i4,j2,k6) + &
                   & ax5*ay2*az6*mesh(1:nvar,i5,j2,k6) + &
                   & ax6*ay2*az6*mesh(1:nvar,i6,j2,k6) + &
                   & ax1*ay3*az6*mesh(1:nvar,i1,j3,k6) + &
                   & ax2*ay3*az6*mesh(1:nvar,i2,j3,k6) + &
                   & ax3*ay3*az6*mesh(1:nvar,i3,j3,k6) + &
                   & ax4*ay3*az6*mesh(1:nvar,i4,j3,k6) + &
                   & ax5*ay3*az6*mesh(1:nvar,i5,j3,k6) + &
                   & ax6*ay3*az6*mesh(1:nvar,i6,j3,k6) + &
                   & ax1*ay4*az6*mesh(1:nvar,i1,j4,k6) + &
                   & ax2*ay4*az6*mesh(1:nvar,i2,j4,k6) + &
                   & ax3*ay4*az6*mesh(1:nvar,i3,j4,k6) + &
                   & ax4*ay4*az6*mesh(1:nvar,i4,j4,k6) + &
                   & ax5*ay4*az6*mesh(1:nvar,i5,j4,k6) + &
                   & ax6*ay4*az6*mesh(1:nvar,i6,j4,k6) + &
                   & ax1*ay5*az6*mesh(1:nvar,i1,j5,k6) + &
                   & ax2*ay5*az6*mesh(1:nvar,i2,j5,k6) + &
                   & ax3*ay5*az6*mesh(1:nvar,i3,j5,k6) + &
                   & ax4*ay5*az6*mesh(1:nvar,i4,j5,k6) + &
                   & ax5*ay5*az6*mesh(1:nvar,i5,j5,k6) + &
                   & ax6*ay5*az6*mesh(1:nvar,i6,j5,k6) + &
                   & ax1*ay6*az6*mesh(1:nvar,i1,j6,k6) + &
                   & ax2*ay6*az6*mesh(1:nvar,i2,j6,k6) + &
                   & ax3*ay6*az6*mesh(1:nvar,i3,j6,k6) + &
                   & ax4*ay6*az6*mesh(1:nvar,i4,j6,k6) + &
                   & ax5*ay6*az6*mesh(1:nvar,i5,j6,k6) + &
                   & ax6*ay6*az6*mesh(1:nvar,i6,j6,k6)

  END DO ! particle

  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'Interpolation scheme/order unknown.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_interp_mesh_particle_var
