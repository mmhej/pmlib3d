!---------------------------------------------------------------------------------!
! pmlib_interp_mesh_mesh.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_interp_mesh_mesh( topo_from,topo_to,mesh_from,mesh_to, &
                                   & ierr,clear )

USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_communication

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN) :: topo_from
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN) :: topo_to
  REAL(MK), DIMENSION(:,:,:,:), POINTER, INTENT(IN)     :: mesh_from
  REAL(MK), DIMENSION(:,:,:,:), POINTER, INTENT(OUT)    :: mesh_to
  INTEGER, INTENT(OUT)                                  :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                         :: clear

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=16)                      :: caller = 'pmlib_interp_mesh_mesh'
  INTEGER                                :: i,j,k, nvar
  REAL(MK),DIMENSION(pmlib_ndim)         :: xmin_to, xmin_from
  REAL(MK),DIMENSION(pmlib_ndim)         :: xmax_to, xmax_from
  REAL(MK),DIMENSION(pmlib_ndim)         :: dx_to, dx_from
  INTEGER,DIMENSION(pmlib_ndim)          :: ncell_to, ncell_from
  INTEGER,DIMENSION(2*pmlib_ndim)        :: nghost_to, nghost_from
  INTEGER,DIMENSION(2*pmlib_ndim)        :: bc_to, bc_from
  REAL(MK)                               :: px, py, pz
  REAL(MK)                               :: reldx

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
  nvar = SIZE(mesh_from,1)
  IF(nvar .NE. SIZE(mesh_to,1))THEN
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
         & 'Dimension of interpolation variables does not correspond.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Get topologies
!---------------------------------------------------------------------------------!
  dx_to     = topo_to(mpi_rank)%dx
  xmin_to   = topo_to(mpi_rank)%xmin
  xmax_to   = topo_to(mpi_rank)%xmax
  ncell_to  = topo_to(mpi_rank)%ncell
  nghost_to = topo_to(mpi_rank)%nghost
  bc_to     = topo_to(mpi_rank)%bound_cond

  dx_from     = topo_from(mpi_rank)%dx
  xmin_from   = topo_from(mpi_rank)%xmin
  xmax_from   = topo_from(mpi_rank)%xmax
  ncell_from  = topo_from(mpi_rank)%ncell
  nghost_from = topo_from(mpi_rank)%nghost
  bc_from     = topo_from(mpi_rank)%bound_cond

  reldx = dx_from(1)/dx_to(1)

!---------------------------------------------------------------------------------!
! Clear mesh if specified
!---------------------------------------------------------------------------------!
  IF(PRESENT(clear) .AND. clear)THEN
    mesh_to = 0.0_MK
  END IF

!---------------------------------------------------------------------------------!
! Patch to parent
!---------------------------------------------------------------------------------!
  IF(reldx .LE. 0.9999_MK)THEN

    IF(pmlib_interpolation_order .EQ. 3)THEN
!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
      c_1_dx = 1.0_MK/dx_to(1)
      c_1_dy = 1.0_MK/dx_to(2)
      c_1_dz = 1.0_MK/dx_to(3)

!---------------------------------------------------------------------------------!
! Interpolate
!---------------------------------------------------------------------------------!
      DO k = 1-nghost_from(5),ncell_from(3)+nghost_from(6)
        pz = xmin_from(3) + (REAL(k - 1,MK) + 0.5_MK) *dx_from(3)

        DO j = 1-nghost_from(3),ncell_from(2)+nghost_from(4)
          py = xmin_from(2) + (REAL(j - 1,MK) + 0.5_MK) *dx_from(2)

          DO i = 1-nghost_from(1),ncell_from(1)+nghost_from(2)
            px = xmin_from(1) + (REAL(i - 1,MK) + 0.5_MK) *dx_from(1)

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
            i2 = NINT( (px - xmin_to(1)) * c_1_dx )
            j2 = NINT( (py - xmin_to(2)) * c_1_dy )
            k2 = NINT( (pz - xmin_to(3)) * c_1_dz )

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

            IF( i4 .GE. 1-nghost_to(1) .AND.                                      &
              & i1 .LE. ncell_to(1)+nghost_to(2) .AND.                            &
              & j4 .GE. 1-nghost_to(3) .AND.                                      &
              & j1 .LE. ncell_to(2)+nghost_to(4) .AND.                            &
              & k4 .GE. 1-nghost_to(5) .AND.                                      &
              & k1 .LE. ncell_to(3)+nghost_to(6) )THEN

!---------------------------------------------------------------------------------!
! Find normalised distance to other cells
!---------------------------------------------------------------------------------!
              dx2 = ( px - xmin_to(1) - dx_to(1)*(REAL(i2,MK) - 0.5_MK) ) * c_1_dx
              dy2 = ( py - xmin_to(2) - dx_to(2)*(REAL(j2,MK) - 0.5_MK) ) * c_1_dy
              dz2 = ( pz - xmin_to(3) - dx_to(3)*(REAL(k2,MK) - 0.5_MK) ) * c_1_dz

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
              ax1 = reldx * ( 0.5_MK * (2.0_MK - dx1)**2 * (1.0_MK - dx1) )
              ax2 = reldx * ( 1.0_MK - 2.5_MK*dx2**2 + 1.5_MK*dx2**3 )
              ax3 = reldx * ( 1.0_MK - 2.5_MK*dx3**2 + 1.5_MK*dx3**3 )
              ax4 = reldx * ( 0.5_MK * (2.0_MK - dx4)**2 * (1.0_MK - dx4) )

              ay1 = reldx * ( 0.5_MK * (2.0_MK - dy1)**2 * (1.0_MK - dy1) )
              ay2 = reldx * ( 1.0_MK - 2.5_MK*dy2**2 + 1.5_MK*dy2**3 )
              ay3 = reldx * ( 1.0_MK - 2.5_MK*dy3**2 + 1.5_MK*dy3**3 )
              ay4 = reldx * ( 0.5_MK * (2.0_MK - dy4)**2 * (1.0_MK - dy4) )

              az1 = reldx * ( 0.5_MK * (2.0_MK - dz1)**2 * (1.0_MK - dz1) )
              az2 = reldx * ( 1.0_MK - 2.5_MK*dz2**2 + 1.5_MK*dz2**3 )
              az3 = reldx * ( 1.0_MK - 2.5_MK*dz3**2 + 1.5_MK*dz3**3 )
              az4 = reldx * ( 0.5_MK * (2.0_MK - dz4)**2 * (1.0_MK - dz4) )

!---------------------------------------------------------------------------------!
! Cancel stencils that are not supported by the mesh
!---------------------------------------------------------------------------------!
              IF(bc_to(1) .EQ. 0)THEN
                IF( i1 .LT. 1-nghost_to(1) ) ax1 = 0.0_MK
                IF( i1 .LT. 1-nghost_to(1) ) i1  = 1-nghost_to(1)
                IF( i2 .LT. 1-nghost_to(1) ) ax2 = 0.0_MK
                IF( i2 .LT. 1-nghost_to(1) ) i2  = 1-nghost_to(1)
                IF( i3 .LT. 1-nghost_to(1) ) ax3 = 0.0_MK
                IF( i3 .LT. 1-nghost_to(1) ) i3  = 1-nghost_to(1)
                IF( i4 .LT. 1-nghost_to(1) ) ax4 = 0.0_MK
                IF( i4 .LT. 1-nghost_to(1) ) i4  = 1-nghost_to(1)
              ELSE
                IF( px .LT. xmin_to(1) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  i1  = 1-nghost_to(1)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  i2  = 1-nghost_to(1)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  i3  = 1-nghost_to(1)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  i4  = 1-nghost_to(1)
                END IF
              END IF

              IF(bc_to(2) .EQ. 0)THEN
                IF(i1 .GT. ncell_to(1)+nghost_to(2)) ax1 = 0.0_MK
                IF(i1 .GT. ncell_to(1)+nghost_to(2)) i1  = ncell_to(1)+nghost_to(2)
                IF(i2 .GT. ncell_to(1)+nghost_to(2)) ax2 = 0.0_MK
                IF(i2 .GT. ncell_to(1)+nghost_to(2)) i2  = ncell_to(1)+nghost_to(2)
                IF(i3 .GT. ncell_to(1)+nghost_to(2)) ax3 = 0.0_MK
                IF(i3 .GT. ncell_to(1)+nghost_to(2)) i3  = ncell_to(1)+nghost_to(2)
                IF(i4 .GT. ncell_to(1)+nghost_to(2)) ax4 = 0.0_MK
                IF(i4 .GT. ncell_to(1)+nghost_to(2)) i4  = ncell_to(1)+nghost_to(2)
              ELSE
                IF( px .GT. xmax_to(1) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  i1  = ncell_to(1)+nghost_to(2)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  i2  = ncell_to(1)+nghost_to(2)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  i3  = ncell_to(1)+nghost_to(2)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  i4  = ncell_to(1)+nghost_to(2)
                END IF
              END IF

              IF(bc_to(3) .EQ. 0)THEN
                IF( j1 .LT. 1-nghost_to(3) ) ay1 = 0.0_MK
                IF( j1 .LT. 1-nghost_to(3) ) j1  = 1-nghost_to(3)
                IF( j2 .LT. 1-nghost_to(3) ) ay2 = 0.0_MK
                IF( j2 .LT. 1-nghost_to(3) ) j2  = 1-nghost_to(3)
                IF( j3 .LT. 1-nghost_to(3) ) ay3 = 0.0_MK
                IF( j3 .LT. 1-nghost_to(3) ) j3  = 1-nghost_to(3)
                IF( j4 .LT. 1-nghost_to(3) ) ay4 = 0.0_MK
                IF( j4 .LT. 1-nghost_to(3) ) j4  = 1-nghost_to(3)
              ELSE
                IF( py .LT. xmin_to(2) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  j1  = 1-nghost_to(3)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  j2  = 1-nghost_to(3)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  j3  = 1-nghost_to(3)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  j4  = 1-nghost_to(3)
                END IF
              END IF

              IF(bc_to(4) .EQ. 0)THEN
                IF(j1 .GT. ncell_to(2)+nghost_to(4)) ay1 = 0.0_MK
                IF(j1 .GT. ncell_to(2)+nghost_to(4)) j1  = ncell_to(2)+nghost_to(4)
                IF(j2 .GT. ncell_to(2)+nghost_to(4)) ay2 = 0.0_MK
                IF(j2 .GT. ncell_to(2)+nghost_to(4)) j2  = ncell_to(2)+nghost_to(4)
                IF(j3 .GT. ncell_to(2)+nghost_to(4)) ay3 = 0.0_MK
                IF(j3 .GT. ncell_to(2)+nghost_to(4)) j3  = ncell_to(2)+nghost_to(4)
                IF(j4 .GT. ncell_to(2)+nghost_to(4)) ay4 = 0.0_MK
                IF(j4 .GT. ncell_to(2)+nghost_to(4)) j4  = ncell_to(2)+nghost_to(4)
              ELSE
                IF( py .GT. xmax_to(2) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  j1  = ncell_to(2)+nghost_to(4)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  j2  = ncell_to(2)+nghost_to(4)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  j3  = ncell_to(2)+nghost_to(4)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  j4  = ncell_to(2)+nghost_to(4)
                END IF
              END IF

              IF(bc_to(5) .EQ. 0)THEN
                IF( k1 .LT. 1-nghost_to(5) ) az1 = 0.0_MK
                IF( k1 .LT. 1-nghost_to(5) ) k1  = 1-nghost_to(5)
                IF( k2 .LT. 1-nghost_to(5) ) az2 = 0.0_MK
                IF( k2 .LT. 1-nghost_to(5) ) k2  = 1-nghost_to(5)
                IF( k3 .LT. 1-nghost_to(5) ) az3 = 0.0_MK
                IF( k3 .LT. 1-nghost_to(5) ) k3  = 1-nghost_to(5)
                IF( k4 .LT. 1-nghost_to(5) ) az4 = 0.0_MK
                IF( k4 .LT. 1-nghost_to(5) ) k4  = 1-nghost_to(5)
              ELSE
                IF( pz .LT. xmin_to(3) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  k1  = 1-nghost_to(5)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  k2  = 1-nghost_to(5)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  k3  = 1-nghost_to(5)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  k4  = 1-nghost_to(5)
                END IF
              END IF

              IF(bc_to(6) .EQ. 0)THEN
                IF(k1 .GT. ncell_to(3)+nghost_to(6)) az1 = 0.0_MK
                IF(k1 .GT. ncell_to(3)+nghost_to(6)) k1  = ncell_to(3)+nghost_to(6)
                IF(k2 .GT. ncell_to(3)+nghost_to(6)) az2 = 0.0_MK
                IF(k2 .GT. ncell_to(3)+nghost_to(6)) k2  = ncell_to(3)+nghost_to(6)
                IF(k3 .GT. ncell_to(3)+nghost_to(6)) az3 = 0.0_MK
                IF(k3 .GT. ncell_to(3)+nghost_to(6)) k3  = ncell_to(3)+nghost_to(6)
                IF(k4 .GT. ncell_to(3)+nghost_to(6)) az4 = 0.0_MK
                IF(k4 .GT. ncell_to(3)+nghost_to(6)) k4  = ncell_to(3)+nghost_to(6)
              ELSE
                IF( pz .GT. xmax_to(3) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  k1  = ncell_to(3)+nghost_to(6)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  k2  = ncell_to(3)+nghost_to(6)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  k3  = ncell_to(3)+nghost_to(6)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  k4  = ncell_to(3)+nghost_to(6)
                END IF
              END IF

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 3D kernel and apply to the field
!---------------------------------------------------------------------------------!
              mesh_to(1:nvar,i1,j1,k1) = mesh_to(1:nvar,i1,j1,k1) & 
                                       & + ax1*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k1) = mesh_to(1:nvar,i2,j1,k1)  & 
                                       & + ax2*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k1) = mesh_to(1:nvar,i3,j1,k1)  & 
                                       & + ax3*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k1) = mesh_to(1:nvar,i4,j1,k1)  & 
                                       & + ax4*ay1*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j2,k1) = mesh_to(1:nvar,i1,j2,k1)  & 
                                       & + ax1*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k1) = mesh_to(1:nvar,i2,j2,k1)  & 
                                       & + ax2*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k1) = mesh_to(1:nvar,i3,j2,k1)  & 
                                       & + ax3*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k1) = mesh_to(1:nvar,i4,j2,k1)  & 
                                       & + ax4*ay2*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k1) = mesh_to(1:nvar,i1,j3,k1)  & 
                                       & + ax1*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k1) = mesh_to(1:nvar,i2,j3,k1)  & 
                                       & + ax2*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k1) = mesh_to(1:nvar,i3,j3,k1)  & 
                                       & + ax3*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k1) = mesh_to(1:nvar,i4,j3,k1)  & 
                                       & + ax4*ay3*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k1) = mesh_to(1:nvar,i1,j4,k1)  & 
                                       & + ax1*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k1) = mesh_to(1:nvar,i2,j4,k1)  & 
                                       & + ax2*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k1) = mesh_to(1:nvar,i3,j4,k1)  & 
                                       & + ax3*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k1) = mesh_to(1:nvar,i4,j4,k1)  & 
                                       & + ax4*ay4*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k2) = mesh_to(1:nvar,i1,j1,k2)  & 
                                       & + ax1*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k2) = mesh_to(1:nvar,i2,j1,k2)  & 
                                       & + ax2*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k2) = mesh_to(1:nvar,i3,j1,k2)  & 
                                       & + ax3*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k2) = mesh_to(1:nvar,i4,j1,k2)  & 
                                       & + ax4*ay1*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j2,k2) = mesh_to(1:nvar,i1,j2,k2)  & 
                                       & + ax1*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k2) = mesh_to(1:nvar,i2,j2,k2)  & 
                                       & + ax2*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k2) = mesh_to(1:nvar,i3,j2,k2)  & 
                                       & + ax3*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k2) = mesh_to(1:nvar,i4,j2,k2)  & 
                                       & + ax4*ay2*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k2) = mesh_to(1:nvar,i1,j3,k2)  & 
                                       & + ax1*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k2) = mesh_to(1:nvar,i2,j3,k2)  & 
                                       & + ax2*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k2) = mesh_to(1:nvar,i3,j3,k2)  & 
                                       & + ax3*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k2) = mesh_to(1:nvar,i4,j3,k2)  & 
                                       & + ax4*ay3*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k2) = mesh_to(1:nvar,i1,j4,k2)  & 
                                       & + ax1*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k2) = mesh_to(1:nvar,i2,j4,k2)  & 
                                       & + ax2*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k2) = mesh_to(1:nvar,i3,j4,k2)  & 
                                       & + ax3*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k2) = mesh_to(1:nvar,i4,j4,k2)  & 
                                       & + ax4*ay4*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k3) = mesh_to(1:nvar,i1,j1,k3)  & 
                                       & + ax1*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k3) = mesh_to(1:nvar,i2,j1,k3)  & 
                                       & + ax2*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k3) = mesh_to(1:nvar,i3,j1,k3)  & 
                                       & + ax3*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k3) = mesh_to(1:nvar,i4,j1,k3)  & 
                                       & + ax4*ay1*az3*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j2,k3) = mesh_to(1:nvar,i1,j2,k3)  & 
                                       & + ax1*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k3) = mesh_to(1:nvar,i2,j2,k3)  & 
                                       & + ax2*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k3) = mesh_to(1:nvar,i3,j2,k3)  & 
                                       & + ax3*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k3) = mesh_to(1:nvar,i4,j2,k3)  & 
                                       & + ax4*ay2*az3*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k3) = mesh_to(1:nvar,i1,j3,k3)  & 
                                       & + ax1*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k3) = mesh_to(1:nvar,i2,j3,k3)  & 
                                       & + ax2*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k3) = mesh_to(1:nvar,i3,j3,k3)  & 
                                       & + ax3*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k3) = mesh_to(1:nvar,i4,j3,k3)  & 
                                       & + ax4*ay3*az3*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k3) = mesh_to(1:nvar,i1,j4,k3)  & 
                                       & + ax1*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k3) = mesh_to(1:nvar,i2,j4,k3)  & 
                                       & + ax2*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k3) = mesh_to(1:nvar,i3,j4,k3)  & 
                                       & + ax3*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k3) = mesh_to(1:nvar,i4,j4,k3)  & 
                                       & + ax4*ay4*az3*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k4) = mesh_to(1:nvar,i1,j1,k4)  & 
                                       & + ax1*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k4) = mesh_to(1:nvar,i2,j1,k4)  & 
                                       & + ax2*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k4) = mesh_to(1:nvar,i3,j1,k4)  & 
                                       & + ax3*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k4) = mesh_to(1:nvar,i4,j1,k4)  & 
                                       & + ax4*ay1*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j2,k4) = mesh_to(1:nvar,i1,j2,k4)  & 
                                       & + ax1*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k4) = mesh_to(1:nvar,i2,j2,k4)  & 
                                       & + ax2*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k4) = mesh_to(1:nvar,i3,j2,k4)  & 
                                       & + ax3*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k4) = mesh_to(1:nvar,i4,j2,k4)  & 
                                       & + ax4*ay2*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k4) = mesh_to(1:nvar,i1,j3,k4)  & 
                                       & + ax1*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k4) = mesh_to(1:nvar,i2,j3,k4)  & 
                                       & + ax2*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k4) = mesh_to(1:nvar,i3,j3,k4)  & 
                                       & + ax3*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k4) = mesh_to(1:nvar,i4,j3,k4)  & 
                                       & + ax4*ay3*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k4) = mesh_to(1:nvar,i1,j4,k4)  & 
                                       & + ax1*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k4) = mesh_to(1:nvar,i2,j4,k4)  & 
                                       & + ax2*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k4) = mesh_to(1:nvar,i3,j4,k4)  & 
                                       & + ax3*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k4) = mesh_to(1:nvar,i4,j4,k4)  & 
                                       & + ax4*ay4*az4*mesh_from(1:nvar,i,j,k)

            END IF
          END DO ! i
        END DO ! j
      END DO ! k

    ELSEIF(pmlib_interpolation_order .EQ. 4)THEN
!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
      c_1_dx = 1.0_MK/dx_to(1)
      c_1_dy = 1.0_MK/dx_to(2)
      c_1_dz = 1.0_MK/dx_to(3)

      c1 = -1.0_MK/24.0_MK
      c2 =  1.0_MK/24.0_MK 
      c3 = -1.0_MK/12.0_MK

!---------------------------------------------------------------------------------!
! Interpolate
!---------------------------------------------------------------------------------!
      DO k = 1-nghost_from(5),ncell_from(3)+nghost_from(6)
        pz = xmin_from(3) + (REAL(k - 1,MK) + 0.5_MK) *dx_from(3)

        DO j = 1-nghost_from(3),ncell_from(2)+nghost_from(4)
          py = xmin_from(2) + (REAL(j - 1,MK) + 0.5_MK) *dx_from(2)

          DO i = 1-nghost_from(1),ncell_from(1)+nghost_from(2)
            px = xmin_from(1) + (REAL(i - 1,MK) + 0.5_MK) *dx_from(1)

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
            i3 = NINT( (px - xmin_to(1)) * c_1_dx )
            j3 = NINT( (py - xmin_to(2)) * c_1_dy )
            k3 = NINT( (pz - xmin_to(3)) * c_1_dz )

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

            IF( i6 .GE. 1-nghost_to(1) .AND.                                      &
              & i1 .LE. ncell_to(1)+nghost_to(2) .AND.                            &
              & j6 .GE. 1-nghost_to(3) .AND.                                      &
              & j1 .LE. ncell_to(2)+nghost_to(4) .AND.                            &
              & k6 .GE. 1-nghost_to(5) .AND.                                      &
              & k1 .LE. ncell_to(3)+nghost_to(6) )THEN

!---------------------------------------------------------------------------------!
! Find normalised distance to other cells
!---------------------------------------------------------------------------------!
              dx3 = ( px - xmin_to(1) - dx_to(1)*(REAL(i3,MK) - 0.5_MK) ) * c_1_dx
              dy3 = ( py - xmin_to(2) - dx_to(2)*(REAL(j3,MK) - 0.5_MK) ) * c_1_dy
              dz3 = ( pz - xmin_to(3) - dx_to(3)*(REAL(k3,MK) - 0.5_MK) ) * c_1_dz

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
              ax1 = reldx * ( c1 * (dx1 - 2.0_MK)* (dx1 - 3.0_MK)**3 & 
                  & * (5.0_MK * dx1 - 8.0_MK) )
              ax2 = reldx * ( c2 * (dx2 - 1.0_MK) * (dx2 - 2.0_MK) &
                  & * (25.0_MK*dx2**3 - 114.0_MK*dx2**2 + 153.0_MK*dx2 - 48.0_MK) )
              ax3 = reldx * ( c3 * (dx3 - 1.0_MK) * (25.0_MK*dx3**4 & 
                  & - 38.0_MK*dx3**3 - 3.0_MK*dx3**2 + 12.0_MK*dx3 + 12.0_MK) )
              ax4 = reldx * ( c3 * (dx4 - 1.0_MK) * (25.0_MK*dx4**4 & 
                  & - 38.0_MK*dx4**3 - 3.0_MK*dx4**2 + 12.0_MK*dx4 + 12.0_MK) )
              ax5 = reldx * ( c2 * (dx5 - 1.0_MK) * (dx5 - 2.0_MK) &
                  & * (25.0_MK*dx5**3 - 114.0_MK*dx5**2 + 153.0_MK*dx5 - 48.0_MK) )
              ax6 = reldx * ( c1 * (dx6 - 2.0_MK)* (dx6 - 3.0_MK)**3 & 
                  & * (5.0_MK * dx6 - 8.0_MK) )

              ay1 = reldx * ( c1 * (dy1 - 2.0_MK)* (dy1 - 3.0_MK)**3 & 
                  & * (5.0_MK * dy1 - 8.0_MK) )
              ay2 = reldx * ( c2 * (dy2 - 1.0_MK) * (dy2 - 2.0_MK) &
                  & * (25.0_MK*dy2**3 - 114.0_MK*dy2**2 + 153.0_MK*dy2 - 48.0_MK) )
              ay3 = reldx * ( c3 * (dy3 - 1.0_MK) * (25.0_MK*dy3**4 & 
                  & - 38.0_MK*dy3**3 - 3.0_MK*dy3**2 + 12.0_MK*dy3 + 12.0_MK) )
              ay4 = reldx * ( c3 * (dy4 - 1.0_MK) * (25.0_MK*dy4**4 & 
                  & - 38.0_MK*dy4**3 - 3.0_MK*dy4**2 + 12.0_MK*dy4 + 12.0_MK) )
              ay5 = reldx * ( c2 * (dy5 - 1.0_MK) * (dy5 - 2.0_MK) &
                  & * (25.0_MK*dy5**3 - 114.0_MK*dy5**2 + 153.0_MK*dy5 - 48.0_MK) )
              ay6 = reldx * ( c1 * (dy6 - 2.0_MK)* (dy6 - 3.0_MK)**3 & 
                  & * (5.0_MK * dy6 - 8.0_MK) )

              az1 = reldx * ( c1 * (dz1 - 2.0_MK)* (dz1 - 3.0_MK)**3 & 
                  & * (5.0_MK * dz1 - 8.0_MK) )
              az2 = reldx * ( c2 * (dz2 - 1.0_MK) * (dz2 - 2.0_MK) &
                  & * (25.0_MK*dz2**3 - 114.0_MK*dz2**2 + 153.0_MK*dz2 - 48.0_MK) )
              az3 = reldx * ( c3 * (dz3 - 1.0_MK) * (25.0_MK*dz3**4 & 
                  & - 38.0_MK*dz3**3 - 3.0_MK*dz3**2 + 12.0_MK*dz3 + 12.0_MK) )
              az4 = reldx * ( c3 * (dz4 - 1.0_MK) * (25.0_MK*dz4**4 & 
                  & - 38.0_MK*dz4**3 - 3.0_MK*dz4**2 + 12.0_MK*dz4 + 12.0_MK) )
              az5 = reldx * ( c2 * (dz5 - 1.0_MK) * (dz5 - 2.0_MK) &
                  & * (25.0_MK*dz5**3 - 114.0_MK*dz5**2 + 153.0_MK*dz5 - 48.0_MK) )
              az6 = reldx * ( c1 * (dz6 - 2.0_MK)* (dz6 - 3.0_MK)**3 & 
                  & * (5.0_MK * dz6 - 8.0_MK) )

!---------------------------------------------------------------------------------!
! Cancel stencils that are not supported by the mesh
!---------------------------------------------------------------------------------!
              IF(bc_to(1) .EQ. 0)THEN
                IF( i1 .LT. 1-nghost_to(1) ) ax1 = 0.0_MK
                IF( i1 .LT. 1-nghost_to(1) ) i1  = 1-nghost_to(1)
                IF( i2 .LT. 1-nghost_to(1) ) ax2 = 0.0_MK
                IF( i2 .LT. 1-nghost_to(1) ) i2  = 1-nghost_to(1)
                IF( i3 .LT. 1-nghost_to(1) ) ax3 = 0.0_MK
                IF( i3 .LT. 1-nghost_to(1) ) i3  = 1-nghost_to(1)
                IF( i4 .LT. 1-nghost_to(1) ) ax4 = 0.0_MK
                IF( i4 .LT. 1-nghost_to(1) ) i4  = 1-nghost_to(1)
                IF( i5 .LT. 1-nghost_to(1) ) ax5 = 0.0_MK
                IF( i5 .LT. 1-nghost_to(1) ) i5  = 1-nghost_to(1)
                IF( i6 .LT. 1-nghost_to(1) ) ax6 = 0.0_MK
                IF( i6 .LT. 1-nghost_to(1) ) i6  = 1-nghost_to(1)
              ELSE
                IF( px .LT. xmin_to(1) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  i1  = 1-nghost_to(1)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  i2  = 1-nghost_to(1)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  i3  = 1-nghost_to(1)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  i4  = 1-nghost_to(1)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  i5  = 1-nghost_to(1)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  i6  = 1-nghost_to(1)
                END IF
              END IF


              IF(bc_to(2) .EQ. 0)THEN
                IF( i1 .GT. ncell_to(1)+nghost_to(2) ) ax1 = 0.0_MK
                IF( i1 .GT. ncell_to(1)+nghost_to(2) ) i1=ncell_to(1)+nghost_to(2)
                IF( i2 .GT. ncell_to(1)+nghost_to(2) ) ax2 = 0.0_MK
                IF( i2 .GT. ncell_to(1)+nghost_to(2) ) i2=ncell_to(1)+nghost_to(2)
                IF( i3 .GT. ncell_to(1)+nghost_to(2) ) ax3 = 0.0_MK
                IF( i3 .GT. ncell_to(1)+nghost_to(2) ) i3=ncell_to(1)+nghost_to(2)
                IF( i4 .GT. ncell_to(1)+nghost_to(2) ) ax4 = 0.0_MK
                IF( i4 .GT. ncell_to(1)+nghost_to(2) ) i4=ncell_to(1)+nghost_to(2)
                IF( i5 .GT. ncell_to(1)+nghost_to(2) ) ax5 = 0.0_MK
                IF( i5 .GT. ncell_to(1)+nghost_to(2) ) i5=ncell_to(1)+nghost_to(2)
                IF( i6 .GT. ncell_to(1)+nghost_to(2) ) ax6 = 0.0_MK
                IF( i6 .GT. ncell_to(1)+nghost_to(2) ) i6=ncell_to(1)+nghost_to(2)
              ELSE
                IF( px .GT. xmax_to(1) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  i1  = ncell_to(1)+nghost_to(2)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  i2  = ncell_to(1)+nghost_to(2)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  i3  = ncell_to(1)+nghost_to(2)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  i4  = ncell_to(1)+nghost_to(2)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  i5  = ncell_to(1)+nghost_to(2)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  i6  = ncell_to(1)+nghost_to(2)
                END IF
              END IF

              IF(bc_to(3) .EQ. 0)THEN
                IF( j1 .LT. 1-nghost_to(3) ) ay1 = 0.0_MK
                IF( j1 .LT. 1-nghost_to(3) ) j1  = 1-nghost_to(3)
                IF( j2 .LT. 1-nghost_to(3) ) ay2 = 0.0_MK
                IF( j2 .LT. 1-nghost_to(3) ) j2  = 1-nghost_to(3)
                IF( j3 .LT. 1-nghost_to(3) ) ay3 = 0.0_MK
                IF( j3 .LT. 1-nghost_to(3) ) j3  = 1-nghost_to(3)
                IF( j4 .LT. 1-nghost_to(3) ) ay4 = 0.0_MK
                IF( j4 .LT. 1-nghost_to(3) ) j4  = 1-nghost_to(3)
                IF( j5 .LT. 1-nghost_to(3) ) ay5 = 0.0_MK
                IF( j5 .LT. 1-nghost_to(3) ) j5  = 1-nghost_to(3)
                IF( j6 .LT. 1-nghost_to(3) ) ay6 = 0.0_MK
                IF( j6 .LT. 1-nghost_to(3) ) j6  = 1-nghost_to(3)
              ELSE
                IF( py .LT. xmin_to(2) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  j1  = 1-nghost_to(3)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  j2  = 1-nghost_to(3)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  j3  = 1-nghost_to(3)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  j4  = 1-nghost_to(3)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  j5  = 1-nghost_to(3)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  j6  = 1-nghost_to(3)
                END IF
              END IF

              IF(bc_to(4) .EQ. 0)THEN
                IF( j1 .GT. ncell_to(2)+nghost_to(4) ) ay1 = 0.0_MK
                IF( j1 .GT. ncell_to(2)+nghost_to(4) ) j1=ncell_to(2)+nghost_to(4)
                IF( j2 .GT. ncell_to(2)+nghost_to(4) ) ay2 = 0.0_MK
                IF( j2 .GT. ncell_to(2)+nghost_to(4) ) j2=ncell_to(2)+nghost_to(4)
                IF( j3 .GT. ncell_to(2)+nghost_to(4) ) ay3 = 0.0_MK
                IF( j3 .GT. ncell_to(2)+nghost_to(4) ) j3=ncell_to(2)+nghost_to(4)
                IF( j4 .GT. ncell_to(2)+nghost_to(4) ) ay4 = 0.0_MK
                IF( j4 .GT. ncell_to(2)+nghost_to(4) ) j4=ncell_to(2)+nghost_to(4)
                IF( j5 .GT. ncell_to(2)+nghost_to(4) ) ay5 = 0.0_MK
                IF( j5 .GT. ncell_to(2)+nghost_to(4) ) j5=ncell_to(2)+nghost_to(4)
                IF( j6 .GT. ncell_to(2)+nghost_to(4) ) ay6 = 0.0_MK
                IF( j6 .GT. ncell_to(2)+nghost_to(4) ) j6=ncell_to(2)+nghost_to(4)
              ELSE
                IF( py .GT. xmax_to(2) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  j1  = ncell_to(2)+nghost_to(4)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  j2  = ncell_to(2)+nghost_to(4)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  j3  = ncell_to(2)+nghost_to(4)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  j4  = ncell_to(2)+nghost_to(4)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  j5  = ncell_to(2)+nghost_to(4)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  j6  = ncell_to(2)+nghost_to(4)
                END IF
              END IF

              IF(bc_to(5) .EQ. 0)THEN
                IF( k1 .LT. 1-nghost_to(5) ) az1 = 0.0_MK
                IF( k1 .LT. 1-nghost_to(5) ) k1  = 1-nghost_to(5)
                IF( k2 .LT. 1-nghost_to(5) ) az2 = 0.0_MK
                IF( k2 .LT. 1-nghost_to(5) ) k2  = 1-nghost_to(5)
                IF( k3 .LT. 1-nghost_to(5) ) az3 = 0.0_MK
                IF( k3 .LT. 1-nghost_to(5) ) k3  = 1-nghost_to(5)
                IF( k4 .LT. 1-nghost_to(5) ) az4 = 0.0_MK
                IF( k4 .LT. 1-nghost_to(5) ) k4  = 1-nghost_to(5)
                IF( k5 .LT. 1-nghost_to(5) ) az5 = 0.0_MK
                IF( k5 .LT. 1-nghost_to(5) ) k5  = 1-nghost_to(5)
                IF( k6 .LT. 1-nghost_to(5) ) az6 = 0.0_MK
                IF( k6 .LT. 1-nghost_to(5) ) k6  = 1-nghost_to(5)
              ELSE
                IF( pz .LT. xmin_to(3) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  k1  = 1-nghost_to(5)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  k2  = 1-nghost_to(5)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  k3  = 1-nghost_to(5)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  k4  = 1-nghost_to(5)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  k5  = 1-nghost_to(5)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  k6  = 1-nghost_to(5)
                END IF
              END IF

              IF(bc_to(6) .EQ. 0)THEN
                IF( k1 .GT. ncell_to(3)+nghost_to(6) ) az1 = 0.0_MK
                IF( k1 .GT. ncell_to(3)+nghost_to(6) ) k1=ncell_to(3)+nghost_to(6)
                IF( k2 .GT. ncell_to(3)+nghost_to(6) ) az2 = 0.0_MK
                IF( k2 .GT. ncell_to(3)+nghost_to(6) ) k2=ncell_to(3)+nghost_to(6)
                IF( k3 .GT. ncell_to(3)+nghost_to(6) ) az3 = 0.0_MK
                IF( k3 .GT. ncell_to(3)+nghost_to(6) ) k3=ncell_to(3)+nghost_to(6)
                IF( k4 .GT. ncell_to(3)+nghost_to(6) ) az4 = 0.0_MK
                IF( k4 .GT. ncell_to(3)+nghost_to(6) ) k4=ncell_to(3)+nghost_to(6)
                IF( k5 .GT. ncell_to(3)+nghost_to(6) ) az5 = 0.0_MK
                IF( k5 .GT. ncell_to(3)+nghost_to(6) ) k5=ncell_to(3)+nghost_to(6)
                IF( k6 .GT. ncell_to(3)+nghost_to(6) ) az6 = 0.0_MK
                IF( k6 .GT. ncell_to(3)+nghost_to(6) ) k6=ncell_to(3)+nghost_to(6)
              ELSE
                IF( pz .GT. xmax_to(3) )THEN
                  ax1 = 0.0_MK
                  ay1 = 0.0_MK
                  az1 = 0.0_MK
                  k1  = ncell_to(3)+nghost_to(6)
                  ax2 = 0.0_MK
                  ay2 = 0.0_MK
                  az2 = 0.0_MK
                  k2  = ncell_to(3)+nghost_to(6)
                  ax3 = 0.0_MK
                  ay3 = 0.0_MK
                  az3 = 0.0_MK
                  k3  = ncell_to(3)+nghost_to(6)
                  ax4 = 0.0_MK
                  ay4 = 0.0_MK
                  az4 = 0.0_MK
                  k4  = ncell_to(3)+nghost_to(6)
                  ax5 = 0.0_MK
                  ay5 = 0.0_MK
                  az5 = 0.0_MK
                  k5  = ncell_to(3)+nghost_to(6)
                  ax6 = 0.0_MK
                  ay6 = 0.0_MK
                  az6 = 0.0_MK
                  k6  = ncell_to(3)+nghost_to(6)
                END IF
              END IF

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 3D kernel and apply to the field
!---------------------------------------------------------------------------------!
              mesh_to(1:nvar,i1,j1,k1) = mesh_to(1:nvar,i1,j1,k1) &
                                       & + ax1*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k1) = mesh_to(1:nvar,i2,j1,k1) &
                                       & + ax2*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k1) = mesh_to(1:nvar,i3,j1,k1) &
                                       & + ax3*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k1) = mesh_to(1:nvar,i4,j1,k1) &
                                       & + ax4*ay1*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k1) = mesh_to(1:nvar,i5,j1,k1) &
                                       & + ax5*ay1*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k1) = mesh_to(1:nvar,i6,j1,k1) &
                                       & + ax6*ay1*az1*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j2,k1) = mesh_to(1:nvar,i1,j2,k1) &
                                       & + ax1*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k1) = mesh_to(1:nvar,i2,j2,k1) &
                                       & + ax2*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k1) = mesh_to(1:nvar,i3,j2,k1) &
                                       & + ax3*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k1) = mesh_to(1:nvar,i4,j2,k1) &
                                       & + ax4*ay2*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k1) = mesh_to(1:nvar,i5,j2,k1) &
                                       & + ax5*ay2*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k1) = mesh_to(1:nvar,i6,j2,k1) &
                                       & + ax6*ay2*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k1) = mesh_to(1:nvar,i1,j3,k1) &
                                       & + ax1*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k1) = mesh_to(1:nvar,i2,j3,k1) &
                                       & + ax2*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k1) = mesh_to(1:nvar,i3,j3,k1) &
                                       & + ax3*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k1) = mesh_to(1:nvar,i4,j3,k1) &
                                       & + ax4*ay3*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k1) = mesh_to(1:nvar,i5,j3,k1) &
                                       & + ax5*ay3*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k1) = mesh_to(1:nvar,i6,j3,k1) &
                                       & + ax6*ay3*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k1) = mesh_to(1:nvar,i1,j4,k1) &
                                       & + ax1*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k1) = mesh_to(1:nvar,i2,j4,k1) &
                                       & + ax2*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k1) = mesh_to(1:nvar,i3,j4,k1) &
                                       & + ax3*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k1) = mesh_to(1:nvar,i4,j4,k1) &
                                       & + ax4*ay4*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k1) = mesh_to(1:nvar,i5,j4,k1) &
                                       & + ax5*ay4*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k1) = mesh_to(1:nvar,i6,j4,k1) &
                                       & + ax6*ay4*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j5,k1) = mesh_to(1:nvar,i1,j5,k1) &
                                       & + ax1*ay5*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k1) = mesh_to(1:nvar,i2,j5,k1) &
                                       & + ax2*ay5*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k1) = mesh_to(1:nvar,i3,j5,k1) &
                                       & + ax3*ay5*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k1) = mesh_to(1:nvar,i4,j5,k1) &
                                       & + ax4*ay5*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k1) = mesh_to(1:nvar,i5,j5,k1) &
                                       & + ax5*ay5*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k1) = mesh_to(1:nvar,i6,j5,k1) &
                                       & + ax6*ay5*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j6,k1) = mesh_to(1:nvar,i1,j6,k1) &
                                       & + ax1*ay6*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k1) = mesh_to(1:nvar,i2,j6,k1) &
                                       & + ax2*ay6*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k1) = mesh_to(1:nvar,i3,j6,k1) &
                                       & + ax3*ay6*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k1) = mesh_to(1:nvar,i4,j6,k1) &
                                       & + ax4*ay6*az1*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k1) = mesh_to(1:nvar,i5,j6,k1) &
                                       & + ax5*ay6*az1*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k1) = mesh_to(1:nvar,i6,j6,k1) &
                                       & + ax6*ay6*az1*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k2) = mesh_to(1:nvar,i1,j1,k2) &
                                       & + ax1*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k2) = mesh_to(1:nvar,i2,j1,k2) &
                                       & + ax2*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k2) = mesh_to(1:nvar,i3,j1,k2) &
                                       & + ax3*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k2) = mesh_to(1:nvar,i4,j1,k2) &
                                       & + ax4*ay1*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k2) = mesh_to(1:nvar,i5,j1,k2) &
                                       & + ax5*ay1*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k2) = mesh_to(1:nvar,i6,j1,k2) &
                                       & + ax6*ay1*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j2,k2) = mesh_to(1:nvar,i1,j2,k2) &
                                       & + ax1*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k2) = mesh_to(1:nvar,i2,j2,k2) &
                                       & + ax2*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k2) = mesh_to(1:nvar,i3,j2,k2) &
                                       & + ax3*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k2) = mesh_to(1:nvar,i4,j2,k2) &
                                       & + ax4*ay2*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k2) = mesh_to(1:nvar,i5,j2,k2) &
                                       & + ax5*ay2*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k2) = mesh_to(1:nvar,i6,j2,k2) &
                                       & + ax6*ay2*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j3,k2) = mesh_to(1:nvar,i1,j3,k2) &
                                       & + ax1*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k2) = mesh_to(1:nvar,i2,j3,k2) &
                                       & + ax2*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k2) = mesh_to(1:nvar,i3,j3,k2) &
                                       & + ax3*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k2) = mesh_to(1:nvar,i4,j3,k2) &
                                       & + ax4*ay3*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k2) = mesh_to(1:nvar,i5,j3,k2) &
                                       & + ax5*ay3*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k2) = mesh_to(1:nvar,i6,j3,k2) &
                                       & + ax6*ay3*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j4,k2) = mesh_to(1:nvar,i1,j4,k2) &
                                       & + ax1*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k2) = mesh_to(1:nvar,i2,j4,k2) &
                                       & + ax2*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k2) = mesh_to(1:nvar,i3,j4,k2) &
                                       & + ax3*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k2) = mesh_to(1:nvar,i4,j4,k2) &
                                       & + ax4*ay4*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k2) = mesh_to(1:nvar,i5,j4,k2) &
                                       & + ax5*ay4*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k2) = mesh_to(1:nvar,i6,j4,k2) &
                                       & + ax6*ay4*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j5,k2) = mesh_to(1:nvar,i1,j5,k2) &
                                       & + ax1*ay5*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k2) = mesh_to(1:nvar,i2,j5,k2) &
                                       & + ax2*ay5*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k2) = mesh_to(1:nvar,i3,j5,k2) &
                                       & + ax3*ay5*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k2) = mesh_to(1:nvar,i4,j5,k2) &
                                       & + ax4*ay5*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k2) = mesh_to(1:nvar,i5,j5,k2) &
                                       & + ax5*ay5*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k2) = mesh_to(1:nvar,i6,j5,k2) &
                                       & + ax6*ay5*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j6,k2) = mesh_to(1:nvar,i1,j6,k2) &
                                       & + ax1*ay6*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k2) = mesh_to(1:nvar,i2,j6,k2) &
                                       & + ax2*ay6*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k2) = mesh_to(1:nvar,i3,j6,k2) &
                                       & + ax3*ay6*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k2) = mesh_to(1:nvar,i4,j6,k2) &
                                       & + ax4*ay6*az2*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k2) = mesh_to(1:nvar,i5,j6,k2) &
                                       & + ax5*ay6*az2*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k2) = mesh_to(1:nvar,i6,j6,k2) &
                                       & + ax6*ay6*az2*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k3) = mesh_to(1:nvar,i1,j1,k3) &
                                       & + ax1*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k3) = mesh_to(1:nvar,i2,j1,k3) &
                                       & + ax2*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k3) = mesh_to(1:nvar,i3,j1,k3) &
                                       & + ax3*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k3) = mesh_to(1:nvar,i4,j1,k3) &
                                       & + ax4*ay1*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k3) = mesh_to(1:nvar,i5,j1,k3) &
                                       & + ax5*ay1*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k3) = mesh_to(1:nvar,i6,j1,k3) &
                                       & + ax6*ay1*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j2,k3) = mesh_to(1:nvar,i1,j2,k3) &
                                       & + ax1*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k3) = mesh_to(1:nvar,i2,j2,k3) &
                                       & + ax2*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k3) = mesh_to(1:nvar,i3,j2,k3) &
                                       & + ax3*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k3) = mesh_to(1:nvar,i4,j2,k3) &
                                       & + ax4*ay2*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k3) = mesh_to(1:nvar,i5,j2,k3) &
                                       & + ax5*ay2*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k3) = mesh_to(1:nvar,i6,j2,k3) &
                                       & + ax6*ay2*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j3,k3) = mesh_to(1:nvar,i1,j3,k3) &
                                       & + ax1*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k3) = mesh_to(1:nvar,i2,j3,k3) &
                                       & + ax2*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k3) = mesh_to(1:nvar,i3,j3,k3) &
                                       & + ax3*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k3) = mesh_to(1:nvar,i4,j3,k3) &
                                       & + ax4*ay3*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k3) = mesh_to(1:nvar,i5,j3,k3) &
                                       & + ax5*ay3*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k3) = mesh_to(1:nvar,i6,j3,k3) &
                                       & + ax6*ay3*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j4,k3) = mesh_to(1:nvar,i1,j4,k3) &
                                       & + ax1*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k3) = mesh_to(1:nvar,i2,j4,k3) &
                                       & + ax2*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k3) = mesh_to(1:nvar,i3,j4,k3) &
                                       & + ax3*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k3) = mesh_to(1:nvar,i4,j4,k3) &
                                       & + ax4*ay4*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k3) = mesh_to(1:nvar,i5,j4,k3) &
                                       & + ax5*ay4*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k3) = mesh_to(1:nvar,i6,j4,k3) &
                                       & + ax6*ay4*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j5,k3) = mesh_to(1:nvar,i1,j5,k3) &
                                       & + ax1*ay5*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k3) = mesh_to(1:nvar,i2,j5,k3) &
                                       & + ax2*ay5*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k3) = mesh_to(1:nvar,i3,j5,k3) &
                                       & + ax3*ay5*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k3) = mesh_to(1:nvar,i4,j5,k3) &
                                       & + ax4*ay5*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k3) = mesh_to(1:nvar,i5,j5,k3) &
                                       & + ax5*ay5*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k3) = mesh_to(1:nvar,i6,j5,k3) &
                                       & + ax6*ay5*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j6,k3) = mesh_to(1:nvar,i1,j6,k3) &
                                       & + ax1*ay6*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k3) = mesh_to(1:nvar,i2,j6,k3) &
                                       & + ax2*ay6*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k3) = mesh_to(1:nvar,i3,j6,k3) &
                                       & + ax3*ay6*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k3) = mesh_to(1:nvar,i4,j6,k3) &
                                       & + ax4*ay6*az3*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k3) = mesh_to(1:nvar,i5,j6,k3) &
                                       & + ax5*ay6*az3*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k3) = mesh_to(1:nvar,i6,j6,k3) &
                                       & + ax6*ay6*az3*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j1,k4) = mesh_to(1:nvar,i1,j1,k4) &
                                       & + ax1*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k4) = mesh_to(1:nvar,i2,j1,k4) &
                                       & + ax2*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k4) = mesh_to(1:nvar,i3,j1,k4) &
                                       & + ax3*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k4) = mesh_to(1:nvar,i4,j1,k4) &
                                       & + ax4*ay1*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k4) = mesh_to(1:nvar,i5,j1,k4) &
                                       & + ax5*ay1*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k4) = mesh_to(1:nvar,i6,j1,k4) &
                                       & + ax6*ay1*az4*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j2,k4) = mesh_to(1:nvar,i1,j2,k4) &
                                       & + ax1*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k4) = mesh_to(1:nvar,i2,j2,k4) &
                                       & + ax2*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k4) = mesh_to(1:nvar,i3,j2,k4) &
                                       & + ax3*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k4) = mesh_to(1:nvar,i4,j2,k4) &
                                       & + ax4*ay2*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k4) = mesh_to(1:nvar,i5,j2,k4) &
                                       & + ax5*ay2*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k4) = mesh_to(1:nvar,i6,j2,k4) &
                                       & + ax6*ay2*az4*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j3,k4) = mesh_to(1:nvar,i1,j3,k4) &
                                       & + ax1*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k4) = mesh_to(1:nvar,i2,j3,k4) &
                                       & + ax2*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k4) = mesh_to(1:nvar,i3,j3,k4) &
                                       & + ax3*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k4) = mesh_to(1:nvar,i4,j3,k4) &
                                       & + ax4*ay3*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k4) = mesh_to(1:nvar,i5,j3,k4) &
                                       & + ax5*ay3*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k4) = mesh_to(1:nvar,i6,j3,k4) &
                                       & + ax6*ay3*az4*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j4,k4) = mesh_to(1:nvar,i1,j4,k4) &
                                       & + ax1*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k4) = mesh_to(1:nvar,i2,j4,k4) &
                                       & + ax2*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k4) = mesh_to(1:nvar,i3,j4,k4) &
                                       & + ax3*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k4) = mesh_to(1:nvar,i4,j4,k4) &
                                       & + ax4*ay4*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k4) = mesh_to(1:nvar,i5,j4,k4) &
                                       & + ax5*ay4*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k4) = mesh_to(1:nvar,i6,j4,k4) &
                                       & + ax6*ay4*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j5,k4) = mesh_to(1:nvar,i1,j5,k4) &
                                       & + ax1*ay5*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k4) = mesh_to(1:nvar,i2,j5,k4) &
                                       & + ax2*ay5*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k4) = mesh_to(1:nvar,i3,j5,k4) &
                                       & + ax3*ay5*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k4) = mesh_to(1:nvar,i4,j5,k4) &
                                       & + ax4*ay5*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k4) = mesh_to(1:nvar,i5,j5,k4) &
                                       & + ax5*ay5*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k4) = mesh_to(1:nvar,i6,j5,k4) &
                                       & + ax6*ay5*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j6,k4) = mesh_to(1:nvar,i1,j6,k4) &
                                       & + ax1*ay6*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k4) = mesh_to(1:nvar,i2,j6,k4) &
                                       & + ax2*ay6*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k4) = mesh_to(1:nvar,i3,j6,k4) &
                                       & + ax3*ay6*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k4) = mesh_to(1:nvar,i4,j6,k4) &
                                       & + ax4*ay6*az4*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k4) = mesh_to(1:nvar,i5,j6,k4) &
                                       & + ax5*ay6*az4*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k4) = mesh_to(1:nvar,i6,j6,k4) &
                                       & + ax6*ay6*az4*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k5) = mesh_to(1:nvar,i1,j1,k5) &
                                       & + ax1*ay1*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k5) = mesh_to(1:nvar,i2,j1,k5) &
                                       & + ax2*ay1*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k5) = mesh_to(1:nvar,i3,j1,k5) &
                                       & + ax3*ay1*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k5) = mesh_to(1:nvar,i4,j1,k5) &
                                       & + ax4*ay1*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k5) = mesh_to(1:nvar,i5,j1,k5) &
                                       & + ax5*ay1*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k5) = mesh_to(1:nvar,i6,j1,k5) &
                                       & + ax6*ay1*az5*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j2,k5) = mesh_to(1:nvar,i1,j2,k5) &
                                       & + ax1*ay2*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k5) = mesh_to(1:nvar,i2,j2,k5) &
                                       & + ax2*ay2*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k5) = mesh_to(1:nvar,i3,j2,k5) &
                                       & + ax3*ay2*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k5) = mesh_to(1:nvar,i4,j2,k5) &
                                       & + ax4*ay2*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k5) = mesh_to(1:nvar,i5,j2,k5) &
                                       & + ax5*ay2*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k5) = mesh_to(1:nvar,i6,j2,k5) &
                                       & + ax6*ay2*az5*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j3,k5) = mesh_to(1:nvar,i1,j3,k5) &
                                       & + ax1*ay3*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k5) = mesh_to(1:nvar,i2,j3,k5) &
                                       & + ax2*ay3*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k5) = mesh_to(1:nvar,i3,j3,k5) &
                                       & + ax3*ay3*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k5) = mesh_to(1:nvar,i4,j3,k5) &
                                       & + ax4*ay3*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k5) = mesh_to(1:nvar,i5,j3,k5) &
                                       & + ax5*ay3*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k5) = mesh_to(1:nvar,i6,j3,k5) &
                                       & + ax6*ay3*az5*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j4,k5) = mesh_to(1:nvar,i1,j4,k5) &
                                       & + ax1*ay4*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k5) = mesh_to(1:nvar,i2,j4,k5) &
                                       & + ax2*ay4*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k5) = mesh_to(1:nvar,i3,j4,k5) &
                                       & + ax3*ay4*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k5) = mesh_to(1:nvar,i4,j4,k5) &
                                       & + ax4*ay4*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k5) = mesh_to(1:nvar,i5,j4,k5) &
                                       & + ax5*ay4*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k5) = mesh_to(1:nvar,i6,j4,k5) &
                                       & + ax6*ay4*az5*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j5,k5) = mesh_to(1:nvar,i1,j5,k5) &
                                       & + ax1*ay5*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k5) = mesh_to(1:nvar,i2,j5,k5) &
                                       & + ax2*ay5*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k5) = mesh_to(1:nvar,i3,j5,k5) &
                                       & + ax3*ay5*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k5) = mesh_to(1:nvar,i4,j5,k5) &
                                       & + ax4*ay5*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k5) = mesh_to(1:nvar,i5,j5,k5) &
                                       & + ax5*ay5*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k5) = mesh_to(1:nvar,i6,j5,k5) &
                                       & + ax6*ay5*az5*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j6,k5) = mesh_to(1:nvar,i1,j6,k5) &
                                       & + ax1*ay6*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k5) = mesh_to(1:nvar,i2,j6,k5) &
                                       & + ax2*ay6*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k5) = mesh_to(1:nvar,i3,j6,k5) &
                                       & + ax3*ay6*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k5) = mesh_to(1:nvar,i4,j6,k5) &
                                       & + ax4*ay6*az5*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k5) = mesh_to(1:nvar,i5,j6,k5) &
                                       & + ax5*ay6*az5*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k5) = mesh_to(1:nvar,i6,j6,k5) &
                                       & + ax6*ay6*az5*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j1,k6) = mesh_to(1:nvar,i1,j1,k6) &
                                       & + ax1*ay1*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j1,k6) = mesh_to(1:nvar,i2,j1,k6) &
                                       & + ax2*ay1*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j1,k6) = mesh_to(1:nvar,i3,j1,k6) &
                                       & + ax3*ay1*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j1,k6) = mesh_to(1:nvar,i4,j1,k6) &
                                       & + ax4*ay1*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j1,k6) = mesh_to(1:nvar,i5,j1,k6) &
                                       & + ax5*ay1*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j1,k6) = mesh_to(1:nvar,i6,j1,k6) &
                                       & + ax6*ay1*az6*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j2,k6) = mesh_to(1:nvar,i1,j2,k6) &
                                       & + ax1*ay2*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j2,k6) = mesh_to(1:nvar,i2,j2,k6) &
                                       & + ax2*ay2*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j2,k6) = mesh_to(1:nvar,i3,j2,k6) &
                                       & + ax3*ay2*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j2,k6) = mesh_to(1:nvar,i4,j2,k6) &
                                       & + ax4*ay2*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j2,k6) = mesh_to(1:nvar,i5,j2,k6) &
                                       & + ax5*ay2*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j2,k6) = mesh_to(1:nvar,i6,j2,k6) &
                                       & + ax6*ay2*az6*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j3,k6) = mesh_to(1:nvar,i1,j3,k6) &
                                       & + ax1*ay3*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j3,k6) = mesh_to(1:nvar,i2,j3,k6) &
                                       & + ax2*ay3*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j3,k6) = mesh_to(1:nvar,i3,j3,k6) &
                                       & + ax3*ay3*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j3,k6) = mesh_to(1:nvar,i4,j3,k6) &
                                       & + ax4*ay3*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j3,k6) = mesh_to(1:nvar,i5,j3,k6) &
                                       & + ax5*ay3*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j3,k6) = mesh_to(1:nvar,i6,j3,k6) &
                                       & + ax6*ay3*az6*mesh_from(1:nvar,i,j,k) 

              mesh_to(1:nvar,i1,j4,k6) = mesh_to(1:nvar,i1,j4,k6) &
                                       & + ax1*ay4*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j4,k6) = mesh_to(1:nvar,i2,j4,k6) &
                                       & + ax2*ay4*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j4,k6) = mesh_to(1:nvar,i3,j4,k6) &
                                       & + ax3*ay4*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j4,k6) = mesh_to(1:nvar,i4,j4,k6) &
                                       & + ax4*ay4*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j4,k6) = mesh_to(1:nvar,i5,j4,k6) &
                                       & + ax5*ay4*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j4,k6) = mesh_to(1:nvar,i6,j4,k6) &
                                       & + ax6*ay4*az6*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j5,k6) = mesh_to(1:nvar,i1,j5,k6) &
                                       & + ax1*ay5*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j5,k6) = mesh_to(1:nvar,i2,j5,k6) &
                                       & + ax2*ay5*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j5,k6) = mesh_to(1:nvar,i3,j5,k6) &
                                       & + ax3*ay5*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j5,k6) = mesh_to(1:nvar,i4,j5,k6) &
                                       & + ax4*ay5*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j5,k6) = mesh_to(1:nvar,i5,j5,k6) &
                                       & + ax5*ay5*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j5,k6) = mesh_to(1:nvar,i6,j5,k6) &
                                       & + ax6*ay5*az6*mesh_from(1:nvar,i,j,k)

              mesh_to(1:nvar,i1,j6,k6) = mesh_to(1:nvar,i1,j6,k6) &
                                       & + ax1*ay6*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i2,j6,k6) = mesh_to(1:nvar,i2,j6,k6) &
                                       & + ax2*ay6*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i3,j6,k6) = mesh_to(1:nvar,i3,j6,k6) &
                                       & + ax3*ay6*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i4,j6,k6) = mesh_to(1:nvar,i4,j6,k6) &
                                       & + ax4*ay6*az6*mesh_from(1:nvar,i,j,k) 
              mesh_to(1:nvar,i5,j6,k6) = mesh_to(1:nvar,i5,j6,k6) &
                                       & + ax5*ay6*az6*mesh_from(1:nvar,i,j,k)
              mesh_to(1:nvar,i6,j6,k6) = mesh_to(1:nvar,i6,j6,k6) &
                                       & + ax6*ay6*az6*mesh_from(1:nvar,i,j,k)

            END IF
          END DO ! i
        END DO ! j
      END DO ! k

    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Interpolation scheme/order unknown.')
      GOTO 9999
    END IF

!---------------------------------------------------------------------------------!
! Map the ghosts cells (overwrite) in order to get consistent values at 
! the border of the subdomains.
!---------------------------------------------------------------------------------!
    DO i = 1,pmlib_ndim
      CALL pmlib_mesh_map_ghost(topo_to,i,nvar,ierr,incl_edges = .TRUE.)
      CALL pmlib_comm_pack(mesh_to,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack(topo_to,mesh_to,1,ierr,clear = .FALSE.)
      CALL pmlib_comm_finalise(ierr)
    END DO

!---------------------------------------------------------------------------------!
! Parent to patch
!---------------------------------------------------------------------------------!
  ELSE IF(reldx .GT. 0.9999_MK)THEN

    IF(pmlib_interpolation_order .EQ. 3)THEN

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
      c_1_dx = 1.0_MK/dx_from(1)
      c_1_dy = 1.0_MK/dx_from(2)
      c_1_dz = 1.0_MK/dx_from(3)

!---------------------------------------------------------------------------------!
! Interpolate
!---------------------------------------------------------------------------------!
      DO k = 1-nghost_to(5),ncell_to(3)+nghost_to(6)
        pz = xmin_to(3) + (REAL(k - 1,MK) + 0.5_MK) *dx_to(3)

        DO j = 1-nghost_to(3),ncell_to(2)+nghost_to(4)
          py = xmin_to(2) + (REAL(j - 1,MK) + 0.5_MK) *dx_to(2)

          DO i = 1-nghost_to(1),ncell_to(1)+nghost_to(2)
            px = xmin_to(1) + (REAL(i - 1,MK) + 0.5_MK) *dx_to(1)

            IF( ( bc_to(1) .EQ. 0 .OR. i .GE. 1 ) .AND.                           &
                ( bc_to(2) .EQ. 0 .OR. i .LE. ncell_to(1) ) .AND.                 &
                ( bc_to(3) .EQ. 0 .OR. j .GE. 1 ) .AND.                           &
                ( bc_to(4) .EQ. 0 .OR. j .LE. ncell_to(2) ) .AND.                 &
                ( bc_to(5) .EQ. 0 .OR. k .GE. 1 ) .AND.                           &
                ( bc_to(6) .EQ. 0 .OR. k .LE. ncell_to(3) ) )THEN

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
              i2 = NINT( (px - xmin_from(1)) * c_1_dx )
              j2 = NINT( (py - xmin_from(2)) * c_1_dy )
              k2 = NINT( (pz - xmin_from(3)) * c_1_dz )

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
              dx2 = ( px - xmin_from(1) - dx_from(1)*(REAL(i2,MK) - 0.5_MK) ) * c_1_dx
              dy2 = ( py - xmin_from(2) - dx_from(2)*(REAL(j2,MK) - 0.5_MK) ) * c_1_dy
              dz2 = ( pz - xmin_from(3) - dx_from(3)*(REAL(k2,MK) - 0.5_MK) ) * c_1_dz

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
              IF(i1 .LT. 1-nghost_from(1)) ax1 = 0.0_MK
              IF(i1 .LT. 1-nghost_from(1)) i1  = 1-nghost_from(1)
              IF(i2 .LT. 1-nghost_from(1)) ax2 = 0.0_MK
              IF(i2 .LT. 1-nghost_from(1)) i2  = 1-nghost_from(1)
              IF(i3 .LT. 1-nghost_from(1)) ax3 = 0.0_MK
              IF(i3 .LT. 1-nghost_from(1)) i3  = 1-nghost_from(1)
              IF(i4 .LT. 1-nghost_from(1)) ax4 = 0.0_MK
              IF(i4 .LT. 1-nghost_from(1)) i4  = 1-nghost_from(1)

              IF(i1 .GT. ncell_from(1)+nghost_from(2)) ax1 = 0.0_MK
              IF(i1 .GT. ncell_from(1)+nghost_from(2))i1=ncell_from(1)+nghost_from(2)
              IF(i2 .GT. ncell_from(1)+nghost_from(2)) ax2 = 0.0_MK
              IF(i2 .GT. ncell_from(1)+nghost_from(2))i2=ncell_from(1)+nghost_from(2)
              IF(i3 .GT. ncell_from(1)+nghost_from(2)) ax3 = 0.0_MK
              IF(i3 .GT. ncell_from(1)+nghost_from(2))i3=ncell_from(1)+nghost_from(2)
              IF(i4 .GT. ncell_from(1)+nghost_from(2)) ax4 = 0.0_MK
              IF(i4 .GT. ncell_from(1)+nghost_from(2))i4=ncell_from(1)+nghost_from(2)

              IF(j1 .LT. 1-nghost_from(3)) ay1 = 0.0_MK
              IF(j1 .LT. 1-nghost_from(3)) j1  = 1-nghost_from(3)
              IF(j2 .LT. 1-nghost_from(3)) ay2 = 0.0_MK
              IF(j2 .LT. 1-nghost_from(3)) j2  = 1-nghost_from(3)
              IF(j3 .LT. 1-nghost_from(3)) ay3 = 0.0_MK
              IF(j3 .LT. 1-nghost_from(3)) j3  = 1-nghost_from(3)
              IF(j4 .LT. 1-nghost_from(3)) ay4 = 0.0_MK
              IF(j4 .LT. 1-nghost_from(3)) j4  = 1-nghost_from(3)

              IF(j1 .GT. ncell_from(2)+nghost_from(4)) ay1 = 0.0_MK
              IF(j1 .GT. ncell_from(2)+nghost_from(4))j1=ncell_from(2)+nghost_from(4)
              IF(j2 .GT. ncell_from(2)+nghost_from(4)) ay2 = 0.0_MK
              IF(j2 .GT. ncell_from(2)+nghost_from(4))j2=ncell_from(2)+nghost_from(4)
              IF(j3 .GT. ncell_from(2)+nghost_from(4)) ay3 = 0.0_MK
              IF(j3 .GT. ncell_from(2)+nghost_from(4))j3=ncell_from(2)+nghost_from(4)
              IF(j4 .GT. ncell_from(2)+nghost_from(4)) ay4 = 0.0_MK
              IF(j4 .GT. ncell_from(2)+nghost_from(4))j4=ncell_from(2)+nghost_from(4)

              IF(k1 .LT. 1-nghost_from(5) ) az1 = 0.0_MK
              IF(k1 .LT. 1-nghost_from(5) ) k1  = 1-nghost_from(5)
              IF(k2 .LT. 1-nghost_from(5) ) az2 = 0.0_MK
              IF(k2 .LT. 1-nghost_from(5) ) k2  = 1-nghost_from(5)
              IF(k3 .LT. 1-nghost_from(5) ) az3 = 0.0_MK
              IF(k3 .LT. 1-nghost_from(5) ) k3  = 1-nghost_from(5)
              IF(k4 .LT. 1-nghost_from(5) ) az4 = 0.0_MK
              IF(k4 .LT. 1-nghost_from(5) ) k4  = 1-nghost_from(5)

              IF(k1 .GT. ncell_from(3)+nghost_from(6)) az1 = 0.0_MK
              IF(k1 .GT. ncell_from(3)+nghost_from(6))k1=ncell_from(3)+nghost_from(6)
              IF(k2 .GT. ncell_from(3)+nghost_from(6)) az2 = 0.0_MK
              IF(k2 .GT. ncell_from(3)+nghost_from(6))k2=ncell_from(3)+nghost_from(6)
              IF(k3 .GT. ncell_from(3)+nghost_from(6)) az3 = 0.0_MK
              IF(k3 .GT. ncell_from(3)+nghost_from(6))k3=ncell_from(3)+nghost_from(6)
              IF(k4 .GT. ncell_from(3)+nghost_from(6)) az4 = 0.0_MK
              IF(k4 .GT. ncell_from(3)+nghost_from(6))k4=ncell_from(3)+nghost_from(6)

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 2D kernel and apply to the field
!---------------------------------------------------------------------------------!
              mesh_to(1:nvar,i,j,k) = mesh_to(1:nvar,i,j,k) + &
                                    & ax1*ay1*az1*mesh_from(1:nvar,i1,j1,k1) + &
                                    & ax2*ay1*az1*mesh_from(1:nvar,i2,j1,k1) + &
                                    & ax3*ay1*az1*mesh_from(1:nvar,i3,j1,k1) + &
                                    & ax4*ay1*az1*mesh_from(1:nvar,i4,j1,k1) + &
                                    & ax1*ay2*az1*mesh_from(1:nvar,i1,j2,k1) + &
                                    & ax2*ay2*az1*mesh_from(1:nvar,i2,j2,k1) + &
                                    & ax3*ay2*az1*mesh_from(1:nvar,i3,j2,k1) + &
                                    & ax4*ay2*az1*mesh_from(1:nvar,i4,j2,k1) + &
                                    & ax1*ay3*az1*mesh_from(1:nvar,i1,j3,k1) + &
                                    & ax2*ay3*az1*mesh_from(1:nvar,i2,j3,k1) + &
                                    & ax3*ay3*az1*mesh_from(1:nvar,i3,j3,k1) + &
                                    & ax4*ay3*az1*mesh_from(1:nvar,i4,j3,k1) + &
                                    & ax1*ay4*az1*mesh_from(1:nvar,i1,j4,k1) + &
                                    & ax2*ay4*az1*mesh_from(1:nvar,i2,j4,k1) + &
                                    & ax3*ay4*az1*mesh_from(1:nvar,i3,j4,k1) + &
                                    & ax4*ay4*az1*mesh_from(1:nvar,i4,j4,k1) + &
                                    & ax1*ay1*az2*mesh_from(1:nvar,i1,j1,k2) + &
                                    & ax2*ay1*az2*mesh_from(1:nvar,i2,j1,k2) + &
                                    & ax3*ay1*az2*mesh_from(1:nvar,i3,j1,k2) + &
                                    & ax4*ay1*az2*mesh_from(1:nvar,i4,j1,k2) + &
                                    & ax1*ay2*az2*mesh_from(1:nvar,i1,j2,k2) + &
                                    & ax2*ay2*az2*mesh_from(1:nvar,i2,j2,k2) + &
                                    & ax3*ay2*az2*mesh_from(1:nvar,i3,j2,k2) + &
                                    & ax4*ay2*az2*mesh_from(1:nvar,i4,j2,k2) + &
                                    & ax1*ay3*az2*mesh_from(1:nvar,i1,j3,k2) + &
                                    & ax2*ay3*az2*mesh_from(1:nvar,i2,j3,k2) + &
                                    & ax3*ay3*az2*mesh_from(1:nvar,i3,j3,k2) + &
                                    & ax4*ay3*az2*mesh_from(1:nvar,i4,j3,k2) + &
                                    & ax1*ay4*az2*mesh_from(1:nvar,i1,j4,k2) + &
                                    & ax2*ay4*az2*mesh_from(1:nvar,i2,j4,k2) + &
                                    & ax3*ay4*az2*mesh_from(1:nvar,i3,j4,k2) + &
                                    & ax4*ay4*az2*mesh_from(1:nvar,i4,j4,k2) + &
                                    & ax1*ay1*az3*mesh_from(1:nvar,i1,j1,k3) + &
                                    & ax2*ay1*az3*mesh_from(1:nvar,i2,j1,k3) + &
                                    & ax3*ay1*az3*mesh_from(1:nvar,i3,j1,k3) + &
                                    & ax4*ay1*az3*mesh_from(1:nvar,i4,j1,k3) + &
                                    & ax1*ay2*az3*mesh_from(1:nvar,i1,j2,k3) + &
                                    & ax2*ay2*az3*mesh_from(1:nvar,i2,j2,k3) + &
                                    & ax3*ay2*az3*mesh_from(1:nvar,i3,j2,k3) + &
                                    & ax4*ay2*az3*mesh_from(1:nvar,i4,j2,k3) + &
                                    & ax1*ay3*az3*mesh_from(1:nvar,i1,j3,k3) + &
                                    & ax2*ay3*az3*mesh_from(1:nvar,i2,j3,k3) + &
                                    & ax3*ay3*az3*mesh_from(1:nvar,i3,j3,k3) + &
                                    & ax4*ay3*az3*mesh_from(1:nvar,i4,j3,k3) + &
                                    & ax1*ay4*az3*mesh_from(1:nvar,i1,j4,k3) + &
                                    & ax2*ay4*az3*mesh_from(1:nvar,i2,j4,k3) + &
                                    & ax3*ay4*az3*mesh_from(1:nvar,i3,j4,k3) + &
                                    & ax4*ay4*az3*mesh_from(1:nvar,i4,j4,k3) + &
                                    & ax1*ay1*az4*mesh_from(1:nvar,i1,j1,k4) + &
                                    & ax2*ay1*az4*mesh_from(1:nvar,i2,j1,k4) + &
                                    & ax3*ay1*az4*mesh_from(1:nvar,i3,j1,k4) + &
                                    & ax4*ay1*az4*mesh_from(1:nvar,i4,j1,k4) + &
                                    & ax1*ay2*az4*mesh_from(1:nvar,i1,j2,k4) + &
                                    & ax2*ay2*az4*mesh_from(1:nvar,i2,j2,k4) + &
                                    & ax3*ay2*az4*mesh_from(1:nvar,i3,j2,k4) + &
                                    & ax4*ay2*az4*mesh_from(1:nvar,i4,j2,k4) + &
                                    & ax1*ay3*az4*mesh_from(1:nvar,i1,j3,k4) + &
                                    & ax2*ay3*az4*mesh_from(1:nvar,i2,j3,k4) + &
                                    & ax3*ay3*az4*mesh_from(1:nvar,i3,j3,k4) + &
                                    & ax4*ay3*az4*mesh_from(1:nvar,i4,j3,k4) + &
                                    & ax1*ay4*az4*mesh_from(1:nvar,i1,j4,k4) + &
                                    & ax2*ay4*az4*mesh_from(1:nvar,i2,j4,k4) + &
                                    & ax3*ay4*az4*mesh_from(1:nvar,i3,j4,k4) + &
                                    & ax4*ay4*az4*mesh_from(1:nvar,i4,j4,k4)

            END IF
          END DO ! i
        END DO ! j
      END DO ! k

    ELSEIF(pmlib_interpolation_order .EQ. 4)THEN
!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
      c1 = -1.0_MK/24.0_MK
      c2 =  1.0_MK/24.0_MK 
      c3 = -1.0_MK/12.0_MK

      c_1_dx = 1.0_MK/dx_from(1)
      c_1_dy = 1.0_MK/dx_from(2)
      c_1_dz = 1.0_MK/dx_from(3)

!---------------------------------------------------------------------------------!
! Interpolate
!---------------------------------------------------------------------------------!
      DO k = 1-nghost_to(5),ncell_to(3)+nghost_to(6)
        pz = xmin_to(3) + (REAL(k - 1,MK) + 0.5_MK) *dx_to(3)

        DO j = 1-nghost_to(3),ncell_to(2)+nghost_to(4)
          py = xmin_to(2) + (REAL(j - 1,MK) + 0.5_MK) *dx_to(2)

          DO i = 1-nghost_to(1),ncell_to(1)+nghost_to(2)
            px = xmin_to(1) + (REAL(i - 1,MK) + 0.5_MK) *dx_to(1)

            IF( ( bc_to(1) .EQ. 0 .OR. i .GE. 1 ) .AND.                           &
                ( bc_to(2) .EQ. 0 .OR. i .LE. ncell_to(1) ) .AND.                 &
                ( bc_to(3) .EQ. 0 .OR. j .GE. 1 ) .AND.                           &
                ( bc_to(4) .EQ. 0 .OR. j .LE. ncell_to(2) ) .AND.                 &
                ( bc_to(5) .EQ. 0 .OR. k .GE. 1 ) .AND.                           &
                ( bc_to(6) .EQ. 0 .OR. k .LE. ncell_to(3) ) )THEN

!---------------------------------------------------------------------------------!
! Find index of south-west-bottom cell
!---------------------------------------------------------------------------------!
              i3 = NINT( (px - xmin_from(1)) * c_1_dx )
              j3 = NINT( (py - xmin_from(2)) * c_1_dy )
              k3 = NINT( (pz - xmin_from(3)) * c_1_dz )

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
              dx3 = (px - xmin_from(1) - dx_from(1)*(REAL(i3,MK) - 0.5_MK))*c_1_dx
              dy3 = (py - xmin_from(2) - dx_from(2)*(REAL(j3,MK) - 0.5_MK))*c_1_dy
              dz3 = (pz - xmin_from(3) - dx_from(3)*(REAL(k3,MK) - 0.5_MK))*c_1_dz

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
              IF(i1 .LT. 1-nghost_from(1)) ax1 = 0.0_MK
              IF(i1 .LT. 1-nghost_from(1)) i1  = 1-nghost_from(1)
              IF(i2 .LT. 1-nghost_from(1)) ax2 = 0.0_MK
              IF(i2 .LT. 1-nghost_from(1)) i2  = 1-nghost_from(1)
              IF(i3 .LT. 1-nghost_from(1)) ax3 = 0.0_MK
              IF(i3 .LT. 1-nghost_from(1)) i3  = 1-nghost_from(1)
              IF(i4 .LT. 1-nghost_from(1)) ax4 = 0.0_MK
              IF(i4 .LT. 1-nghost_from(1)) i4  = 1-nghost_from(1)
              IF(i5 .LT. 1-nghost_from(1)) ax5 = 0.0_MK
              IF(i5 .LT. 1-nghost_from(1)) i5  = 1-nghost_from(1)
              IF(i6 .LT. 1-nghost_from(1)) ax6 = 0.0_MK
              IF(i6 .LT. 1-nghost_from(1)) i6  = 1-nghost_from(1)

              IF(i1 .GT. ncell_from(1)+nghost_from(2)) ax1 = 0.0_MK
              IF(i1 .GT. ncell_from(1)+nghost_from(2))i1=ncell_from(1)+nghost_from(2)
              IF(i2 .GT. ncell_from(1)+nghost_from(2)) ax2 = 0.0_MK
              IF(i2 .GT. ncell_from(1)+nghost_from(2))i2=ncell_from(1)+nghost_from(2)
              IF(i3 .GT. ncell_from(1)+nghost_from(2)) ax3 = 0.0_MK
              IF(i3 .GT. ncell_from(1)+nghost_from(2))i3=ncell_from(1)+nghost_from(2)
              IF(i4 .GT. ncell_from(1)+nghost_from(2)) ax4 = 0.0_MK
              IF(i4 .GT. ncell_from(1)+nghost_from(2))i4=ncell_from(1)+nghost_from(2)
              IF(i5 .GT. ncell_from(1)+nghost_from(2)) ax5 = 0.0_MK
              IF(i5 .GT. ncell_from(1)+nghost_from(2))i5=ncell_from(1)+nghost_from(2)
              IF(i6 .GT. ncell_from(1)+nghost_from(2)) ax6 = 0.0_MK
              IF(i6 .GT. ncell_from(1)+nghost_from(2))i6=ncell_from(1)+nghost_from(2)

              IF(j1 .LT. 1-nghost_from(3)) ay1 = 0.0_MK
              IF(j1 .LT. 1-nghost_from(3)) j1  = 1-nghost_from(3)
              IF(j2 .LT. 1-nghost_from(3)) ay2 = 0.0_MK
              IF(j2 .LT. 1-nghost_from(3)) j2  = 1-nghost_from(3)
              IF(j3 .LT. 1-nghost_from(3)) ay3 = 0.0_MK
              IF(j3 .LT. 1-nghost_from(3)) j3  = 1-nghost_from(3)
              IF(j4 .LT. 1-nghost_from(3)) ay4 = 0.0_MK
              IF(j4 .LT. 1-nghost_from(3)) j4  = 1-nghost_from(3)
              IF(j5 .LT. 1-nghost_from(3)) ay5 = 0.0_MK
              IF(j5 .LT. 1-nghost_from(3)) j5  = 1-nghost_from(3)
              IF(j6 .LT. 1-nghost_from(3)) ay6 = 0.0_MK
              IF(j6 .LT. 1-nghost_from(3)) j6  = 1-nghost_from(3)

              IF(j1 .GT. ncell_from(2)+nghost_from(4)) ay1 = 0.0_MK
              IF(j1 .GT. ncell_from(2)+nghost_from(4))j1=ncell_from(2)+nghost_from(4)
              IF(j2 .GT. ncell_from(2)+nghost_from(4)) ay2 = 0.0_MK
              IF(j2 .GT. ncell_from(2)+nghost_from(4))j2=ncell_from(2)+nghost_from(4)
              IF(j3 .GT. ncell_from(2)+nghost_from(4)) ay3 = 0.0_MK
              IF(j3 .GT. ncell_from(2)+nghost_from(4))j3=ncell_from(2)+nghost_from(4)
              IF(j4 .GT. ncell_from(2)+nghost_from(4)) ay4 = 0.0_MK
              IF(j4 .GT. ncell_from(2)+nghost_from(4))j4=ncell_from(2)+nghost_from(4)
              IF(j5 .GT. ncell_from(2)+nghost_from(4)) ay5 = 0.0_MK
              IF(j5 .GT. ncell_from(2)+nghost_from(4))j5=ncell_from(2)+nghost_from(4)
              IF(j6 .GT. ncell_from(2)+nghost_from(4)) ay6 = 0.0_MK
              IF(j6 .GT. ncell_from(2)+nghost_from(4))j6=ncell_from(2)+nghost_from(4)

              IF(k1 .LT. 1-nghost_from(5)) az1 = 0.0_MK
              IF(k1 .LT. 1-nghost_from(5)) k1  = 1-nghost_from(5)
              IF(k2 .LT. 1-nghost_from(5)) az2 = 0.0_MK
              IF(k2 .LT. 1-nghost_from(5)) k2  = 1-nghost_from(5)
              IF(k3 .LT. 1-nghost_from(5)) az3 = 0.0_MK
              IF(k3 .LT. 1-nghost_from(5)) k3  = 1-nghost_from(5)
              IF(k4 .LT. 1-nghost_from(5)) az4 = 0.0_MK
              IF(k4 .LT. 1-nghost_from(5)) k4  = 1-nghost_from(5)
              IF(k5 .LT. 1-nghost_from(5)) az5 = 0.0_MK
              IF(k5 .LT. 1-nghost_from(5)) k5  = 1-nghost_from(5)
              IF(k6 .LT. 1-nghost_from(5)) az6 = 0.0_MK
              IF(k6 .LT. 1-nghost_from(5)) k6  = 1-nghost_from(5)

              IF(k1 .GT. ncell_from(3)+nghost_from(6)) az1 = 0.0_MK
              IF(k1 .GT. ncell_from(3)+nghost_from(6))k1=ncell_from(3)+nghost_from(6)
              IF(k2 .GT. ncell_from(3)+nghost_from(6)) az2 = 0.0_MK
              IF(k2 .GT. ncell_from(3)+nghost_from(6))k2=ncell_from(3)+nghost_from(6)
              IF(k3 .GT. ncell_from(3)+nghost_from(6)) az3 = 0.0_MK
              IF(k3 .GT. ncell_from(3)+nghost_from(6))k3=ncell_from(3)+nghost_from(6)
              IF(k4 .GT. ncell_from(3)+nghost_from(6)) az4 = 0.0_MK
              IF(k4 .GT. ncell_from(3)+nghost_from(6))k4=ncell_from(3)+nghost_from(6)
              IF(k5 .GT. ncell_from(3)+nghost_from(6)) az5 = 0.0_MK
              IF(k5 .GT. ncell_from(3)+nghost_from(6))k5=ncell_from(3)+nghost_from(6)
              IF(k6 .GT. ncell_from(3)+nghost_from(6)) az6 = 0.0_MK
              IF(k6 .GT. ncell_from(3)+nghost_from(6))k6=ncell_from(3)+nghost_from(6)

!---------------------------------------------------------------------------------!
! Combine 1D kernels into the 2D kernel and apply to the field
!---------------------------------------------------------------------------------!
              mesh_to(1:nvar,i,j,k) = mesh_to(1:nvar,i,j,k) +                  &
                                    & ax1*ay1*az1*mesh_from(1:nvar,i1,j1,k1) + &
                                    & ax2*ay1*az1*mesh_from(1:nvar,i2,j1,k1) + &
                                    & ax3*ay1*az1*mesh_from(1:nvar,i3,j1,k1) + &
                                    & ax4*ay1*az1*mesh_from(1:nvar,i4,j1,k1) + &
                                    & ax5*ay1*az1*mesh_from(1:nvar,i5,j1,k1) + &
                                    & ax6*ay1*az1*mesh_from(1:nvar,i6,j1,k1) + &
                                    & ax1*ay2*az1*mesh_from(1:nvar,i1,j2,k1) + &
                                    & ax2*ay2*az1*mesh_from(1:nvar,i2,j2,k1) + &
                                    & ax3*ay2*az1*mesh_from(1:nvar,i3,j2,k1) + &
                                    & ax4*ay2*az1*mesh_from(1:nvar,i4,j2,k1) + &
                                    & ax5*ay2*az1*mesh_from(1:nvar,i5,j2,k1) + &
                                    & ax6*ay2*az1*mesh_from(1:nvar,i6,j2,k1) + &
                                    & ax1*ay3*az1*mesh_from(1:nvar,i1,j3,k1) + &
                                    & ax2*ay3*az1*mesh_from(1:nvar,i2,j3,k1) + &
                                    & ax3*ay3*az1*mesh_from(1:nvar,i3,j3,k1) + &
                                    & ax4*ay3*az1*mesh_from(1:nvar,i4,j3,k1) + &
                                    & ax5*ay3*az1*mesh_from(1:nvar,i5,j3,k1) + &
                                    & ax6*ay3*az1*mesh_from(1:nvar,i6,j3,k1) + &
                                    & ax1*ay4*az1*mesh_from(1:nvar,i1,j4,k1) + &
                                    & ax2*ay4*az1*mesh_from(1:nvar,i2,j4,k1) + &
                                    & ax3*ay4*az1*mesh_from(1:nvar,i3,j4,k1) + &
                                    & ax4*ay4*az1*mesh_from(1:nvar,i4,j4,k1) + &
                                    & ax5*ay4*az1*mesh_from(1:nvar,i5,j4,k1) + &
                                    & ax6*ay4*az1*mesh_from(1:nvar,i6,j4,k1) + &
                                    & ax1*ay5*az1*mesh_from(1:nvar,i1,j5,k1) + &
                                    & ax2*ay5*az1*mesh_from(1:nvar,i2,j5,k1) + &
                                    & ax3*ay5*az1*mesh_from(1:nvar,i3,j5,k1) + &
                                    & ax4*ay5*az1*mesh_from(1:nvar,i4,j5,k1) + &
                                    & ax5*ay5*az1*mesh_from(1:nvar,i5,j5,k1) + &
                                    & ax6*ay5*az1*mesh_from(1:nvar,i6,j5,k1) + &
                                    & ax1*ay6*az1*mesh_from(1:nvar,i1,j6,k1) + &
                                    & ax2*ay6*az1*mesh_from(1:nvar,i2,j6,k1) + &
                                    & ax3*ay6*az1*mesh_from(1:nvar,i3,j6,k1) + &
                                    & ax4*ay6*az1*mesh_from(1:nvar,i4,j6,k1) + &
                                    & ax5*ay6*az1*mesh_from(1:nvar,i5,j6,k1) + &
                                    & ax6*ay6*az1*mesh_from(1:nvar,i6,j6,k1) + &
                                    & ax1*ay1*az2*mesh_from(1:nvar,i1,j1,k2) + &
                                    & ax2*ay1*az2*mesh_from(1:nvar,i2,j1,k2) + &
                                    & ax3*ay1*az2*mesh_from(1:nvar,i3,j1,k2) + &
                                    & ax4*ay1*az2*mesh_from(1:nvar,i4,j1,k2) + &
                                    & ax5*ay1*az2*mesh_from(1:nvar,i5,j1,k2) + &
                                    & ax6*ay1*az2*mesh_from(1:nvar,i6,j1,k2) + &
                                    & ax1*ay2*az2*mesh_from(1:nvar,i1,j2,k2) + &
                                    & ax2*ay2*az2*mesh_from(1:nvar,i2,j2,k2) + &
                                    & ax3*ay2*az2*mesh_from(1:nvar,i3,j2,k2) + &
                                    & ax4*ay2*az2*mesh_from(1:nvar,i4,j2,k2) + &
                                    & ax5*ay2*az2*mesh_from(1:nvar,i5,j2,k2) + &
                                    & ax6*ay2*az2*mesh_from(1:nvar,i6,j2,k2) + &
                                    & ax1*ay3*az2*mesh_from(1:nvar,i1,j3,k2) + &
                                    & ax2*ay3*az2*mesh_from(1:nvar,i2,j3,k2) + &
                                    & ax3*ay3*az2*mesh_from(1:nvar,i3,j3,k2) + &
                                    & ax4*ay3*az2*mesh_from(1:nvar,i4,j3,k2) + &
                                    & ax5*ay3*az2*mesh_from(1:nvar,i5,j3,k2) + &
                                    & ax6*ay3*az2*mesh_from(1:nvar,i6,j3,k2) + &
                                    & ax1*ay4*az2*mesh_from(1:nvar,i1,j4,k2) + &
                                    & ax2*ay4*az2*mesh_from(1:nvar,i2,j4,k2) + &
                                    & ax3*ay4*az2*mesh_from(1:nvar,i3,j4,k2) + &
                                    & ax4*ay4*az2*mesh_from(1:nvar,i4,j4,k2) + &
                                    & ax5*ay4*az2*mesh_from(1:nvar,i5,j4,k2) + &
                                    & ax6*ay4*az2*mesh_from(1:nvar,i6,j4,k2) + &
                                    & ax1*ay5*az2*mesh_from(1:nvar,i1,j5,k2) + &
                                    & ax2*ay5*az2*mesh_from(1:nvar,i2,j5,k2) + &
                                    & ax3*ay5*az2*mesh_from(1:nvar,i3,j5,k2) + &
                                    & ax4*ay5*az2*mesh_from(1:nvar,i4,j5,k2) + &
                                    & ax5*ay5*az2*mesh_from(1:nvar,i5,j5,k2) + &
                                    & ax6*ay5*az2*mesh_from(1:nvar,i6,j5,k2) + &
                                    & ax1*ay6*az2*mesh_from(1:nvar,i1,j6,k2) + &
                                    & ax2*ay6*az2*mesh_from(1:nvar,i2,j6,k2) + &
                                    & ax3*ay6*az2*mesh_from(1:nvar,i3,j6,k2) + &
                                    & ax4*ay6*az2*mesh_from(1:nvar,i4,j6,k2) + &
                                    & ax5*ay6*az2*mesh_from(1:nvar,i5,j6,k2) + &
                                    & ax6*ay6*az2*mesh_from(1:nvar,i6,j6,k2) + &
                                    & ax1*ay1*az3*mesh_from(1:nvar,i1,j1,k3) + &
                                    & ax2*ay1*az3*mesh_from(1:nvar,i2,j1,k3) + &
                                    & ax3*ay1*az3*mesh_from(1:nvar,i3,j1,k3) + &
                                    & ax4*ay1*az3*mesh_from(1:nvar,i4,j1,k3) + &
                                    & ax5*ay1*az3*mesh_from(1:nvar,i5,j1,k3) + &
                                    & ax6*ay1*az3*mesh_from(1:nvar,i6,j1,k3) + &
                                    & ax1*ay2*az3*mesh_from(1:nvar,i1,j2,k3) + &
                                    & ax2*ay2*az3*mesh_from(1:nvar,i2,j2,k3) + &
                                    & ax3*ay2*az3*mesh_from(1:nvar,i3,j2,k3) + &
                                    & ax4*ay2*az3*mesh_from(1:nvar,i4,j2,k3) + &
                                    & ax5*ay2*az3*mesh_from(1:nvar,i5,j2,k3) + &
                                    & ax6*ay2*az3*mesh_from(1:nvar,i6,j2,k3) + &
                                    & ax1*ay3*az3*mesh_from(1:nvar,i1,j3,k3) + &
                                    & ax2*ay3*az3*mesh_from(1:nvar,i2,j3,k3) + &
                                    & ax3*ay3*az3*mesh_from(1:nvar,i3,j3,k3) + &
                                    & ax4*ay3*az3*mesh_from(1:nvar,i4,j3,k3) + &
                                    & ax5*ay3*az3*mesh_from(1:nvar,i5,j3,k3) + &
                                    & ax6*ay3*az3*mesh_from(1:nvar,i6,j3,k3) + &
                                    & ax1*ay4*az3*mesh_from(1:nvar,i1,j4,k3) + &
                                    & ax2*ay4*az3*mesh_from(1:nvar,i2,j4,k3) + &
                                    & ax3*ay4*az3*mesh_from(1:nvar,i3,j4,k3) + &
                                    & ax4*ay4*az3*mesh_from(1:nvar,i4,j4,k3) + &
                                    & ax5*ay4*az3*mesh_from(1:nvar,i5,j4,k3) + &
                                    & ax6*ay4*az3*mesh_from(1:nvar,i6,j4,k3) + &
                                    & ax1*ay5*az3*mesh_from(1:nvar,i1,j5,k3) + &
                                    & ax2*ay5*az3*mesh_from(1:nvar,i2,j5,k3) + &
                                    & ax3*ay5*az3*mesh_from(1:nvar,i3,j5,k3) + &
                                    & ax4*ay5*az3*mesh_from(1:nvar,i4,j5,k3) + &
                                    & ax5*ay5*az3*mesh_from(1:nvar,i5,j5,k3) + &
                                    & ax6*ay5*az3*mesh_from(1:nvar,i6,j5,k3) + &
                                    & ax1*ay6*az3*mesh_from(1:nvar,i1,j6,k3) + &
                                    & ax2*ay6*az3*mesh_from(1:nvar,i2,j6,k3) + &
                                    & ax3*ay6*az3*mesh_from(1:nvar,i3,j6,k3) + &
                                    & ax4*ay6*az3*mesh_from(1:nvar,i4,j6,k3) + &
                                    & ax5*ay6*az3*mesh_from(1:nvar,i5,j6,k3) + &
                                    & ax6*ay6*az3*mesh_from(1:nvar,i6,j6,k3) + &
                                    & ax1*ay1*az4*mesh_from(1:nvar,i1,j1,k4) + &
                                    & ax2*ay1*az4*mesh_from(1:nvar,i2,j1,k4) + &
                                    & ax3*ay1*az4*mesh_from(1:nvar,i3,j1,k4) + &
                                    & ax4*ay1*az4*mesh_from(1:nvar,i4,j1,k4) + &
                                    & ax5*ay1*az4*mesh_from(1:nvar,i5,j1,k4) + &
                                    & ax6*ay1*az4*mesh_from(1:nvar,i6,j1,k4) + &
                                    & ax1*ay2*az4*mesh_from(1:nvar,i1,j2,k4) + &
                                    & ax2*ay2*az4*mesh_from(1:nvar,i2,j2,k4) + &
                                    & ax3*ay2*az4*mesh_from(1:nvar,i3,j2,k4) + &
                                    & ax4*ay2*az4*mesh_from(1:nvar,i4,j2,k4) + &
                                    & ax5*ay2*az4*mesh_from(1:nvar,i5,j2,k4) + &
                                    & ax6*ay2*az4*mesh_from(1:nvar,i6,j2,k4) + &
                                    & ax1*ay3*az4*mesh_from(1:nvar,i1,j3,k4) + &
                                    & ax2*ay3*az4*mesh_from(1:nvar,i2,j3,k4) + &
                                    & ax3*ay3*az4*mesh_from(1:nvar,i3,j3,k4) + &
                                    & ax4*ay3*az4*mesh_from(1:nvar,i4,j3,k4) + &
                                    & ax5*ay3*az4*mesh_from(1:nvar,i5,j3,k4) + &
                                    & ax6*ay3*az4*mesh_from(1:nvar,i6,j3,k4) + &
                                    & ax1*ay4*az4*mesh_from(1:nvar,i1,j4,k4) + &
                                    & ax2*ay4*az4*mesh_from(1:nvar,i2,j4,k4) + &
                                    & ax3*ay4*az4*mesh_from(1:nvar,i3,j4,k4) + &
                                    & ax4*ay4*az4*mesh_from(1:nvar,i4,j4,k4) + &
                                    & ax5*ay4*az4*mesh_from(1:nvar,i5,j4,k4) + &
                                    & ax6*ay4*az4*mesh_from(1:nvar,i6,j4,k4) + &
                                    & ax1*ay5*az4*mesh_from(1:nvar,i1,j5,k4) + &
                                    & ax2*ay5*az4*mesh_from(1:nvar,i2,j5,k4) + &
                                    & ax3*ay5*az4*mesh_from(1:nvar,i3,j5,k4) + &
                                    & ax4*ay5*az4*mesh_from(1:nvar,i4,j5,k4) + &
                                    & ax5*ay5*az4*mesh_from(1:nvar,i5,j5,k4) + &
                                    & ax6*ay5*az4*mesh_from(1:nvar,i6,j5,k4) + &
                                    & ax1*ay6*az4*mesh_from(1:nvar,i1,j6,k4) + &
                                    & ax2*ay6*az4*mesh_from(1:nvar,i2,j6,k4) + &
                                    & ax3*ay6*az4*mesh_from(1:nvar,i3,j6,k4) + &
                                    & ax4*ay6*az4*mesh_from(1:nvar,i4,j6,k4) + &
                                    & ax5*ay6*az4*mesh_from(1:nvar,i5,j6,k4) + &
                                    & ax6*ay6*az4*mesh_from(1:nvar,i6,j6,k4) + &
                                    & ax1*ay1*az5*mesh_from(1:nvar,i1,j1,k5) + &
                                    & ax2*ay1*az5*mesh_from(1:nvar,i2,j1,k5) + &
                                    & ax3*ay1*az5*mesh_from(1:nvar,i3,j1,k5) + &
                                    & ax4*ay1*az5*mesh_from(1:nvar,i4,j1,k5) + &
                                    & ax5*ay1*az5*mesh_from(1:nvar,i5,j1,k5) + &
                                    & ax6*ay1*az5*mesh_from(1:nvar,i6,j1,k5) + &
                                    & ax1*ay2*az5*mesh_from(1:nvar,i1,j2,k5) + &
                                    & ax2*ay2*az5*mesh_from(1:nvar,i2,j2,k5) + &
                                    & ax3*ay2*az5*mesh_from(1:nvar,i3,j2,k5) + &
                                    & ax4*ay2*az5*mesh_from(1:nvar,i4,j2,k5) + &
                                    & ax5*ay2*az5*mesh_from(1:nvar,i5,j2,k5) + &
                                    & ax6*ay2*az5*mesh_from(1:nvar,i6,j2,k5) + &
                                    & ax1*ay3*az5*mesh_from(1:nvar,i1,j3,k5) + &
                                    & ax2*ay3*az5*mesh_from(1:nvar,i2,j3,k5) + &
                                    & ax3*ay3*az5*mesh_from(1:nvar,i3,j3,k5) + &
                                    & ax4*ay3*az5*mesh_from(1:nvar,i4,j3,k5) + &
                                    & ax5*ay3*az5*mesh_from(1:nvar,i5,j3,k5) + &
                                    & ax6*ay3*az5*mesh_from(1:nvar,i6,j3,k5) + &
                                    & ax1*ay4*az5*mesh_from(1:nvar,i1,j4,k5) + &
                                    & ax2*ay4*az5*mesh_from(1:nvar,i2,j4,k5) + &
                                    & ax3*ay4*az5*mesh_from(1:nvar,i3,j4,k5) + &
                                    & ax4*ay4*az5*mesh_from(1:nvar,i4,j4,k5) + &
                                    & ax5*ay4*az5*mesh_from(1:nvar,i5,j4,k5) + &
                                    & ax6*ay4*az5*mesh_from(1:nvar,i6,j4,k5) + &
                                    & ax1*ay5*az5*mesh_from(1:nvar,i1,j5,k5) + &
                                    & ax2*ay5*az5*mesh_from(1:nvar,i2,j5,k5) + &
                                    & ax3*ay5*az5*mesh_from(1:nvar,i3,j5,k5) + &
                                    & ax4*ay5*az5*mesh_from(1:nvar,i4,j5,k5) + &
                                    & ax5*ay5*az5*mesh_from(1:nvar,i5,j5,k5) + &
                                    & ax6*ay5*az5*mesh_from(1:nvar,i6,j5,k5) + &
                                    & ax1*ay6*az5*mesh_from(1:nvar,i1,j6,k5) + &
                                    & ax2*ay6*az5*mesh_from(1:nvar,i2,j6,k5) + &
                                    & ax3*ay6*az5*mesh_from(1:nvar,i3,j6,k5) + &
                                    & ax4*ay6*az5*mesh_from(1:nvar,i4,j6,k5) + &
                                    & ax5*ay6*az5*mesh_from(1:nvar,i5,j6,k5) + &
                                    & ax6*ay6*az5*mesh_from(1:nvar,i6,j6,k5) + &
                                    & ax1*ay1*az6*mesh_from(1:nvar,i1,j1,k6) + &
                                    & ax2*ay1*az6*mesh_from(1:nvar,i2,j1,k6) + &
                                    & ax3*ay1*az6*mesh_from(1:nvar,i3,j1,k6) + &
                                    & ax4*ay1*az6*mesh_from(1:nvar,i4,j1,k6) + &
                                    & ax5*ay1*az6*mesh_from(1:nvar,i5,j1,k6) + &
                                    & ax6*ay1*az6*mesh_from(1:nvar,i6,j1,k6) + &
                                    & ax1*ay2*az6*mesh_from(1:nvar,i1,j2,k6) + &
                                    & ax2*ay2*az6*mesh_from(1:nvar,i2,j2,k6) + &
                                    & ax3*ay2*az6*mesh_from(1:nvar,i3,j2,k6) + &
                                    & ax4*ay2*az6*mesh_from(1:nvar,i4,j2,k6) + &
                                    & ax5*ay2*az6*mesh_from(1:nvar,i5,j2,k6) + &
                                    & ax6*ay2*az6*mesh_from(1:nvar,i6,j2,k6) + &
                                    & ax1*ay3*az6*mesh_from(1:nvar,i1,j3,k6) + &
                                    & ax2*ay3*az6*mesh_from(1:nvar,i2,j3,k6) + &
                                    & ax3*ay3*az6*mesh_from(1:nvar,i3,j3,k6) + &
                                    & ax4*ay3*az6*mesh_from(1:nvar,i4,j3,k6) + &
                                    & ax5*ay3*az6*mesh_from(1:nvar,i5,j3,k6) + &
                                    & ax6*ay3*az6*mesh_from(1:nvar,i6,j3,k6) + &
                                    & ax1*ay4*az6*mesh_from(1:nvar,i1,j4,k6) + &
                                    & ax2*ay4*az6*mesh_from(1:nvar,i2,j4,k6) + &
                                    & ax3*ay4*az6*mesh_from(1:nvar,i3,j4,k6) + &
                                    & ax4*ay4*az6*mesh_from(1:nvar,i4,j4,k6) + &
                                    & ax5*ay4*az6*mesh_from(1:nvar,i5,j4,k6) + &
                                    & ax6*ay4*az6*mesh_from(1:nvar,i6,j4,k6) + &
                                    & ax1*ay5*az6*mesh_from(1:nvar,i1,j5,k6) + &
                                    & ax2*ay5*az6*mesh_from(1:nvar,i2,j5,k6) + &
                                    & ax3*ay5*az6*mesh_from(1:nvar,i3,j5,k6) + &
                                    & ax4*ay5*az6*mesh_from(1:nvar,i4,j5,k6) + &
                                    & ax5*ay5*az6*mesh_from(1:nvar,i5,j5,k6) + &
                                    & ax6*ay5*az6*mesh_from(1:nvar,i6,j5,k6) + &
                                    & ax1*ay6*az6*mesh_from(1:nvar,i1,j6,k6) + &
                                    & ax2*ay6*az6*mesh_from(1:nvar,i2,j6,k6) + &
                                    & ax3*ay6*az6*mesh_from(1:nvar,i3,j6,k6) + &
                                    & ax4*ay6*az6*mesh_from(1:nvar,i4,j6,k6) + &
                                    & ax5*ay6*az6*mesh_from(1:nvar,i5,j6,k6) + &
                                    & ax6*ay6*az6*mesh_from(1:nvar,i6,j6,k6)
            END IF
          END DO ! i
        END DO ! j
      END DO ! k

    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller,'Interpolation scheme/order unknown.')
      GOTO 9999
    END IF

!---------------------------------------------------------------------------------!
! Now map the ghosts cells (overwrite) in order to get consistent values at 
! the border of the subdomains.
!---------------------------------------------------------------------------------!
    DO i = 1,pmlib_ndim
      CALL pmlib_mesh_map_ghost(topo_to,i,nvar,ierr,incl_edges = .FALSE.)
      CALL pmlib_comm_pack(mesh_to,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack(topo_to,mesh_to,0,ierr,clear=.FALSE.)
      CALL pmlib_comm_finalise(ierr)
    END DO

  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'Cell size ratio not supported.')
    GOTO 9999
  END IF


!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_interp_mesh_mesh
