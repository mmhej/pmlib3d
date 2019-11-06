!---------------------------------------------------------------------------------!
! part_deriv_vorticity.f
!---------------------------------------------------------------------------------!
! This routines computes the vorticity equation.
!---------------------------------------------------------------------------------!
SUBROUTINE part_deriv_vorticity(topo,mesh,ierr)

USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_communication

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
  INTEGER                                      :: i,j,k,l
  REAL(MK)                                     :: px,py,pz

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER                                      :: nvar

  REAL(MK)                                     :: w1, w2, w3
  REAL(MK)                                     :: ddwdx, ddwdy, ddwdz
  REAL(MK),DIMENSION(ndim,ndim)                :: S, dudx

  REAL(MK)                                     :: j1, j2, j3, r, q
  REAL(MK)                                     :: modeD

  REAL(MK)                                     :: c_1_dx2, c_1_2dx2
  REAL(MK)                                     :: c_1_dy2, c_1_2dy2
  REAL(MK)                                     :: c_1_dz2, c_1_2dz2
  REAL(MK)                                     :: c1, c2, c3, c4

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  nvar = 3

!---------------------------------------------------------------------------------!
! Calculate right-hand-side of the vorticity equation
!---------------------------------------------------------------------------------!
  ncell  = topo(rank)%ncell
  dx     = topo(rank)%dx

!---------------------------------------------------------------------------------!
! 2nd order FD
!---------------------------------------------------------------------------------!
  IF (pmlib_fd_order .EQ. 2) THEN

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients for the finite difference stencils
!---------------------------------------------------------------------------------!
    c_1_dx2 = 1.0_MK/(dx(1)**2)
    c_1_dy2 = 1.0_MK/(dx(2)**2)
    c_1_dz2 = 1.0_MK/(dx(3)**2)

    c_1_2dx2 = 1.0_MK/(2.0_MK*dx(1))
    c_1_2dy2 = 1.0_MK/(2.0_MK*dx(2))
    c_1_2dz2 = 1.0_MK/(2.0_MK*dx(3))

!---------------------------------------------------------------------------------!
! Calculate the particle derivative
!---------------------------------------------------------------------------------!
    DO k = 1,ncell(3)
      DO j = 1,ncell(2)
        DO i = 1,ncell(1)

          w1 = mesh%vort(1,i,j,k)
          w2 = mesh%vort(2,i,j,k)
          w3 = mesh%vort(3,i,j,k)

          DO l = 1,ndim
! x-derivatives
            ddwdx     = ( 1.0_MK * mesh%vort(l,i-1,j,k) &
                      & - 2.0_MK * mesh%vort(l,i  ,j,k) &
                      & + 1.0_MK * mesh%vort(l,i+1,j,k) ) * c_1_dx2
            dudx(l,1) = ( - 1.0_MK * mesh%vel(l,i-1,j,k) & 
                      & + 1.0_MK * mesh%vel(l,i+1,j,k) ) &
                      & * c_1_2dx2
! y-derivatives
            ddwdy     = ( 1.0_MK * mesh%vort(l,i,j-1,k) &
                      & - 2.0_MK * mesh%vort(l,i,j  ,k) &
                      & + 1.0_MK * mesh%vort(l,i,j+1,k) ) * c_1_dy2
            dudx(l,2) = ( - 1.0_MK * mesh%vel(l,i,j-1,k) &
                      & + 1.0_MK * mesh%vel(l,i,j+1,k) ) &
                      & * c_1_2dy2
! z-derivatives
            ddwdz     = ( 1.0_MK * mesh%vort(l,i,j,k-1) &
                      & - 2.0_MK * mesh%vort(l,i,j,k  ) &
                      & + 1.0_MK * mesh%vort(l,i,j,k+1) ) * c_1_dz2
            dudx(l,3) = ( - 1.0_MK * mesh%vel(l,i,j,k-1) &
                    & + 1.0_MK * mesh%vel(l,i,j,k+1) ) &
                    & * c_1_2dz2

! The 3D vorticity equation
            mesh%dvort(l,i,j,k) =  &
                 & visc * (ddwdx + ddwdy + ddwdz) + &
                 & w1*dudx(l,1) + w2*dudx(l,2) + w3*dudx(l,3)
          END DO

        END DO
      END DO
    END DO

!---------------------------------------------------------------------------------!
! 4nd order FD
!---------------------------------------------------------------------------------!
  ELSE IF (pmlib_fd_order .EQ. 4) THEN

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients for the finite difference stencils
!---------------------------------------------------------------------------------!
    c_1_dx2 = 1.0_MK/(dx(1)**2)
    c_1_dy2 = 1.0_MK/(dx(2)**2)
    c_1_dz2 = 1.0_MK/(dx(3)**2)

    c_1_2dx2 = 1.0_MK/(2.0_MK*dx(1))
    c_1_2dy2 = 1.0_MK/(2.0_MK*dx(2))
    c_1_2dz2 = 1.0_MK/(2.0_MK*dx(3))

    c1 = - 1.0_MK/12.0_MK
    c2 =   4.0_MK/3.0_MK

    c3 = 1.0_MK/6.0_MK
    c4 = 4.0_MK/3.0_MK

!---------------------------------------------------------------------------------!
! Calculate the particle derivative
!---------------------------------------------------------------------------------!
    DO k = 1,ncell(3)
      DO j = 1,ncell(2)
        DO i = 1,ncell(1)

          w1 = mesh%vort(1,i,j,k)
          w2 = mesh%vort(2,i,j,k)
          w3 = mesh%vort(3,i,j,k)

          DO l = 1,ndim
! x-derivative
            ddwdx     = (     c1 * mesh%vort(l,i-2,j,k) &
                      & +     c2 * mesh%vort(l,i-1,j,k) &
                      & - 2.5_MK * mesh%vort(l,i  ,j,k) &
                      & +     c2 * mesh%vort(l,i+1,j,k) &
                      & +     c1 * mesh%vort(l,i+2,j,k) &
                      & ) * c_1_dx2
            dudx(l,1) = (     c3 * mesh%vel(l,i-2,j,k) &
                      & -     c4 * mesh%vel(l,i-1,j,k) & 
                      & +     c4 * mesh%vel(l,i+1,j,k) & 
                      & -     c3 * mesh%vel(l,i+2,j,k) & 
                      & ) * c_1_2dx2
! y-derivative
            ddwdy     = (     c1 * mesh%vort(l,i,j-2,k) &
                      & +     c2 * mesh%vort(l,i,j-1,k) &
                      & - 2.5_MK * mesh%vort(l,i,j  ,k) &
                      & +     c2 * mesh%vort(l,i,j+1,k) &
                      & +     c1 * mesh%vort(l,i,j+2,k) &
                      & ) * c_1_dy2
            dudx(l,2) = (     c3 * mesh%vel(l,i,j-2,k) &
                      & -     c4 * mesh%vel(l,i,j-1,k) & 
                      & +     c4 * mesh%vel(l,i,j+1,k) & 
                      & -     c3 * mesh%vel(l,i,j+2,k) & 
                      & ) * c_1_2dy2
! z-derivative
            ddwdz     = (     c1 * mesh%vort(l,i,j,k-2) &
                      & +     c2 * mesh%vort(l,i,j,k-1) &
                      & - 2.5_MK * mesh%vort(l,i,j,k  ) &
                      & +     c2 * mesh%vort(l,i,j,k+1) &
                      & +     c1 * mesh%vort(l,i,j,k+2) &
                      & ) * c_1_dz2
            dudx(l,3) = (     c3 * mesh%vel(l,i,j,k-2) &
                    & -       c4 * mesh%vel(l,i,j,k-1) & 
                    & +       c4 * mesh%vel(l,i,j,k+1) & 
                    & -       c3 * mesh%vel(l,i,j,k+2) & 
                    & ) * c_1_2dz2

! The 3D vorticity equation
            mesh%dvort(l,i,j,k) =  &
                 & visc * (ddwdx + ddwdy + ddwdz) + &
                 & w1*dudx(l,1) + w2*dudx(l,2) + w3*dudx(l,3)
          END DO

        END DO
      END DO
    END DO 

  ENDIF

!---------------------------------------------------------------------------------!
! Map ghost cells
!---------------------------------------------------------------------------------!
  DO i = 1,ndim
    CALL pmlib_mesh_map_ghost(topo,i,nvort,ierr,incl_edges = .FALSE.)
    CALL pmlib_comm_pack(mesh%dvort,ierr)
    CALL pmlib_comm_send(ierr)
    CALL pmlib_comm_unpack(topo,mesh%dvort,0,ierr,clear=.FALSE.)
    CALL pmlib_comm_finalise(ierr)
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE part_deriv_vorticity
