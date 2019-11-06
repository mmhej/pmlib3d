!---------------------------------------------------------------------------------!
! pmlib_poisson_setup.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_poisson_setup(patch,topo_all,mesh,ierr)

USE pmlib_mod_patch
USE pmlib_mod_mesh
USE pmlib_mod_topology
USE pmlib_mod_fourier

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_patch), INTENT(IN)                          :: patch
  TYPE(class_topology_all), INTENT(INOUT)                :: topo_all
  TYPE(class_mesh), INTENT(INOUT)                        :: mesh
  INTEGER, INTENT(OUT)                                   :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                         :: caller = 'pmlib_poisson_setup'
  REAL(MK)                                  :: PI
  INTEGER                                   :: i, j, k, ii, jj, kk

  REAL(MK)                                  :: px, py, pz
  REAL(MK)                                  :: r, rho, eps

  INTEGER, DIMENSION(pmlib_ndim)            :: ncell, patch_bc
  REAL(MK), DIMENSION(pmlib_ndim)           :: xmin, xmax, dx 
  INTEGER, DIMENSION(pmlib_ndim)            :: fextend_x, fextend_y, fextend_z
  INTEGER, DIMENSION(pmlib_ndim)            :: nextend_x, nextend_y, nextend_z

  REAL(MK), DIMENSION(pmlib_ndim)           :: patch_len

  REAL(MK)                                  :: c_1_4pi, c_sq2_4pieps
  REAL(MK)                                  :: c_1_sqrt2, c_1_sqrt2pi
  REAL(MK)                                  :: c1, c2, c3, c4
  REAL(MK)                                  :: a
  REAL(MK)                                  :: c_1_a,  c_1_a2, c_1_a3
  REAL(MK)                                  :: c_1_a5, c_1_a7, c_1_a9

  REAL(MK), DIMENSION(:,:,:,:),POINTER      :: pkernel_real => NULL()

  INTEGER, DIMENSION(pmlib_ndim)            :: fft_seq  

  TYPE(class_topology),DIMENSION(:),POINTER :: topo_first => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_third => NULL()

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! De-allocate Poisson kernel
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh%pkernel_fft) )THEN 
    DEALLOCATE(mesh%pkernel_fft,stat=ierr)
  END IF
  IF(ierr .NE. 0)THEN
    CALL pmlib_write(mpi_rank,caller,'error in de-allocating the Poisson kernel.')
    GOTO 9999
  END IF


!  IF( patch%bound_cond(1) .OR. 1 .AND. &
!    & patch%bound_cond(2) .OR. 1 .AND. &
!    & patch%bound_cond(3) .OR. 1 )THEN

!---------------------------------------------------------------------------------!
! Set up temporary kernel topologies
!---------------------------------------------------------------------------------!
  DO i = 1,pmlib_ndim
    IF(patch%bound_cond(i) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(i) = 2
        nextend_x(i) = 0

        fextend_y(i) = 2
        nextend_y(i) = 0

        fextend_z(i) = 2
        nextend_z(i) = 0
!      ELSE
!        fextend_x(i) = 0
!        nextend_x(i) = pmlib_patch_ext

!        fextend_y(i) = 0
!        nextend_y(i) = pmlib_patch_ext

!        fextend_z(i) = 0
!        nextend_z(i) = pmlib_patch_ext
!      END IF
    ELSE
      fextend_x(i) = 0
      nextend_x(i) = 0

      fextend_y(i) = 0
      nextend_y(i) = 0

      fextend_z(i) = 0
      nextend_z(i) = 0
    END IF
  END DO

  CALL pmlib_topology_kernel( patch,topo_all%xpencil,1,ierr, &
                            & fextend = fextend_x, nextend = nextend_x )
  CALL pmlib_topology_kernel( patch,topo_all%ypencil,2,ierr, &
                            & fextend = fextend_y, nextend = nextend_y )
  CALL pmlib_topology_kernel( patch,topo_all%zpencil,3,ierr, &
                            & fextend = fextend_z, nextend = nextend_z )

  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to set up pencil topologies.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Determine the fft sequence
!---------------------------------------------------------------------------------!
  ALLOCATE( topo_first(0:mpi_nproc-1), topo_third(0:mpi_nproc-1) )

  fft_seq = patch%fft_seq

  IF(fft_seq(1) .EQ. 1)THEN
    DO i = 0,mpi_nproc-1
      topo_first(i) = topo_all%xpencil(i)
    END DO
  ELSEIF(fft_seq(1) .EQ. 2)THEN
    DO i = 0,mpi_nproc-1
      topo_first(i) = topo_all%ypencil(i)
    END DO
  ELSEIF(fft_seq(1) .EQ. 3)THEN
    DO i = 0,mpi_nproc-1
      topo_first(i) = topo_all%zpencil(i)
    END DO
  ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller, 'fft direction unknown.')
      GOTO 9999
  END IF

  IF(fft_seq(3) .EQ. 1)THEN
    DO i = 0,mpi_nproc-1
      topo_third(i) = topo_all%xpencil(i)
    END DO
  ELSEIF(fft_seq(3) .EQ. 2)THEN
    DO i = 0,mpi_nproc-1
      topo_third(i) = topo_all%ypencil(i)
    END DO
  ELSEIF(fft_seq(3) .EQ. 3)THEN
    DO i = 0,mpi_nproc-1
      topo_third(i) = topo_all%zpencil(i)
    END DO
  ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller, 'fft direction unknown.')
      GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Allocate kernel array
!---------------------------------------------------------------------------------!
  dx     = topo_first(mpi_rank)%dx
  ncell  = topo_first(mpi_rank)%ncell
  xmin   = topo_first(mpi_rank)%xmin

  IF(pmlib_poisson_kernel .EQ. 1)THEN
    ALLOCATE( pkernel_real(1,ncell(1),ncell(2),ncell(3)), stat=ierr)
  ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
    ALLOCATE( pkernel_real(3,ncell(1),ncell(2),ncell(3)), stat=ierr)
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'Poisson kernel unknown.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients for the kernel functions
!---------------------------------------------------------------------------------!
  eps = pmlib_regularisation_radius*MAXVAL(dx)

  c_1_4pi      = 0.25_MK/PI
  c_1_sqrt2    = 1.0_MK/SQRT(2.0_MK)
  c_1_sqrt2pi  = 1.0_MK/SQRT(2.0_MK*PI)
  c_sq2_4pieps = SQRT(2.0_MK)/(4.0_MK*PI**(1.5_MK)*eps)

  IF(patch%level .NE. 1)THEN
    a = 2.0_MK**( patch%level - patch%parent(1) ) 
    c_1_a  = 1/a
    c_1_a2 = 1/a**2
    c_1_a3 = 1/a**3
    c_1_a5 = 1/a**5
    c_1_a7 = 1/a**7
    c_1_a9 = 1/a**9
  END IF

  IF(pmlib_poisson_order .EQ. 8)THEN
    IF(pmlib_poisson_kernel .EQ. 1)THEN
      c1 = - 2.0_MK/3.0_MK
      c2 =   1.0_MK/24.0_MK
    ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
      c1 =   89.0_MK/24.0_MK
      c2 = - 20.0_MK/24.0_MK
      c3 =   1.0_MK/24.0_MK
    END IF
  END IF

  IF(pmlib_poisson_order .EQ. 10)THEN
    IF(pmlib_poisson_kernel .EQ. 1)THEN
      c1 = - 233.0_MK/192.0_MK
      c2 =   29.0_MK/192.0_MK
      c3 = - 1.0_MK/192.0_MK
    ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN
      c1 = 1027.0_MK/192.0_MK
      c2 = - 349.0_MK/192.0_MK 
      c3 =   35.0_MK/192.0_MK
      c4 = - 1.0_MK/192.0_MK
    END IF
  END IF

!---------------------------------------------------------------------------------!
! Compute patch length
!---------------------------------------------------------------------------------!
  DO i = 1,pmlib_ndim
    IF( patch%bound_cond(i) .EQ. 0 )THEN
      patch_len(i) = dx(i)*REAL(2*(patch%ncell(i) + 2*patch%nghost), MK)
    ELSE
      patch_len(i) = dx(i)*REAL(patch%ncell(i), MK)
    END IF
  END DO

!---------------------------------------------------------------------------------!
! Compute Poisson kernel
! Do summation of images in case of periodic boundary conditions
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! 0th order kernel
!---------------------------------------------------------------------------------!
  IF (pmlib_poisson_order .EQ. 0) THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 1.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + c_1_4pi/r
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) &
                       & - c_1_4pi * px/r**3
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) &
                       & - c_1_4pi * py/r**3
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) &
                       & - c_1_4pi * pz/r**3
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller, &
        & 'Regularisation must be used for multi-patch simulations.')
      GOTO 9999
    END IF

!---------------------------------------------------------------------------------!
! 2nd order kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 2)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_sq2_4pieps
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_1_4pi/r * ERF(c_1_sqrt2 * rho)
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & - c_1_4pi * px/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & - c_1_4pi * py/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & - c_1_4pi * pz/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_sq2_4pieps * ( 1.0_MK - c_1_a)
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & + c_1_4pi/r *                                            &
                       & ( ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & - c_1_4pi * px/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho * c_1_a                               &
                       &   ) * EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi       &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)                     &
                       & - c_1_4pi * py/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho * c_1_a                               &
                       &   ) * EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi       &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)                     &
                       & - c_1_4pi * pz/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho * c_1_a                               &
                       &   ) * EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi       &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      END IF
    END IF

!---------------------------------------------------------------------------------!
! 4th order kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 4)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 1.5_MK * c_sq2_4pieps
              ELSE
                pkernel_real(1,i,j,k) =  pkernel_real(1,i,j,k)           &
                       & + c_1_4pi/r *                                   &
                       & ( rho * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi      &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & - c_1_4pi * px/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +          rho**3                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & - c_1_4pi * py/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +          rho**3                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & - c_1_4pi * pz/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +          rho**3                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 1.5_MK * c_sq2_4pieps * ( 1.0_MK - c_1_a)
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & + c_1_4pi/r *                                            &
                       & ( rho * EXP(-0.5_MK*rho**2)*c_1_sqrt2pi                  &
                       & - rho * EXP(-0.5_MK*rho**2*c_1_a2)*c_1_sqrt2pi * c_1_a   &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & - c_1_4pi * px/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &     +          rho**3                                    &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho    * c_1_a                            &
                       &     +          rho**3 * c_1_a3                           &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)                     &
                       & - c_1_4pi * py/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &     +          rho**3                                    &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho    * c_1_a                            &
                       &     +          rho**3 * c_1_a3                           &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)                     &
                       & - c_1_4pi * pz/r**3 *                                    &
                       & ( ( - 2.0_MK * rho                                       &
                       &     +          rho**3                                    &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.0_MK * rho    * c_1_a                            &
                       &     +          rho**3 * c_1_a3                           &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      END IF
    END IF
!---------------------------------------------------------------------------------!
! 6th order kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 6)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 1.875_MK * c_sq2_4pieps
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_1_4pi/r *                                   &
                       & ((   1.75_MK * rho                              &
                       &    - 0.25_MK * rho**3                           &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & - c_1_4pi * px/r**3 *                           &
                       & (( - 2.00_MK * rho                              &
                       &    + 2.25_MK * rho**3                           &
                       &    - 0.25_MK * rho**5                           &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & - c_1_4pi * py/r**3 *                           &
                       & (( - 2.00_MK * rho                              &
                       &    + 2.25_MK * rho**3                           &
                       &    - 0.25_MK * rho**5                           &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & - c_1_4pi * pz/r**3 *                           &
                       & (( - 2.00_MK * rho                              &
                       &    + 2.25_MK * rho**3                           &
                       &    - 0.25_MK * rho**5                           &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 1.875_MK * c_sq2_4pieps * ( 1.0_MK - c_1_a)
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & + c_1_4pi/r *                                            &
                       & ( (   1.75_MK * rho                                      &
                       &     - 0.25_MK * rho**3                                   &
                       &   ) * EXP(-0.5_MK*rho**2)*c_1_sqrt2pi                    &
                       & - (   1.75_MK * rho    * c_1_a                           &
                       &     - 0.25_MK * rho**3 * c_1_a3                          &
                       &   ) * EXP(-0.5_MK*rho**2*c_1_a2)*c_1_sqrt2pi             &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & - c_1_4pi * px/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     + 2.25_MK * rho**3                                   &
                       &     - 0.25_MK * rho**5                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     + 2.25_MK * rho**3 * c_1_a3                          &
                       &     - 0.25_MK * rho**5 * c_1_a5                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)                     &
                       & - c_1_4pi * py/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     + 2.25_MK * rho**3                                   &
                       &     - 0.25_MK * rho**5                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     + 2.25_MK * rho**3 * c_1_a3                          &
                       &     - 0.25_MK * rho**5 * c_1_a5                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)                     &
                       & - c_1_4pi * pz/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     + 2.25_MK * rho**3                                   &
                       &     - 0.25_MK * rho**5                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     + 2.25_MK * rho**3 * c_1_a3                          &
                       &     - 0.25_MK * rho**5 * c_1_a5                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      END IF
    END IF

!---------------------------------------------------------------------------------!
! 8th order kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 8)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 2.1875_MK * c_sq2_4pieps
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_1_4pi/r *                                   &
                       & ((   2.375_MK * rho                             &
                       &    +       c1 * rho**3                          &
                       &    +       c2 * rho**5                          &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & - c_1_4pi * px/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +     c1 * rho**3                            &
                       &    +     c2 * rho**5                            &
                       &    +     c3 * rho**7                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & - c_1_4pi * py/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +     c1 * rho**3                            &
                       &    +     c2 * rho**5                            &
                       &    +     c3 * rho**7                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & - c_1_4pi * pz/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +     c1 * rho**3                            &
                       &    +     c2 * rho**5                            &
                       &    +     c3 * rho**7                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 2.1875_MK * c_sq2_4pieps * ( 1.0_MK - c_1_a)
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & + c_1_4pi/r *                                            &
                       & ( (   2.375_MK * rho                                     &
                       &     +       c1 * rho**3                                  &
                       &     +       c2 * rho**5                                  &
                       &   ) * EXP(-0.5_MK*rho**2)*c_1_sqrt2pi                    &
                       & - (   2.375_MK * rho * c_1_a                             &
                       &     +       c1 * rho**3 * c_1_a3                         &
                       &     +       c2 * rho**5 * c_1_a5                         &
                       &   ) * EXP(-0.5_MK*rho**2*c_1_a2)*c_1_sqrt2pi             &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & - c_1_4pi * px/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)                     &
                       & - c_1_4pi * py/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)                     &
                       & - c_1_4pi * pz/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      END IF
    END IF

!---------------------------------------------------------------------------------!
! 10th order kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 10)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 2.4609375_MK * c_sq2_4pieps
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + c_1_4pi/r *                                   &
                       & ((   2.921875_MK * rho                          &
                       &    +          c1 * rho**3                       &
                       &    +          c2 * rho**5                       &
                       &    +          c3 * rho**7                       &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & - c_1_4pi * px/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +     c1 * rho**3                            &
                       &    +     c2 * rho**5                            &
                       &    +     c3 * rho**7                            &
                       &    +     c4 * rho**9                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & - c_1_4pi * py/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    + c1 * rho**3                                &
                       &    + c2 * rho**5                                &
                       &    + c3 * rho**7                                &
                       &    + c4 * rho**9                                &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & - c_1_4pi * pz/r**3 *                           &
                       & (( - 2.0_MK * rho                               &
                       &    +     c1 * rho**3                            &
                       &    +     c2 * rho**5                            &
                       &    +     c3 * rho**7                            &
                       &    +     c4 * rho**9                            &
                       &  ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi         &
                       & + ERF(c_1_sqrt2 * rho) )
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 2.4609375_MK * c_sq2_4pieps * ( 1.0_MK - c_1_a)
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & + c_1_4pi/r *                                            &
                       & ( (   2.921875_MK * rho                                  &
                       &     +       c1 * rho**3                                  &
                       &     +       c2 * rho**5                                  &
                       &     +       c3 * rho**7                                  &
                       &   ) * EXP(-0.5_MK*rho**2)*c_1_sqrt2pi                    &
                       & - (   2.921875_MK * rho * c_1_a                          &
                       &     +       c1 * rho**3 * c_1_a3                         &
                       &     +       c2 * rho**5 * c_1_a5                         &
                       &     +       c3 * rho**7 * c_1_a7                         &
                       &   ) * EXP(-0.5_MK*rho**2*c_1_a2)*c_1_sqrt2pi             &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK
        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = r/eps

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)                     &
                       & - c_1_4pi * px/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &     +      c4 * rho**9                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &     +      c4 * rho**9 * c_1_a9                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)                     &
                       & - c_1_4pi * py/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &     +      c4 * rho**9                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &     +      c4 * rho**9 * c_1_a9                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)                     &
                       & - c_1_4pi * pz/r**3 *                                    &
                       & ( ( - 2.00_MK * rho                                      &
                       &     +      c1 * rho**3                                   &
                       &     +      c2 * rho**5                                   &
                       &     +      c3 * rho**7                                   &
                       &     +      c4 * rho**9                                   &
                       &   ) * EXP(-0.5_MK *rho**2) * c_1_sqrt2pi                 &
                       & - ( - 2.00_MK * rho    * c_1_a                           &
                       &     +      c1 * rho**3 * c_1_a3                          &
                       &     +      c2 * rho**5 * c_1_a5                          &
                       &     +      c3 * rho**7 * c_1_a7                          &
                       &     +      c4 * rho**9 * c_1_a9                          &
                       &   )* EXP(-0.5_MK *rho**2 * c_1_a2 ) * c_1_sqrt2pi        &
                       & + ERF(c_1_sqrt2 * rho) - ERF(c_1_sqrt2 * rho * c_1_a) )
              END IF

            END DO !i
          END DO !j
        END DO !k

      END IF
    END IF

!---------------------------------------------------------------------------------!
! Bessel kernel
!---------------------------------------------------------------------------------!
  ELSEIF(pmlib_poisson_order .EQ. 100)THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Base patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    IF(patch%level .EQ. 1)THEN
      IF(pmlib_poisson_kernel .EQ. 1)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = PI*r/MAXVAL(dx)

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + 1.0_MK/(2.0_MK * PI * MAXVAL(dx))
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + SI(rho)/(2.0_MK * PI**2 * r)
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      ELSEIF(pmlib_poisson_kernel .EQ. 2)THEN

        pkernel_real = 0.0_MK 
        DO kk = -patch%bound_cond(3),patch%bound_cond(3)
          DO jj = -patch%bound_cond(2),patch%bound_cond(2)
            DO ii = -patch%bound_cond(1),patch%bound_cond(1)

        DO k = 1,ncell(3)
          pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3) &
             & + patch_len(3) * REAL(kk, MK)

          DO j = 1,ncell(2)
            py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2) &
               & + patch_len(2) * REAL(jj, MK)

            DO i = 1,ncell(1)
              px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1) &
                 & + patch_len(1) * REAL(ii, MK)

              r   = SQRT(px*px + py*py + pz*pz)      
              rho = PI*r/MAXVAL(dx)

              IF(r .EQ. 0.0_MK)THEN
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k) + 0.0_MK
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k) + 0.0_MK
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k) + 0.0_MK
              ELSE
                pkernel_real(1,i,j,k) = pkernel_real(1,i,j,k)            &
                       & + px * ( SIN(rho)-SI(rho) )/(2.0_MK * PI**2 * r**3)
                pkernel_real(2,i,j,k) = pkernel_real(2,i,j,k)            &
                       & + py * ( SIN(rho)-SI(rho) )/(2.0_MK * PI**2 * r**3)
                pkernel_real(3,i,j,k) = pkernel_real(3,i,j,k)            &
                       & + pz * ( SIN(rho)-SI(rho) )/(2.0_MK * PI**2 * r**3)
              END IF

            END DO !i
          END DO !j
        END DO !k

            END DO !ii
          END DO !jj
        END DO !kk

      END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Refinement patch kernel
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ELSE
      ierr = -1
      CALL pmlib_write(mpi_rank,caller, &
        & 'Gaussian regularisation must be used for multi-patch simulations.')
      GOTO 9999
    END IF
!---------------------------------------------------------------------------------!
! Unknown kernel
!---------------------------------------------------------------------------------!
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller, &
         & 'Integration kernel unknown - change POISSON ORDER.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Fourier transform field
!---------------------------------------------------------------------------------!
  CALL pmlib_fft( topo_all,topo_first,pkernel_real, &
                & topo_third,mesh%pkernel_fft, &
                & fft_seq,ierr,fft_shift = .TRUE.)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to Fourier transform.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Deallocate real integration kernel
!---------------------------------------------------------------------------------!
  DEALLOCATE(pkernel_real,stat=ierr)
  DEALLOCATE(topo_first,stat=ierr)
  DEALLOCATE(topo_third,stat=ierr)

!---------------------------------------------------------------------------------!
! Initiate partially extended pencil topologies for the Poisson solver 
! cf. Hockney:1988
!---------------------------------------------------------------------------------!
  DEALLOCATE(topo_all%xpencil, topo_all%ypencil, topo_all%zpencil,stat=ierr) 
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller, 'Failed to deallocate topologies.')
    GO TO 9999
  END IF

!  END IF ! patch%bound_cond

  fextend_x = 0
  nextend_x = 0
  fextend_y = 0
  nextend_y = 0
  fextend_z = 0
  nextend_z = 0

  IF(fft_seq(2) .EQ. 1)THEN
    IF(patch%bound_cond(1) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(1) = 2
        nextend_x(1) = 0

        fextend_y(1) = 2
        nextend_y(1) = 0

        fextend_z(1) = 2
        nextend_z(1) = 0
!      ELSE
!        fextend_x(1) = 0
!        nextend_x(1) = pmlib_patch_ext

!        fextend_y(1) = 0
!        nextend_y(1) = pmlib_patch_ext

!        fextend_z(1) = 0
!        nextend_z(1) = pmlib_patch_ext
!      END IF
    END IF
  ELSEIF(fft_seq(2) .EQ. 2)THEN
    IF(patch%bound_cond(2) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(2) = 2
        nextend_x(2) = 0

        fextend_y(2) = 2
        nextend_y(2) = 0

        fextend_z(2) = 2
        nextend_z(2) = 0
!      ELSE
!        fextend_x(2) = 0
!        nextend_x(2) = pmlib_patch_ext

!        fextend_y(2) = 0
!        nextend_y(2) = pmlib_patch_ext

!        fextend_z(2) = 0
!        nextend_z(2) = pmlib_patch_ext
!      END IF
    END IF
  ELSEIF(fft_seq(2) .EQ. 3)THEN
    IF(patch%bound_cond(3) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(3) = 2
        nextend_x(3) = 0

        fextend_y(3) = 2
        nextend_y(3) = 0

        fextend_z(3) = 2
        nextend_z(3) = 0
!      ELSE
!        fextend_x(3) = 0
!        nextend_x(3) = pmlib_patch_ext

!        fextend_y(3) = 0
!        nextend_y(3) = pmlib_patch_ext

!        fextend_z(3) = 0
!        nextend_z(3) = pmlib_patch_ext
!      END IF
    END IF
  END IF

  IF(fft_seq(1) .EQ. 1)THEN
    IF(patch%bound_cond(1) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(1) = 2
        nextend_x(1) = 0

        fextend_y(1) = 2
        nextend_y(1) = 0

        fextend_z(1) = 2
        nextend_z(1) = 0
!      ELSE
!        fextend_x(1) = 0
!        nextend_x(1) = pmlib_patch_ext

!        fextend_y(1) = 0
!        nextend_y(1) = pmlib_patch_ext

!        fextend_z(1) = 0
!        nextend_z(1) = pmlib_patch_ext
!      END IF
      fextend_x(2) = 0
      nextend_x(2) = 0

      fextend_x(3) = 0
      nextend_x(3) = 0
    END IF
  ELSEIF(fft_seq(1) .EQ. 2)THEN
    IF(patch%bound_cond(2) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(2) = 2
        nextend_x(2) = 0

        fextend_y(2) = 2
        nextend_y(2) = 0

        fextend_z(2) = 2
        nextend_z(2) = 0
!      ELSE
!        fextend_x(2) = 0
!        nextend_x(2) = pmlib_patch_ext

!        fextend_y(2) = 0
!        nextend_y(2) = pmlib_patch_ext

!        fextend_z(2) = 0
!        nextend_z(2) = pmlib_patch_ext
!      END IF
      fextend_y(1) = 0
      nextend_y(1) = 0

      fextend_y(3) = 0
      nextend_y(3) = 0
    END IF
  ELSEIF(fft_seq(1) .EQ. 3)THEN
    IF(patch%bound_cond(3) .EQ. 0)THEN
!      IF(patch%level .EQ. 1)THEN
        fextend_x(3) = 2
        nextend_x(3) = 0

        fextend_y(3) = 2
        nextend_y(3) = 0

        fextend_z(3) = 2
        nextend_z(3) = 0
!      ELSE
!        fextend_x(3) = 0
!        nextend_x(3) = pmlib_patch_ext

!        fextend_y(3) = 0
!        nextend_y(3) = pmlib_patch_ext

!        fextend_z(3) = 0
!        nextend_z(3) = pmlib_patch_ext
!      END IF
      fextend_z(1) = 0
      nextend_z(1) = 0

      fextend_z(2) = 0
      nextend_z(2) = 0
    END IF
  END IF

  CALL pmlib_topology_pencil(patch,topo_all%xpencil,1,ierr, &
                            & fextend = fextend_x, nextend = nextend_x)
  CALL pmlib_topology_pencil(patch,topo_all%ypencil,2,ierr, &
                            & fextend = fextend_y, nextend = nextend_y)
  CALL pmlib_topology_pencil(patch,topo_all%zpencil,3,ierr, &
                            & fextend = fextend_z, nextend = nextend_z)

  IF (ierr.NE.0) THEN
    CALL pmlib_write(mpi_rank,caller, 'Failed to create pencil topologies.')
    GO TO 9999
  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_poisson_setup



REAL(MK) FUNCTION SI(x)

  INTEGER, PARAMETER                 :: MK = KIND(1D0)
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  REAL(MK), INTENT(IN)               :: x

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                           :: PI
  REAL(MK),PARAMETER                 :: eul = 0.5772156649015328_MK
  INTEGER,PARAMETER                  :: maxit = 1000
  REAL(MK),PARAMETER                 :: eps = 1.0E-10_MK
  REAL(MK),PARAMETER                 :: fpmin = 1.0E-30_MK
  REAL(MK), PARAMETER                :: fpmax = 1.0E30_MK
  REAL(MK), PARAMETER                :: tmin = 2.0

  INTEGER     :: i,k

  REAL(MK)    :: a,err,fact,sign,sum,sumc,sums,t,term
  COMPLEX(MK) :: h,b,c,d,del
  LOGICAL     :: odd


  PI = ACOS(-1.0_MK)

  t = ABS(x)
  IF(t == 0.0)THEN
    SI = 0.0
    RETURN
  END IF

! Evaluate continued fraction by modified Lentz’s method (§5.2).
  IF(t > tmin)THEN
    b = CMPLX(1.0_MK,t, KIND=MK)
    c = fpmax
    d = 1.0_MK/b
    h = d
    DO i = 2,maxit
      a = -(i-1)**2
      b = b+2.0_MK
      d = 1.0_MK/(a*d+b)  ! Denominators cannot be zero.
      c = b+a/c
      del = c*d
      h = h*del

      IF( ABS( REAL(del-1.0_MK,MK) ) + ABS( AIMAG(del-1.0_MK) ) <= eps ) EXIT

    END DO

    IF(i > maxit) WRITE(*,*)'continued fraction failed in cisi'

      h = CMPLX( COS(t), -SIN(t) ,KIND=MK) * h

      SI = 0.5_MK*PI + AIMAG(h)
    ELSE
! Evaluate both series simultaneously.
    IF(t < SQRT(fpmin))THEN
! Special case: avoid failure of convergence test because of underflow.
      sumc = 0.0_MK
      sums = t
    ELSE
      sum  = 0.0_MK
      sums = 0.0_MK
      sumc = 0.0_MK
      sign = 1.0_MK
      fact = 1.0_MK
      odd  = .TRUE.
      DO k = 1,maxit
        fact = fact*t/k
        term = fact/k
        sum = sum + sign*term
        err = term/abs(sum)
        IF(odd)THEN
          sign = -sign
          sums = sum
          sum  = sumc
        ELSE
          sumc = sum
          sum  = sums
        END IF

        IF (err < EPS) EXIT

        odd = .NOT. odd

      END DO

      IF(k > maxit) WRITE(*,*)'MAXIT exceeded in cisi'
    END IF

    SI = sums
  END IF

  IF(x < 0.0_MK) SI=-SI

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
RETURN
END FUNCTION SI
