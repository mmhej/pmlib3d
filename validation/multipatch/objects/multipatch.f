# 1 "sources/multipatch.f"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "sources/multipatch.f"
!---------------------------------------------------------------------------------!
! multipatch.f
! version: 3D multi-core
!---------------------------------------------------------------------------------!
! <<HEADER>>
!---------------------------------------------------------------------------------!
PROGRAM multipatch

USE mod_output

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_poisson
USE pmlib_mod_interpolation
USE pmlib_mod_output
USE pmlib_mod_communication

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER :: ndim = pmlib_ndim
  INTEGER, PARAMETER :: nvort = pmlib_nvort

  CHARACTER(LEN=7), PARAMETER :: caller = 'poisson'
  REAL(MK),PARAMETER :: PI = ACOS(-1.0_MK)
  INTEGER :: ierr

  INTEGER :: ilevel,ipatch
  INTEGER :: i,j,k
  REAL(MK) :: px,py,pz

  REAL(MK) :: r,r0,c

  CHARACTER(LEN=256) :: nodename
  INTEGER :: name_len

  REAL(MK),DIMENSION(ndim) :: xmin, dx
  INTEGER,DIMENSION(ndim) :: ncell
  INTEGER,DIMENSION(2*ndim) :: nghost
  INTEGER,DIMENSION(2) :: parent

  CHARACTER(LEN=256) :: pltfile

  REAL(MK) :: max_vel, max_vort, max_err
  REAL(MK) :: error_vel, error_vort
  REAL(MK) :: diff_vel, diff_vort
  REAL(MK),DIMENSION(ndim) :: vel, vel_inf, vort

  INTEGER :: n, nn
  REAL(MK),DIMENSION(:),POINTER :: global_error_vel
  REAL(MK),DIMENSION(:),POINTER :: global_error_vort
  REAL(MK),DIMENSION(:),POINTER :: dxs
  INTEGER,DIMENSION(:),POINTER :: nx

  REAL(MK) :: max_error_vel
  REAL(MK) :: global_int_vel, global_int_vort
  REAL(MK) :: int_vel, int_vort

  REAL(MK) :: rho, phi, theta
  REAL(MK) :: vort_mag, vel_mag
  REAL(MK) :: sigma, circ

  LOGICAL :: regularise_vort

  INTEGER :: comm
  INTEGER :: rank
  INTEGER :: nproc

  REAL(MK),DIMENSION(:,:,:,:),POINTER :: pltfld => NULL()

!---------------------------------------------------------------------------------!
! Construct classes
!---------------------------------------------------------------------------------!
  TYPE(class_patch),DIMENSION(:,:),POINTER :: patch => NULL()
  TYPE(class_topology_all),DIMENSION(:,:),POINTER :: topo_all => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE :: mesh => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE :: mean => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE :: pmesh => NULL()
  TYPE(class_mesh),DIMENSION(:,:),POINTER,SAVE :: diff => NULL()

!---------------------------------------------------------------------------------!
! Initiate program
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Start MPI
!---------------------------------------------------------------------------------!
  comm = MPI_COMM_WORLD
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(comm,rank,ierr)
  CALL MPI_COMM_SIZE(comm,nproc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(nodename, name_len, ierr)
  WRITE(*,*) 'rank: ',rank, ' starting on "',nodename(1:name_len),'"'

  mpi_comm = comm
  mpi_rank = rank
  mpi_nproc = nproc

!---------------------------------------------------------------------------------!
! Set up patch structure
!---------------------------------------------------------------------------------!
  pmlib_nlevel = 2
  pmlib_max_npatch = 1

  ALLOCATE( pmlib_npatch(pmlib_nlevel) )
! pmlib_npatch = (/ 1 /)
  pmlib_npatch = (/ 1, 1 /)
! pmlib_npatch = (/ 1, 4 /)

! pmlib_npatch = (/ 1, 0, 1 /)
! pmlib_npatch = (/ 1, 0, 4 /)

!---------------------------------------------------------------------------------!
! Set up patch and topology
!---------------------------------------------------------------------------------!
  ALLOCATE( patch(pmlib_nlevel, pmlib_max_npatch ) , stat=ierr )
  ALLOCATE( topo_all(pmlib_nlevel, pmlib_max_npatch ), stat=ierr )
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to allocate patch structure.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Set up parameters
!---------------------------------------------------------------------------------!
  vel_inf = (/ 0.0_MK, 0.0_MK, 0.0_MK /)

  patch(1,1)%level = 1
  patch(1,1)%parent = (/ 0, 0 /)
  patch(1,1)%ptype = 1
  patch(1,1)%xmin = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
  patch(1,1)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
  patch(1,1)%bound_cond = (/ 0, 0, 0 /)
  patch(1,1)%nghost = 4

  IF(pmlib_nlevel .GT. 1)THEN
    IF( pmlib_npatch(2) .EQ. 1 )THEN
      patch(2,1)%level = 2
      patch(2,1)%parent = (/ 1, 1 /)
      patch(2,1)%ptype = 1
      patch(2,1)%xmin = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
      patch(2,1)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
      patch(2,1)%bound_cond = (/ 0, 0, 0 /)
      patch(2,1)%nghost = 8
    ELSE IF( pmlib_npatch(2) .EQ. 4 )THEN
      patch(2,1)%level = 2
      patch(2,1)%parent = (/ 1, 1 /)
      patch(2,1)%ptype = 1
      patch(2,1)%xmin = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
      patch(2,1)%xmax = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
      patch(2,1)%bound_cond = (/ 0, 0, 0 /)
      patch(2,1)%nghost = 32

      patch(2,2)%level = 2
      patch(2,2)%parent = (/ 1, 1 /)
      patch(2,2)%ptype = 1
      patch(2,2)%xmin = (/ -1.0_MK, 0.0_MK, -1.0_MK /)
      patch(2,2)%xmax = (/ 1.0_MK, 1.0_MK, 0.0_MK /)
      patch(2,2)%bound_cond = (/ 0, 0, 0 /)
      patch(2,2)%nghost = 32

      patch(2,3)%level = 2
      patch(2,3)%parent = (/ 1, 1 /)
      patch(2,3)%ptype = 1
      patch(2,3)%xmin = (/ -1.0_MK, -1.0_MK, 0.0_MK /)
      patch(2,3)%xmax = (/ 1.0_MK, 0.0_MK, 1.0_MK /)
      patch(2,3)%bound_cond = (/ 0, 0, 0 /)
      patch(2,3)%nghost = 32

      patch(2,4)%level = 2
      patch(2,4)%parent = (/ 1, 1 /)
      patch(2,4)%ptype = 1
      patch(2,4)%xmin = (/ -1.0_MK, 0.0_MK, 0.0_MK /)
      patch(2,4)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
      patch(2,4)%bound_cond = (/ 0, 0, 0 /)
      patch(2,4)%nghost = 32
    END IF
  END IF
  IF(pmlib_nlevel .GT. 2)THEN
    IF( pmlib_npatch(3) .EQ. 1 )THEN
      patch(3,1)%level = 3
      patch(3,1)%parent = (/ 1, 1 /)
      patch(3,1)%ptype = 1
      patch(3,1)%xmin = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
      patch(3,1)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
      patch(3,1)%bound_cond = (/ 0, 0, 0 /)
      patch(3,1)%nghost = 8
    ELSE IF( pmlib_npatch(3) .EQ. 4 )THEN
      patch(3,1)%level = 3
      patch(3,1)%parent = (/ 1, 1 /)
      patch(3,1)%ptype = 1
      patch(3,1)%xmin = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
      patch(3,1)%xmax = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
      patch(3,1)%bound_cond = (/ 0, 0, 0 /)
      patch(3,1)%nghost = 32

      patch(3,2)%level = 3
      patch(3,2)%parent = (/ 1, 1 /)
      patch(3,2)%ptype = 1
      patch(3,2)%xmin = (/ -1.0_MK, 0.0_MK, -1.0_MK /)
      patch(3,2)%xmax = (/ 1.0_MK, 1.0_MK, 0.0_MK /)
      patch(3,2)%bound_cond = (/ 0, 0, 0 /)
      patch(3,2)%nghost = 32

      patch(3,3)%level = 3
      patch(3,3)%parent = (/ 1, 1 /)
      patch(3,3)%ptype = 1
      patch(3,3)%xmin = (/ -1.0_MK, -1.0_MK, 0.0_MK /)
      patch(3,3)%xmax = (/ 1.0_MK, 0.0_MK, 1.0_MK /)
      patch(3,3)%bound_cond = (/ 0, 0, 0 /)
      patch(3,3)%nghost = 32

      patch(3,4)%level = 3
      patch(3,4)%parent = (/ 1, 1 /)
      patch(3,4)%ptype = 1
      patch(3,4)%xmin = (/ -1.0_MK, 0.0_MK, 0.0_MK /)
      patch(3,4)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
      patch(3,4)%bound_cond = (/ 0, 0, 0 /)
      patch(3,4)%nghost = 32
    END IF
  END IF

  regularise_vort = .TRUE.

  pmlib_regularisation_radius = 2.0_MK
  pmlib_poisson_order = 2
  pmlib_poisson_kernel = 1
  pmlib_time_integration_order = 3
  pmlib_fd_order = 4
  pmlib_interpolation_order = 4

!---------------------------------------------------------------------------------!
! Check that pmlib parameters have been set
!---------------------------------------------------------------------------------!
  CALL pmlib_parameters_check(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'pmlib parameters are not set correctly.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Setup convergence parameters
!---------------------------------------------------------------------------------!
  nn = 2
  ALLOCATE( nx(nn), dxs(nn), global_error_vel(nn), global_error_vort(nn) )
! nx = (/ 32, 64, 128 /)
  nx = (/ 16, 32 /)
! nx = 100

!---------------------------------------------------------------------------------!
! Loop through different resolutions
!---------------------------------------------------------------------------------!
  DO n = 1,nn
    patch(1,1)%dx = 1.0_MK/REAL(nx(n),MK)
    IF(rank .EQ. 0) WRITE(*,*) 'Resolution: ', nx(n)

!---------------------------------------------------------------------------------!
! Adjust patches
!---------------------------------------------------------------------------------!
    CALL pmlib_patch_adjust(patch,ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to adjust patches.')
      GOTO 9999
    ENDIF

! enforce fft sequence
! patch%fft_seq = (/ 1 , 2 , 3 /)
! patch%fft_seq = (/ 3 , 1 , 2 /)
! patch%fft_seq = (/ 2 , 3 , 1 /)

! patch%fft_seq = (/ 2 , 1 , 3 /)
! patch%fft_seq = (/ 3 , 2 , 1 /)
! patch%fft_seq = (/ 1 , 3 , 2 /)

!---------------------------------------------------------------------------------!
! Create main cuboid topology
!---------------------------------------------------------------------------------!
   DO ilevel = 1,pmlib_nlevel
    DO ipatch = 1,pmlib_npatch(ilevel)

      CALL pmlib_topology_cuboid( patch(ilevel,ipatch), &
                                & topo_all(ilevel,ipatch)%cuboid,ierr )

      IF(ilevel .GT. 1)THEN
        CALL pmlib_topology_parent( topo_all(ilevel,ipatch)%cuboid, &
                                  & topo_all(ilevel,ipatch)%parent,ierr )
      END IF

    END DO
  END DO
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to create topology.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Allocate mesh set
!---------------------------------------------------------------------------------!
  CALL pmlib_mesh_allocate( topo_all,mesh,ierr, &
                          & vort = .TRUE., dvort = .TRUE., &
                          & vel = .TRUE., mask = .TRUE.)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to allocate mesh.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Setup Greens function
!---------------------------------------------------------------------------------!
  DO ilevel = 1,pmlib_nlevel
    DO ipatch = 1,pmlib_npatch(ilevel)
      CALL pmlib_poisson_setup( patch(ilevel,ipatch),topo_all(ilevel,ipatch), &
                              & mesh(ilevel,ipatch),ierr)
    END DO !ipatch
  END DO !ilevel
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Setup mask array
!---------------------------------------------------------------------------------!
  DO ilevel = 1,pmlib_nlevel
    DO ipatch = 1,pmlib_npatch(ilevel)

      mesh(ilevel,ipatch)%patch_mask = .TRUE.

    END DO !ipatch
  END DO !ilevel

  DO ilevel = pmlib_nlevel,2,-1
    DO ipatch = 1,pmlib_npatch(ilevel)

      parent = topo_all(ilevel,ipatch)%cuboid(rank)%parent

      pmlib_mesh_mean(ilevel,ipatch)%patch_mask = .FALSE.

! Map mean to parent vorticity (overwrite)
      CALL pmlib_mesh_map( topo_all(ilevel,ipatch)%parent, &
                         & topo_all(parent(1),parent(2))%cuboid,pmlib_nvort,ierr, &
                         & map_ghost=.FALSE.)
      CALL pmlib_comm_pack(pmlib_mesh_mean(ilevel,ipatch)%patch_mask,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack( topo_all(parent(1),parent(2))%cuboid, &
                            & mesh(parent(1),parent(2))%patch_mask, &
                            & 0,ierr,clear=.FALSE.)
      CALL pmlib_comm_finalise(ierr)

      IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller, &
                        & 'Failed to decompose vorticity on patches.')
        GOTO 9999
      END IF

    END DO !ipatch
  END DO !ilevel

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  c = 10.0_MK
  r0 = 0.5_MK

  circ = 7.5_MK
  sigma = 0.2_MK*r0

  DO ilevel = 1,pmlib_nlevel
    DO ipatch = 1,pmlib_npatch(ilevel)
      dx = topo_all(ilevel,ipatch)%cuboid(rank)%dx
      xmin = topo_all(ilevel,ipatch)%cuboid(rank)%xmin
      ncell = topo_all(ilevel,ipatch)%cuboid(rank)%ncell
      nghost = topo_all(ilevel,ipatch)%cuboid(rank)%nghost

      mesh(ilevel,ipatch)%vort = 0.0_MK
      mesh(ilevel,ipatch)%vel = 0.0_MK

      DO k = 1-nghost(5),ncell(3)+nghost(6)
        pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)

        DO j = 1-nghost(3),ncell(2)+nghost(4)
          py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)

          DO i = 1-nghost(1),ncell(1)+nghost(2)
            px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

            rho = SQRT(pz*pz + py*py)
            phi = SQRT((rho - r0)**2 + px**2)
            theta = ATAN2(pz,py)

            IF( phi .LT. r0) THEN

              vort_mag = -EXP(- c*r0**2/(2.0_MK*r0*rho - rho**2 - px**2)) * &
                       & ( 4.0_MK*c**2*r0**4*px**2*rho**2 &
                       & - 16.0_MK*r0**4*rho**4 &
                       & + 32.0_MK*r0**3*rho**5 &
                       & - 24.0_MK*r0**2*rho**6 &
                       & + 8.0_MK*r0*rho**7 &
                       & - 4.0_MK*rho**6*px**2 &
                       & - 6.0_MK*rho**4*px**4 &
                       & - 4.0_MK*rho**2*px**6 &
                       & - 8.0_MK*c*r0**5*rho**3 &
                       & + 8.0_MK*c*r0**4*rho**4 &
                       & - 6.0_MK*c*r0**3*rho**5 &
                       & + 4.0_MK*c**2*r0**6*rho**2 &
                       & - 8.0_MK*c**2*r0**5*rho**3 &
                       & + 4.0_MK*c**2*r0**4*rho**4 &
                       & + 2.0_MK*c*r0**2*rho**6 &
                       & + 32.0_MK*r0**3*rho**3*px**2 &
                       & - 48.0_MK*r0**2*rho**4*px**2 &
                       & - 24.0_MK*r0**2*rho**2*px**4 &
                       & + 24.0_MK*r0*rho**5*px**2 &
                       & + 24.0_MK*r0*rho**3*px**4 &
                       & + 8.0_MK*r0*rho*px**6 &
                       & + 2.0_MK*c*r0**3*rho*px**4 &
                       & + 2.0_MK*c*r0**2*rho**2*px**4 &
                       & - 4.0_MK*c*r0**3*rho**3*px**2 &
                       & + 4.0_MK*c*r0**2*rho**4*px**2 &
                       & - rho**8 - px**8) &
                       & * (2.0_MK*r0*rho - rho**2 - px**2)**(-4)*rho**(-2)

              mesh(ilevel,ipatch)%vort(1,i,j,k) = 0.0_MK
              mesh(ilevel,ipatch)%vort(2,i,j,k) = - SIN(theta)*vort_mag
              mesh(ilevel,ipatch)%vort(3,i,j,k) = COS(theta)*vort_mag

            ELSE
              mesh(ilevel,ipatch)%vort(1,i,j,k) = 0.0_MK
              mesh(ilevel,ipatch)%vort(2,i,j,k) = 0.0_MK
              mesh(ilevel,ipatch)%vort(3,i,j,k) = 0.0_MK
            ENDIF

          END DO !i
        END DO !j
      END DO !k

    END DO !ipatch
  END DO !ilevel

!---------------------------------------------------------------------------------!
! Solve Poisson equation
!---------------------------------------------------------------------------------!
  CALL pmlib_poisson_solve( patch,topo_all,mesh,ierr, &
                          & reg_vort = .TRUE. )
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to solve Poisson equation.')
    GOTO 9999
  END IF

!---------------------------------------------------------------------------------!
! Calculate error
!---------------------------------------------------------------------------------!
! IF( ASSOCIATED( pltfld ) ) THEN
! DEALLOCATE( pltfld )
! END IF
! ALLOCATE( pltfld( 1, 1-nghost(1):ncell(1)+nghost(2), &
! & 1-nghost(3):ncell(2)+nghost(4), &
! & 1-nghost(5):ncell(3)+nghost(6) ) )

  error_vel = 0.0_MK
  int_vel = 0.0_MK

  error_vort = 0.0_MK
  int_vort = 0.0_MK

! max_vel = 3.3539478686618872E-4_MK
! max_vel = EXP(-c) ! stream
! max_vort = 4.0_MK*c*EXP(-c) ! vort
! max_err = 0.0_MK

  DO ilevel = 1,pmlib_nlevel
    DO ipatch = 1,pmlib_npatch(ilevel)

      xmin = topo_all(ilevel,ipatch)%cuboid(rank)%xmin
      dx = topo_all(ilevel,ipatch)%cuboid(rank)%dx
      ncell = topo_all(ilevel,ipatch)%cuboid(rank)%ncell
      nghost = topo_all(ilevel,ipatch)%cuboid(rank)%nghost

! DO k = 1-nghost(5),ncell(3)+nghost(6)
      DO k = 1,ncell(3)
        pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
! DO j = 1-nghost(3),ncell(2)+nghost(4)
        DO j = 1,ncell(2)
          py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
! DO i = 1-nghost(1),ncell(1)+nghost(2)
          DO i = 1,ncell(1)
            px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

            IF( mesh(ilevel,ipatch)%patch_mask(i,j,k) )THEN

            rho = SQRT(pz*pz+py*py)
            phi = SQRT((rho - r0)**2 + px**2)
            theta = ATAN2(pz,py)

            IF (phi .LT. r0) THEN
! velocity check
              vel_mag = 2.0_MK*c*r0**2*px * &
                  & EXP(-c*r0**2/(2.0_MK*r0*rho - rho**2 - px**2)) &
                  & *(2.0_MK*r0*rho - rho**2 - px**2)**(-2)

              vel(1) = EXP(-c*r0**2/(2.0_MK*r0*rho - rho**2 - px**2)) * &
                     & ( 4.0_MK*r0**2 *rho**2 &
                     & - 4.0_MK*r0*rho**3 &
                     & - 4.0_MK*r0*rho*px**2 &
                     & + rho**4 + 2.0_MK*rho**2*px**2 &
                     & + px**4 + 2.0_MK*c*r0**3*rho &
                     & - 2.0_MK*c*r0**2*rho**2) &
                     & *(2.0_MK*r0*rho - rho**2 - px**2)**(-2)*rho**(-1)
              vel(2) = COS(theta)*vel_mag
              vel(3) = SIN(theta)*vel_mag

              diff_vel = ( vel(1) - mesh(ilevel,ipatch)%vel(1,i,j,k) )**2 + &
                       & ( vel(2) - mesh(ilevel,ipatch)%vel(2,i,j,k) )**2 + &
                       & ( vel(3) - mesh(ilevel,ipatch)%vel(3,i,j,k) )**2

              vel_mag = vel(1)**2 + vel(2)**2 + vel(3)**2

! stream function check
! vel_mag = EXP(-c/(1 - (rho**2 + r0**2 + px**2 - 2*r0*rho )/r0**2))

! vel(1) = 0.0_MK
! vel(2) = - SIN(theta) *vel_mag
! vel(3) = COS(theta) *vel_mag

! vorticity check
              vort_mag = -EXP(- c*r0**2/(2.0_MK*r0*rho - rho**2 - px**2)) * &
                       & ( 4.0_MK*c**2*r0**4*px**2*rho**2 &
                       & - 16.0_MK*r0**4*rho**4 &
                       & + 32.0_MK*r0**3*rho**5 &
                       & - 24.0_MK*r0**2*rho**6 &
                       & + 8.0_MK*r0*rho**7 &
                       & - 4.0_MK*rho**6*px**2 &
                       & - 6.0_MK*rho**4*px**4 &
                       & - 4.0_MK*rho**2*px**6 &
                       & - 8.0_MK*c*r0**5*rho**3 &
                       & + 8.0_MK*c*r0**4*rho**4 &
                       & - 6.0_MK*c*r0**3*rho**5 &
                       & + 4.0_MK*c**2*r0**6*rho**2 &
                       & - 8.0_MK*c**2*r0**5*rho**3 &
                       & + 4.0_MK*c**2*r0**4*rho**4 &
                       & + 2.0_MK*c*r0**2*rho**6 &
                       & + 32.0_MK*r0**3*rho**3*px**2 &
                       & - 48.0_MK*r0**2*rho**4*px**2 &
                       & - 24.0_MK*r0**2*rho**2*px**4 &
                       & + 24.0_MK*r0*rho**5*px**2 &
                       & + 24.0_MK*r0*rho**3*px**4 &
                       & + 8.0_MK*r0*rho*px**6 &
                       & + 2.0_MK*c*r0**3*rho*px**4 &
                       & + 2.0_MK*c*r0**2*rho**2*px**4 &
                       & - 4.0_MK*c*r0**3*rho**3*px**2 &
                       & + 4.0_MK*c*r0**2*rho**4*px**2 &
                       & - rho**8 - px**8) &
                       & * (2.0_MK*r0*rho - rho**2 - px**2)**(-4)*rho**(-2)

              vort(1) = 0.0_MK
              vort(2) = - SIN(theta)*vort_mag
              vort(3) = COS(theta)*vort_mag

              diff_vort = ( vort(1) - mesh(ilevel,ipatch)%vort(1,i,j,k) )**2 + &
                        & ( vort(2) - mesh(ilevel,ipatch)%vort(2,i,j,k) )**2 + &
                        & ( vort(3) - mesh(ilevel,ipatch)%vort(3,i,j,k) )**2

              vort_mag = vort(1)**2 + vort(2)**2 + vort(3)**2

            ELSE
              vel_mag = 0.0_MK
              vort_mag = 0.0_MK

              diff_vel = (0.0_MK - mesh(ilevel,ipatch)%vel(1,i,j,k))**2 + &
                       & (0.0_MK - mesh(ilevel,ipatch)%vel(2,i,j,k))**2 + &
                       & (0.0_MK - mesh(ilevel,ipatch)%vel(3,i,j,k))**2

              diff_vort = (0.0_MK - mesh(ilevel,ipatch)%vort(1,i,j,k))**2 + &
                        & (0.0_MK - mesh(ilevel,ipatch)%vort(2,i,j,k))**2 + &
                        & (0.0_MK - mesh(ilevel,ipatch)%vort(3,i,j,k))**2
            END IF


            mesh(ilevel,ipatch)%vort(1,i,j,k) = diff_vel
            mesh(ilevel,ipatch)%vort(2,i,j,k) = 0.0_MK
            mesh(ilevel,ipatch)%vort(3,i,j,k) = 0.0_MK

            error_vel = error_vel + diff_vel * dx(1) * dx(2) * dx(3)
            int_vel = int_vel + vel_mag * dx(1) * dx(2) * dx(3)

            error_vort = error_vort + diff_vort * dx(1) * dx(2) * dx(3)
            int_vort = int_vort + vort_mag * dx(1) * dx(2) * dx(3)

            END IF

          END DO !i
        END DO !j
      END DO !k

    END DO !ipatch
  END DO !ilevel

!---------------------------------------------------------------------------------!
! Gather and calculate global error
!---------------------------------------------------------------------------------!
  CALL MPI_REDUCE(error_vel,global_error_vel(n),1,mpi_prec_real, &
       & MPI_SUM,0, comm, ierr)
  CALL MPI_REDUCE(int_vel,global_int_vel,1,mpi_prec_real, &
       & MPI_SUM,0, comm, ierr)

  CALL MPI_REDUCE(error_vort,global_error_vort(n),1,mpi_prec_real, &
       & MPI_SUM,0, comm, ierr)
  CALL MPI_REDUCE(int_vort,global_int_vort,1,mpi_prec_real, &
       & MPI_SUM,0, comm, ierr)

  IF(rank .EQ. 0)THEN
    global_error_vel(n) = SQRT(global_error_vel(n)/global_int_vel)
    global_error_vort(n) = SQRT(global_error_vort(n)/global_int_vort)
    dxs(n) = MAXVAL(patch(pmlib_nlevel,1)%dx)
  END IF

  END DO

  IF(rank .EQ. 0)THEN
    WRITE(*,*)'Convergence test: '

    DO n = 1,nn
      WRITE(*,'(2E20.12)') dxs(n), global_error_vel(n)!, global_error_vort(n)
    END DO

    OPEN(10,file = 'convergence.dat')
    DO n = 1,nn
      WRITE(10,'(2E20.12)') dxs(n), global_error_vel(n)!, global_error_vort(n)
    END DO
    CALL SYSTEM('gnuplot convergence.gnu')
  END IF


! PLOT
! CALL output_mesh( 'mesh',patch(ilevel,ipatch), &
! & topo_all(ilevel,ipatch)%cuboid, &
! & mesh(ilevel,ipatch),ierr)

!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
! DEALLOCATE( pltfld )

!---------------------------------------------------------------------------------!
! Finalise MPI
!---------------------------------------------------------------------------------!
  CALL MPI_finalize(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to finalize MPI.')
  END IF
  GOTO 1111

!---------------------------------------------------------------------------------!
! Return·
!---------------------------------------------------------------------------------!
 9999 CONTINUE
  CALL MPI_ABORT(comm,ierr)

 1111 CONTINUE
END PROGRAM multipatch
