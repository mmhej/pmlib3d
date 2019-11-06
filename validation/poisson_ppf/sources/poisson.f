!---------------------------------------------------------------------------------!
! poisson.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
! <<HEADER>>
!---------------------------------------------------------------------------------!
PROGRAM poisson

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_poisson
USE pmlib_mod_output

USE mod_poisson

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER                           :: ndim  = pmlib_ndim
  INTEGER, PARAMETER                           :: nvort = pmlib_nvort

  CHARACTER(LEN=7), PARAMETER                  :: caller = 'poisson'
  REAL(MK),PARAMETER                           :: PI = ACOS(-1.0_MK)
  INTEGER                                      :: ierr

  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  REAL(MK)                                     :: r,r0,c

  CHARACTER(LEN=256)                           :: nodename
  INTEGER                                      :: name_len

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  CHARACTER(LEN=256)                           :: pltfile

  REAL(MK)                                     :: max_vel, max_vort, vel_max
  REAL(MK)                                     :: error_vel, error_vort
  REAL(MK)                                     :: diff_vel, diff_vort
  REAL(MK),DIMENSION(ndim)                     :: vel, vel_inf, vort

  INTEGER                                      :: n, nn
  REAL(MK),DIMENSION(:),POINTER                :: global_error_vel
  REAL(MK),DIMENSION(:),POINTER                :: global_error_vort
  REAL(MK),DIMENSION(:),POINTER                :: dxs
  INTEGER,DIMENSION(:),POINTER                 :: nx

  REAL(MK)                                     :: global_int_vel, global_int_vort
  REAL(MK)                                     :: int_vel, int_vort

  REAL(MK)                                     :: rho, phi, theta
  REAL(MK)                                     :: vort_mag, vel_mag
  REAL(MK)                                     :: sigma, circ

  LOGICAL                                      :: regularise_vort

  INTEGER                                      :: comm
  INTEGER                                      :: rank
  INTEGER                                      :: nproc

  REAL(MK),DIMENSION(:,:,:,:),POINTER          :: pltfld => NULL()

!---------------------------------------------------------------------------------!
! Construct classes
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                            :: patch
  TYPE(class_topology_all)                     :: topo_all
  TYPE(class_mesh),SAVE                        :: mesh

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

  mpi_comm  = comm
  mpi_rank  = rank
  mpi_nproc = nproc

!---------------------------------------------------------------------------------!
! Set up parameters
!---------------------------------------------------------------------------------!
  vel_inf = (/ 0.0_MK, 0.0_MK, 0.0_MK /)

  patch%level  = 1
  patch%parent = (/ 0, 0 /)
  patch%ptype  = 1
  patch%xmin   = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
  patch%xmax   = (/  1.0_MK,  1.0_MK,  1.0_MK /)
  patch%bound_cond = (/ 1, 1, 0 /)
  patch%nghost     = 4

  regularise_vort = .TRUE.

  pmlib_regularisation_radius  = 2.0_MK
  pmlib_poisson_order          = 0
  pmlib_poisson_kernel         = 1
  pmlib_time_integration_order = 3
  pmlib_fd_order               = 4
  pmlib_interpolation_order    = 4

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
!  nx = (/ 32 /) 
  nx = (/ 16, 32 /) 
!  nx = (/ 16, 32, 64 /)
!  nx = (/ 16, 32, 64, 128 /)
!  nx = 32

!---------------------------------------------------------------------------------!
! Loop through different resolutions
!---------------------------------------------------------------------------------!
  DO n = 1,nn
    patch%dx = 1.0_MK/REAL(nx(n),MK)
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
    patch%fft_seq = (/ 1 , 2 , 3 /)
!    patch%fft_seq = (/ 3 , 1 , 2 /)
!    patch%fft_seq = (/ 2 , 3 , 1 /)

!    patch%fft_seq = (/ 2 , 1 , 3 /)
!    patch%fft_seq = (/ 3 , 2 , 1 /)
!    patch%fft_seq = (/ 1 , 3 , 2 /)

!---------------------------------------------------------------------------------!
! Topology create
!---------------------------------------------------------------------------------!
    CALL pmlib_topology_cuboid(patch,topo_all%cuboid,ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to create topology.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Allocate mesh set
!---------------------------------------------------------------------------------!
    CALL pmlib_mesh_allocate( topo_all%cuboid,mesh,ierr, &
                            & vort = .TRUE., dvort = .TRUE., &
                            & vel = .TRUE., mask = .FALSE.)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to allocate mesh.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Setup Greens function
!---------------------------------------------------------------------------------!
    IF( patch%bound_cond(1) .EQ. 1 .AND. &
      & patch%bound_cond(2) .EQ. 1 .AND. &
      & patch%bound_cond(3) .EQ. 0)THEN
      CALL poisson_setup_ppf(patch,topo_all,mesh,ierr)
    ELSE
      CALL pmlib_poisson_setup(patch,topo_all,mesh,ierr)
    END IF
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Setup initial vorticity field
!---------------------------------------------------------------------------------!
  dx     = topo_all%cuboid(rank)%dx
  xmin   = topo_all%cuboid(rank)%xmin
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

  max_vort = 0.0_MK
  r0       = 1.0_MK
  c        = 10.0_MK

  mesh%vort = 0.0_MK
  mesh%vel = 0.0_MK

  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)

    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)

      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)


        IF ( ABS(pz) .LT. r0) THEN
          mesh%vort(1,i,j,k) = ( EXP(c*pz**2/((pz - 1.0_MK)*(pz + 1.0_MK))) &
                             &  * SIN(PI*px) &
                             &  * SIN(PI*py) &
                             &  * (  2.0_MK*PI**2* pz**8 &
                             &    -  8.0_MK*PI**2* pz**6 &
                             &    + 12.0_MK*PI**2* pz**4 &
                             &    -  6.0_MK    *c* pz**4 &
                             &    -  8.0_MK*PI**2* pz**2 &
                             &    -  4.0_MK *c**2* pz**2 &
                             &    +  4.0_MK    *c* pz**2 &
                             &    +  2.0_MK*PI**2        &
                             &    +  2.0_MK    *c )      &
                             & )/( (pz - 1.0_MK)**4 * (pz + 1.0_MK)**4 )
          mesh%vort(2,i,j,k) = 0.0_MK
          mesh%vort(3,i,j,k) = 0.0_MK
        ELSE
          mesh%vort(1,i,j,k) = 0.0_MK
          mesh%vort(2,i,j,k) = 0.0_MK
          mesh%vort(3,i,j,k) = 0.0_MK
        ENDIF
        IF( ABS( mesh%vort(1,i,j,k) ) .GT. max_vort ) THEN
          max_vort = ABS( mesh%vort(1,i,j,k) )
        END IF

      END DO !i
    END DO !j
  END DO !k

!---------------------------------------------------------------------------------!
! Solve poisson
!---------------------------------------------------------------------------------!
  CALL pmlib_poisson_solve( patch,topo_all,mesh,ierr, &
                          & reg_vort = .TRUE., reproj = .FALSE.)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to solve the Poisson equation.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Calculate error
!---------------------------------------------------------------------------------!
!  IF( ASSOCIATED( pltfld ) ) THEN
!    DEALLOCATE( pltfld )
!  END IF  
!  ALLOCATE( pltfld(1,1-nghost(1):ncell(1)+nghost(2), &
!                   &  1-nghost(3):ncell(2)+nghost(4), &
!                   &  1-nghost(5):ncell(3)+nghost(6) ) )

  error_vel = 0.0_MK
  int_vel   = 0.0_MK

  error_vort = 0.0_MK
  int_vort   = 0.0_MK

!  max_vel = 0.3139696888_MK*c*EXP(-1.023017903_MK*c) ! velocity
!  max_vel = 3.3539478686618872E-4_MK
!  max_vel = EXP(-c) ! stream
!  max_vort = 4.0_MK*c*EXP(-c) ! vort

  xmin   = topo_all%cuboid(rank)%xmin
  dx     = topo_all%cuboid(rank)%dx
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost

!  DO k = 1-nghost(5),ncell(3)+nghost(6)
  DO k = 1,ncell(3)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
!    DO j = 1-nghost(3),ncell(2)+nghost(4)
    DO j = 1,ncell(2)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
!      DO i = 1-nghost(1),ncell(1)+nghost(2)
      DO i = 1,ncell(1)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)

!        rho   = SQRT(px*px+py*py)
!        phi   = SQRT((rho - r0)**2 + pz**2)
!        theta = ATAN2(py,px)

        rho   = SQRT(pz*pz+py*py)
        phi   = SQRT((rho - r0)**2 + px**2)
        theta = ATAN2(pz,py)

        IF (ABS(pz) .LT. r0) THEN
! velocity check
          vel(1) = 0.0_MK

          vel(2) = - 2.0_MK*c*pz*EXP(c*pz**2/((pz - 1.0_MK)*(pz + 1.0_MK))) &
                 &  * SIN(PI*px)* SIN(PI*py)/( (pz - 1.0_MK)**2 * (pz + 1.0_MK)**2 )

          vel(3) = - PI*EXP(c*pz**2/((pz - 1.0_MK)*(pz + 1.0_MK)))*SIN(PI*px)*COS(PI*py)

          diff_vel = ( vel(1) - mesh%vel(1,i,j,k) )**2 + &
                   & ( vel(2) - mesh%vel(2,i,j,k) )**2 + &
                   & ( vel(3) - mesh%vel(3,i,j,k) )**2

          vel_mag = vel(1)**2 + vel(2)**2 + vel(3)**2

! stream function check
!          vort_mag = EXP(-c/(1 - (rho**2 + r0**2 + px**2 - 2*r0*rho )/r0**2))

!          vel(1) =   0.0_MK
!          vel(2) = - SIN(theta) *vort_mag 
!          vel(3) =   COS(theta) *vort_mag

! vorticity check
          vort(1) = ( EXP(c*pz**2/((pz - 1.0_MK)*(pz + 1.0_MK))) &
                  &  * SIN(PI*px) &
                  &  * SIN(PI*py) &
                  &  * (  2.0_MK*PI**2* pz**8 &
                  &    -  8.0_MK*PI**2* pz**6 &
                  &    + 12.0_MK*PI**2* pz**4 &
                  &    -  6.0_MK    *c* pz**4 &
                  &    -  8.0_MK*PI**2* pz**2 &
                  &    -  4.0_MK *c**2* pz**2 &
                  &    +  4.0_MK    *c* pz**2 &
                  &    +  2.0_MK*PI**2        &
                  &    +  2.0_MK    *c )      &
                  & )/( (pz - 1.0_MK)**4 * (pz + 1.0_MK)**4 )
          vort(2) = 0.0_MK
          vort(3) = 0.0_MK

          diff_vort = ( vort(1) - mesh%vort(1,i,j,k) )**2 + &
                    & ( vort(2) - mesh%vort(2,i,j,k) )**2 + &
                    & ( vort(3) - mesh%vort(3,i,j,k) )**2

          vort_mag = vort(1)**2 + vort(2)**2 + vort(3)**2

        ELSE
          vel_mag  = 0.0_MK
          vort_mag = 0.0_MK

          diff_vel = (0.0_MK - mesh%vel(1,i,j,k))**2 + &
                   & (0.0_MK - mesh%vel(2,i,j,k))**2 + &
                   & (0.0_MK - mesh%vel(3,i,j,k))**2

          diff_vort = (0.0_MK - mesh%vort(1,i,j,k))**2 + &
                    & (0.0_MK - mesh%vort(2,i,j,k))**2 + &
                    & (0.0_MK - mesh%vort(3,i,j,k))**2
        END IF


!        pltfld(1,i,j,k) = SQRT( vel(1)**2 + vel(2)**2 + vel(3)**2 )

        error_vel = error_vel   +  diff_vel * dx(1) * dx(2) * dx(3)
        int_vel   = int_vel     +   vel_mag * dx(1) * dx(2) * dx(3)

        error_vort = error_vort + diff_vort * dx(1) * dx(2) * dx(3)
        int_vort   = int_vort   +  vort_mag * dx(1) * dx(2) * dx(3)

      END DO !i
    END DO !j
  END DO !k 

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
    global_error_vel(n)  = SQRT(global_error_vel(n)/global_int_vel)
    global_error_vort(n) = SQRT(global_error_vort(n)/global_int_vort)
    dxs(n) = MAXVAL(patch%dx)
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
!  WRITE(pltfile,'(A)') 'err'
!  CALL pmlib_output_vtk(topology(rank,1,1),pltfld,pltfile,ierr,ghost = .TRUE.)

!  WRITE(pltfile,'(A)') 'vort'
!  CALL pmlib_output_vtk(topology(rank),mesh%vort,pltfile,ierr,ghost = .TRUE.)

!  WRITE(pltfile,'(A)') 'out'
!  CALL pmlib_output_vtk(topology(rank),mesh%vel,pltfile,ierr,ghost = .TRUE.)


!---------------------------------------------------------------------------------!
! Deallocate local pointers
!---------------------------------------------------------------------------------!
!  DEALLOCATE( pltfld )

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
END PROGRAM poisson
