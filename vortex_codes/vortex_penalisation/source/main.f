!---------------------------------------------------------------------------------!
! pmlib client: vortex_penalisation
!---------------------------------------------------------------------------------!
PROGRAM vortex_penalisation

USE mod_parameters
USE mod_input
USE mod_output
USE mod_flowcases
USE mod_particle_derivatives
USE mod_penalisation
USE mod_diagnostics

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_fourier
USE pmlib_mod_mesh
USE pmlib_mod_particles
USE pmlib_mod_regularise
USE pmlib_mod_remesh
USE pmlib_mod_repatch
USE pmlib_mod_interpolation
!USE pmlib_mod_poisson
USE pmlib_mod_output
USE pmlib_mod_visualise

USE poisson_solver_module

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER                                      :: IARGC
  CHARACTER(LEN=6), PARAMETER                  :: caller = 'main'
  INTEGER,DIMENSION(8)                         :: timedate
  CHARACTER(LEN=256)                           :: infile
  CHARACTER(LEN=256)                           :: outfile
  CHARACTER(LEN=256)                           :: inarg
  REAL(MK)                                     :: input_trunc
  LOGICAL                                      :: file_exist

  INTEGER                                      :: ierr
  INTEGER                                      :: i,j,k
  CHARACTER(LEN=256)                           :: nodename
  INTEGER                                      :: name_len
  INTEGER                                      :: itimestage
  INTEGER                                      :: extend
  REAL(MK)                                     :: dtime_in
  REAL(MK),DIMENSION(ndim)                     :: foc
  REAL(MK)                                     :: dist,azi,ele,zoom,vmax
  REAL(MK)                                     :: timing_start, timing_stop

  REAL(MK),DIMENSION(ndim)                     :: xmin, dx
  INTEGER,DIMENSION(ndim)                      :: ncell, offset
  INTEGER,DIMENSION(2*ndim)                    :: nghost

!---------------------------------------------------------------------------------!
! Construct classes
!---------------------------------------------------------------------------------!
  TYPE(class_patch)                          :: patch
  TYPE(class_topology_all)                   :: topo_all

  TYPE(class_mesh),SAVE                      :: mesh
  TYPE(class_particles),SAVE                 :: part

!---------------------------------------------------------------------------------!
! Initialise program
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
  CALL MPI_BARRIER(comm,ierr)

  mpi_comm  = comm
  mpi_rank  = rank
  mpi_nproc = nproc

!---------------------------------------------------------------------------------!
! Display starting date and a cool logo
!---------------------------------------------------------------------------------!
  CALL MPI_BARRIER(comm,ierr)
  IF (rank .EQ. 0) THEN
    CALL DATE_AND_TIME(values = timedate)
    WRITE(*,'(A)')""
    WRITE(*,'(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)') &
     & 'STARTING SIMULATION ON: ', &
     & timedate(3),'-',timedate(2),'-',timedate(1), &
     & '  at  ',timedate(5),':',timedate(6),':',timedate(7)
    WRITE(*,'(A,I4)') 'Number of processors initiated: ', nproc
    WRITE(*,*) ""
  ENDIF
  CALL MPI_BARRIER(comm,ierr)

!---------------------------------------------------------------------------------!
! Check input setup file
!---------------------------------------------------------------------------------!
  IF (IARGC() .EQ. 1) THEN
    CALL GETARG(1,infile)
    INQUIRE(FILE = TRIM(infile), EXIST = file_exist)
    IF(.NOT. file_exist)THEN
      ierr = -1
      CALL pmlib_write(rank,caller,'Setup file does not exist.')
      GOTO 9999
    END IF
  ELSE
    ierr = - 1
    CALL pmlib_write(rank,caller,'Specify a setup file.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Read setup file
!---------------------------------------------------------------------------------!
  CALL input_setup(infile,patch,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to read setup file.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Check that pmlib parameters have been set
!---------------------------------------------------------------------------------!
  CALL pmlib_parameters_check(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'pmlib parameters are not set correctly.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Setup initial parameters
!---------------------------------------------------------------------------------!
  itime      = 0
  time       = 0.0_MK
  idiag      = 0
  IF(dtime .GT. 1.0E-10_MK)THEN
    dtime_in = dtime
  ELSE
    dtime_in = 0.1_MK
  END IF
!---------------------------------------------------------------------------------!
! Adjust patches
!---------------------------------------------------------------------------------!
  CALL pmlib_patch_adjust(patch,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to adjust patches.')
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
! Create cuboid topology
!---------------------------------------------------------------------------------!
  CALL pmlib_topology_cuboid(patch,topo_all%cuboid,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to create topology.')
    GOTO 9999
  ENDIF

  IF(restart)THEN
!---------------------------------------------------------------------------------!
! Check restart file
!---------------------------------------------------------------------------------!
    INQUIRE(FILE = TRIM(restart_file), EXIST = file_exist)
    IF(.NOT. file_exist)THEN
      CALL pmlib_write(rank,caller,'Particle file does not exist.')
      GOTO 9999
    END IF

!---------------------------------------------------------------------------------!
! Read particles from file
!---------------------------------------------------------------------------------!
    CALL pmlib_write(rank,caller,'Re-starting',verb=.TRUE.)
    CALL input_particles( TRIM(restart_file), &
                               & patch,topo_all,mesh,part,0.0_MK,ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed input particles.')
      GOTO 9999
    ENDIF

    ntime = itime + ntime
  ELSE

!---------------------------------------------------------------------------------!
! Allocate mesh set
!---------------------------------------------------------------------------------!
    CALL pmlib_mesh_allocate(topo_all%cuboid,mesh,ierr, &
                            & vort = .TRUE., dvort = .TRUE., &
                            & vel = .TRUE., mask = .FALSE.)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to allocate mesh.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Setup Poisson solver
!---------------------------------------------------------------------------------!
!    CALL pmlib_poisson_setup(patch,topo_all,mesh,ierr)
!    IF (ierr .NE. 0) THEN
!      CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
!      GOTO 9999
!    ENDIF

! Initialise 2 poisson solvers
    CALL poisson_solver_initialise( 2 )

! Store cuboid partition to solver 1
    ALLOCATE( poisson_solver(1)%partition( 0:nproc-1 ) ) 
    offset = (/ 1, 1, 1 /)
    DO i = 0,nproc-1
      poisson_solver(1)%partition(i)%ncell = topo_all%cuboid(i)%ncell
      poisson_solver(1)%partition(i)%icell = topo_all%cuboid(i)%icell-offset
      poisson_solver(1)%partition(i)%dx    = topo_all%cuboid(i)%dx
    END DO

! Setup solver 1
    CALL poisson_solver_setup3d(1,patch%ncell,patch%bound_cond,patch%dx)
    CALL poisson_solver_set_return_curl(1,.TRUE.) ! specify lhs operator

!---------------------------------------------------------------------------------!
! Set up flow case
!---------------------------------------------------------------------------------!
    SELECT CASE ( TRIM(flowcase) )
      CASE ('sphere')
        CALL flowcases_sphere(patch,ierr)
      CASE DEFAULT
        ierr = -1
        CALL pmlib_write(rank,caller, 'Flowcase unknown.')
        GOTO 9999
    END SELECT
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to set up flow case.')
      GOTO 9999
    ENDIF

! CALL diagnostics_spectral(patch,topo_all,mesh,ierr)

!---------------------------------------------------------------------------------!
! Create particles
!---------------------------------------------------------------------------------!
    CALL pmlib_remesh(patch,topo_all%cuboid,mesh,part,0.0_MK,0.0_MK,ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to create vorticity particles.')
      GOTO 9999
    ENDIF
  END IF

!---------------------------------------------------------------------------------!
! Output job parameters
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0)THEN
    WRITE(*,'(A)')'(hagall): Job parameters'
    WRITE(*,'(2A)')  'Flowcase:                      ', TRIM(flowcase)
    IF(visc .EQ. 0.0_MK)THEN
      WRITE(*,'(A)')'Reynolds number:               infinity'
    ELSE
      WRITE(*,'(A,I7)')'Reynolds number:         ', NINT(1.0_MK/visc)
    END IF
    WRITE(*,'(A,I7)')'Resolution ncell/L:      ', NINT(1.0_MK/patch%dx(1))
    WRITE(*,'(A,F8.4)')'Time step size:              ', dtime_in
    WRITE(*,'(A,F8.4)')'Time parameter:              ', ctime
    WRITE(*,'(A,L )')'Re-projection:                 ', reproject
    WRITE(*,'(A,L )')'Regularise vorticity:          ', regularise_vort
    WRITE(*,'(A,I7,F8.4)')'Re-mesh interval/trunc:  ', iremesh, remesh_trunc
    WRITE(*,'(A,I7,F8.4)')'Re-patch interval/trunc: ', irepatch, repatch_trunc
    WRITE(*,'(A,I7)')'Regularisation order:    ', pmlib_poisson_order
    WRITE(*,'(A,I7)')'Time integration order:  ', pmlib_time_integration_order
    WRITE(*,'(A,I7)')'Finite difference order: ', pmlib_fd_order
    WRITE(*,'(A,I7)')'Interpolation order:     ', pmlib_interpolation_order
  END IF

!---------------------------------------------------------------------------------!
! Delete ABORT file if exists
!---------------------------------------------------------------------------------!
  INQUIRE(FILE='ABORT',EXIST=file_exist)
  IF (file_exist .AND. rank .EQ. 0) THEN
    OPEN(10,FILE='ABORT')
    CLOSE(10,STATUS='delete')
  ENDIF

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
! Simulate
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
  IF(rank .EQ. 0)THEN
    WRITE(*,'(A)')''
    WRITE(*,'(A)')'SIMULATE'
    WRITE(*,'(A)')''
  END IF

  DO WHILE(itime .LE. ntime)
!---------------------------------------------------------------------------------!
! Set time step toggles
!---------------------------------------------------------------------------------!
    remesh  = .FALSE.
    repatch = .FALSE.
    IF(iremesh .NE. 0)THEN
      IF (MOD(itime,iremesh) .EQ. 0) THEN
        remesh = .TRUE.
      END IF
    END IF
    IF(irepatch .NE. 0)THEN
      IF (MOD(itime,irepatch) .EQ. 0) THEN
        repatch = .TRUE.
        remesh  = .TRUE.
      END IF
    END IF
    INQUIRE(FILE='ABORT',EXIST=abort)

!---------------------------------------------------------------------------------!
! Set time step size
!---------------------------------------------------------------------------------!
    IF(ctime .NE. 0.0_MK)THEN
      dtime = MIN( dtime_in, ctime/vort_max )
    END IF
    IF (dtime .LT. 1.0E-10_MK) THEN
      CALL pmlib_write(rank,caller,'time step is too small.')
      GOTO 9999
    ENDIF
    IF (dtime .GT. 1.0E0_MK) THEN
      CALL pmlib_write(rank,caller,'time step is too large.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! output time step data
!---------------------------------------------------------------------------------!
    IF(rank .EQ. 0) THEN
      WRITE(*,'(A,I5,A,F8.3,A,F6.4,A,F6.2)') ' Timestep: ',itime, &
           & '   Time: ', time,'   dt: ', dtime, '   max. vort: ', vort_max
    ENDIF

!---------------------------------------------------------------------------------!
! Do time integration stages
!---------------------------------------------------------------------------------!
    DO itimestage = 1,pmlib_time_integration_order

!---------------------------------------------------------------------------------!
! Interpolate particle vorticity to mesh
!---------------------------------------------------------------------------------!
      CALL pmlib_interp_particle_mesh( topo_all%cuboid, part%pos, part%vort, &
                                     & mesh%vort, ierr, clear = .TRUE.)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(rank,caller, &
                        & 'Failed to interpolate particle vorticity to mesh.')
        GOTO 9999
      ENDIF

!---------------------------------------------------------------------------------!
! Re-patch if specified
!---------------------------------------------------------------------------------!
      IF( repatch .AND. itimestage .EQ. 1 )THEN
        extend = NINT( 0.25_MK/MAXVAL(patch%dx) )

        CALL pmlib_remesh( patch,topo_all%cuboid,mesh,part,0.0_MK,0.0_MK,ierr )
        CALL pmlib_repatch( patch,topo_all,part,mesh, & 
                          & repatch_trunc*vort_max ,extend,ierr )
        CALL pmlib_interp_particle_mesh( topo_all%cuboid,part%pos,part%vort, &
                                       & mesh%vort, ierr, clear = .TRUE.)
        IF (ierr .NE. 0) THEN
          CALL pmlib_write(rank,caller,'Failed to re-patch.')
          GOTO 9999
        ENDIF
      END IF

!---------------------------------------------------------------------------------!
! Calculate right-hand side of the trajectory equation
!---------------------------------------------------------------------------------!
      IF(reproject .OR. (itime .EQ. 0 .AND. itimestage .EQ. 1) )THEN
        CALL poisson_solver_set_reprojection(1,.TRUE.)
        CALL poisson_solver_push(1,offset,mesh%vort)
        CALL poisson_solver_solve3d(1)
        CALL poisson_solver_pull(1,offset,mesh%vel,mesh%vort)
! Map ghost cells
        DO i = 1,ndim
          CALL pmlib_mesh_map_ghost(topo_all%cuboid,i,6,ierr,incl_edges = .FALSE.)
          CALL pmlib_comm_pack(mesh%vort,ierr)
          CALL pmlib_comm_pack(mesh%vel,ierr)
          CALL pmlib_comm_send(ierr)
          CALL pmlib_comm_unpack(topo_all%cuboid,mesh%vel,0,ierr,clear=.FALSE.)
          CALL pmlib_comm_unpack(topo_all%cuboid,mesh%vort,0,ierr,clear=.FALSE.)
          CALL pmlib_comm_finalise(ierr)
        END DO
      ELSE
        CALL poisson_solver_set_reprojection(1,.FALSE.)
        CALL poisson_solver_push(1,offset,mesh%vort)
        CALL poisson_solver_solve3d(1)
        CALL poisson_solver_pull(1,offset,mesh%vel,mesh%vort)
! Map ghost cells
        DO i = 1,ndim
          CALL pmlib_mesh_map_ghost(topo_all%cuboid,i,3,ierr,incl_edges = .FALSE.)
          CALL pmlib_comm_pack(mesh%vel,ierr)
          CALL pmlib_comm_send(ierr)
          CALL pmlib_comm_unpack(topo_all%cuboid,mesh%vel,0,ierr,clear=.FALSE.)
          CALL pmlib_comm_finalise(ierr)
        END DO
      END IF
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
        GOTO 9999
      ENDIF
! add free stream velocity
      IF( vel_inf(1) .NE. 0.0_MK .OR. & 
        & vel_inf(2) .NE. 0.0_MK .OR. & 
        & vel_inf(3) .NE. 0.0_MK)THEN
        nghost = topo_all%cuboid(rank)%nghost
        ncell  = topo_all%cuboid(rank)%ncell
        DO k = 1-nghost(5),ncell(3)+nghost(6)
          DO j = 1-nghost(3),ncell(2)+nghost(4)
            DO i = 1-nghost(1),ncell(1)+nghost(2)
              mesh%vel(1,i,j,k) = mesh%vel(1,i,j,k) + vel_inf(1)
              mesh%vel(2,i,j,k) = mesh%vel(2,i,j,k) + vel_inf(2)
              mesh%vel(3,i,j,k) = mesh%vel(3,i,j,k) + vel_inf(3)
            END DO
          END DO
        END DO
      END IF

!---------------------------------------------------------------------------------!
! Enforce solid body by iterative penalisation
!---------------------------------------------------------------------------------!
      IF( itimestage .EQ. 1 .AND. penalisation)THEN
        CALL penalisation_solve(patch,topo_all,mesh,ierr)
        IF (ierr .NE. 0) THEN
          CALL pmlib_write(rank,caller, 'Failed to enforce solid body.')
          GOTO 9999
        ENDIF
      END IF

!---------------------------------------------------------------------------------!
! Calculate particle derivative of the vorticity
!---------------------------------------------------------------------------------!
      CALL part_deriv_vorticity(topo_all%cuboid,mesh,ierr)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(rank,caller, &
              & 'Failed to calculate particle derivative of the vorticity.')
        GOTO 9999
      ENDIF

!---------------------------------------------------------------------------------!
! Interpolate to particles
!---------------------------------------------------------------------------------!
      IF( (remesh .OR. itime .EQ. 0) .AND. itimestage .EQ. 1 )THEN
!      IF( remesh .AND. itimestage .EQ. 1 )THEN
        CALL pmlib_remesh( patch,topo_all%cuboid,mesh,part,remesh_trunc, & 
                         & remesh_trunc/dtime,ierr, & 
                         & dvort=.TRUE.,vel=.TRUE., mask = .FALSE.)
      ELSE
        CALL pmlib_interp_mesh_particle( topo_all%cuboid,mesh,part,ierr, &
                             & vort=.FALSE.,dvort=.TRUE.,vel=.TRUE., &
                             & clear = .TRUE.)
      END IF
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(rank,caller, &
                        & 'Failed to interpolate rhs to particles.')
        GOTO 9999
      END IF

!---------------------------------------------------------------------------------!
! Diagnostics 
! (here to ensure that the velocity, vorticity and position correspond)
!---------------------------------------------------------------------------------!
      IF(itimestage .EQ. 1)THEN
        SELECT CASE ( TRIM(flowcase) )
          CASE ('sphere')
            CALL diagnostics_sphere(topo_all%cuboid,part,ierr)
        END SELECT
      END IF

!---------------------------------------------------------------------------------!
! Advance particles - do time stepping
!---------------------------------------------------------------------------------!
      CALL pmlib_particles_advance( patch,topo_all%cuboid,part, &
                                  & itimestage,dtime,ierr,euler_vort=.FALSE.)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(rank,caller,'Failed to advance particles.')
        GOTO 9999
      END IF

    END DO !time stages

!---------------------------------------------------------------------------------!
! Plot
!---------------------------------------------------------------------------------!
    IF (iplot_mesh .NE. 0)THEN
      IF( MOD(itime,iplot_mesh) .EQ. 0)THEN

        WRITE(outfile,'(A,I5.5)') './output/mesh_I',itime
        CALL output_mesh(outfile,patch,topo_all%cuboid,mesh,ierr)

      END IF
    END IF
! Particles
    IF (iplot_part .NE. 0)THEN
      IF( MOD(itime,iplot_part) .EQ. 0 )THEN
! Set camera view
        foc  = vort_centroid 
        dist = 4.0_MK
        ele  = 5.0_MK
        zoom = 2.0_MK
        IF(itime .EQ. 0)THEN
          vmax = vort_max
        END IF
        azi  = 40.0_MK
!        WRITE(outfile,'(A,I5.5)') 'view_I',itime
!        WRITE(outfile,'(A)') 'view'
        WRITE(outfile,'(A,I5.5,A,I3.3)') './output/view_I',itime,'_A',0
        CALL pmlib_visualise_particles(outfile,topo_all%cuboid,part,vmax, &
                           & foc,dist,azi,ele,zoom,ierr,scl = 0.5_MK)

      END IF
    END IF

!---------------------------------------------------------------------------------!
! Output particles
!---------------------------------------------------------------------------------!
    IF (ioutput_part .NE. 0)THEN
      IF( MOD(itime,ioutput_part) .EQ. 0 .OR. abort)THEN

        WRITE(outfile,'(A,I5.5)') './output/part_I',itime
        CALL output_particles(outfile,topo_all%cuboid,part,ierr)

      END IF
    END IF

!---------------------------------------------------------------------------------!
! Increment time
!---------------------------------------------------------------------------------!
    itime = itime + 1
    time = time + dtime

!---------------------------------------------------------------------------------!
! If abort break the time loop
!---------------------------------------------------------------------------------!
    IF(abort)THEN
      CALL pmlib_write(rank,caller,'Abort file found - aborting.',verb=.TRUE.)
      GO TO 1111
    END IF


  END DO !itime

 1111 CONTINUE

!---------------------------------------------------------------------------------!
! Finalise poisson solver
!---------------------------------------------------------------------------------!
  CALL poisson_solver_finalise(1)
  CALL poisson_solver_finalise(2)

!---------------------------------------------------------------------------------!
! Finalise MPI
!---------------------------------------------------------------------------------!
  CALL MPI_finalize(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to finalize MPI.')
  END IF
  GOTO 2222

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
  CALL MPI_ABORT(comm,ierr)

 2222 CONTINUE
END PROGRAM vortex_penalisation
