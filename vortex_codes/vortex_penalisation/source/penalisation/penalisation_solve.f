!---------------------------------------------------------------------------------!
! penalisation_solve.f
!---------------------------------------------------------------------------------!
! This routine enforces a solid body by iterative Brinkman penalisation
!---------------------------------------------------------------------------------!
SUBROUTINE penalisation_solve(patch,topo_all,mesh,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_poisson
USE pmlib_mod_communication
USE pmlib_mod_output

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
  CHARACTER(LEN=18)                     :: caller = 'penalisation'
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  REAL(MK),DIMENSION(ndim)                     :: xmin, xmax, dx
  INTEGER,DIMENSION(ndim)                      :: ncell, offset
  INTEGER,DIMENSION(2*ndim)                    :: nghost

  REAL(MK),DIMENSION(ndim)                     :: force, dforce, sum_dforce

  REAL(MK)                                     :: sum_vel_pen, vel_pen
  REAL(MK)                                     :: sum_enst_pen, enst_pen
  REAL(MK),DIMENSION(nvort)                    :: dvort

  INTEGER                                      :: npen
  REAL(MK)                                     :: residual_vel, res_vel
  REAL(MK)                                     :: residual_force
  REAL(MK)                                     :: residual_enst, res_enst

  REAL(MK)                                     :: c_1_2dx2, c_1_2dy2, c_1_2dz2
  REAL(MK)                                     :: c_1_dt
  REAL(MK)                                     :: c1, c2

  REAL(MK)                           :: dvdx, dwdx, dudy, dwdy, dudz, dvdz

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr        = 0
  npen        = 0
  force       = 0.0_MK
  sum_vel_pen = 0.0_MK
  offset = (/1,1,1/)

!---------------------------------------------------------------------------------!
! Map velocity to penalisation mesh (overwrite)
!---------------------------------------------------------------------------------!
  CALL pmlib_mesh_map( topo_all%cuboid,topo_penal%cuboid,ndim,ierr, &
                     & map_ghost = .TRUE.)
  CALL pmlib_comm_pack(mesh%vel,ierr)
  CALL pmlib_comm_send(ierr)
  CALL pmlib_comm_unpack( topo_penal%cuboid,mesh_penal%vel_pen,0,ierr, &
                        & clear = .TRUE.)
  CALL pmlib_comm_finalise(ierr)

!---------------------------------------------------------------------------------!
! Map ghosts
!---------------------------------------------------------------------------------!
  DO i = 1,ndim
    CALL pmlib_mesh_map_ghost(topo_penal%cuboid,i,ndim,ierr,incl_edges = .FALSE.)
    CALL pmlib_comm_pack(mesh_penal%vel_pen,ierr)
    CALL pmlib_comm_send(ierr)
    CALL pmlib_comm_unpack( topo_penal%cuboid,mesh_penal%vel_pen,0, &
                          & ierr,clear=.FALSE.)
    CALL pmlib_comm_finalise(ierr)
  END DO

!---------------------------------------------------------------------------------!
! Calculate initial penalisation velocity
!---------------------------------------------------------------------------------!
  dx     = topo_penal%cuboid(rank)%dx
  xmin   = topo_penal%cuboid(rank)%xmin
  ncell  = topo_penal%cuboid(rank)%ncell
  nghost = topo_penal%cuboid(rank)%nghost
  DO k = 1-nghost(5),ncell(3)+nghost(6)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        mesh_penal%vel_pen(1,i,j,k) = mesh_penal%solid_vel(1,i,j,k) & 
                                  & - mesh_penal%vel_pen(1,i,j,k)
        mesh_penal%vel_pen(2,i,j,k) = mesh_penal%solid_vel(2,i,j,k) &
                                  & - mesh_penal%vel_pen(2,i,j,k)
        mesh_penal%vel_pen(3,i,j,k) = mesh_penal%solid_vel(3,i,j,k) &
                                  & - mesh_penal%vel_pen(3,i,j,k)

        mesh_penal%vel(1,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * mesh_penal%vel_pen(1,i,j,k)
        mesh_penal%vel(2,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * mesh_penal%vel_pen(2,i,j,k)
        mesh_penal%vel(3,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * mesh_penal%vel_pen(3,i,j,k)

        sum_vel_pen = sum_vel_pen + SQRT( mesh_penal%vel(1,i,j,k)**2 &
                                      & + mesh_penal%vel(2,i,j,k)**2 &
                                      & + mesh_penal%vel(3,i,j,k)**2 )

      END DO
    END DO
  END DO

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients for the finite difference stencils
!---------------------------------------------------------------------------------!
  c_1_2dx2 = 1.0_MK/(2.0_MK*dx(1))
  c_1_2dy2 = 1.0_MK/(2.0_MK*dx(2))
  c_1_2dz2 = 1.0_MK/(2.0_MK*dx(3))

  c_1_dt = 1.0_MK/dtime

  c1 = 1.0_MK/6.0_MK
  c2 = 4.0_MK/3.0_MK

!---------------------------------------------------------------------------------!
! Iterative penalisation
!---------------------------------------------------------------------------------!
  DO WHILE( ( residual_force .GT. penal_tolerance .AND. &  
          & npen .LT. penal_max ) .OR. npen .EQ. 0)
    npen = npen + 1

!---------------------------------------------------------------------------------!
! Calculate the penaliation vorticity
!---------------------------------------------------------------------------------!
    sum_dforce   = 0.0_MK
    res_enst     = 0.0_MK
    sum_enst_pen = 0.0_MK
! 2nd order FD
    IF (pmlib_fd_order .EQ. 2) THEN
      DO k = 1,ncell(3)
        pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
        DO j = 1,ncell(2)
          py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
          DO i = 1,ncell(1)
            px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)
! x-derivatives
            dvdx = ( - 1.0_MK * mesh_penal%vel(2,i-1,j,k) & 
                   & + 1.0_MK * mesh_penal%vel(2,i+1,j,k) ) &
                   & * c_1_2dx2
            dwdx = ( - 1.0_MK * mesh_penal%vel(3,i-1,j,k) & 
                   & + 1.0_MK * mesh_penal%vel(3,i+1,j,k) ) &
                   & * c_1_2dx2
! y-derivatives
            dudy = ( - 1.0_MK * mesh_penal%vel(1,i,j-1,k) &
                   & + 1.0_MK * mesh_penal%vel(1,i,j+1,k) ) &
                   & * c_1_2dy2
            dwdy = ( - 1.0_MK * mesh_penal%vel(3,i,j-1,k) &
                   & + 1.0_MK * mesh_penal%vel(3,i,j+1,k) ) &
                   & * c_1_2dy2
! z-derivatives
            dudz = ( - 1.0_MK * mesh_penal%vel(1,i,j,k-1) &
                   & + 1.0_MK * mesh_penal%vel(1,i,j,k+1) ) &
                   & * c_1_2dz2
            dvdz = ( - 1.0_MK * mesh_penal%vel(2,i,j,k-1) &
                   & + 1.0_MK * mesh_penal%vel(2,i,j,k+1) ) &
                   & * c_1_2dz2

! The penalisation vorticity
            dvort(1) = penal_relaxation * (dwdy - dvdz)
            dvort(2) = penal_relaxation * (dudz - dwdx)
            dvort(3) = penal_relaxation * (dvdx - dudy)

            mesh_penal%vort(1,i,j,k) = mesh_penal%vort(1,i,j,k) + dvort(1)
            mesh_penal%vort(2,i,j,k) = mesh_penal%vort(2,i,j,k) + dvort(2)
            mesh_penal%vort(3,i,j,k) = mesh_penal%vort(3,i,j,k) + dvort(3)

!---------------------------------------------------------------------------------!
! Calculate the added penalisation force
!---------------------------------------------------------------------------------!
            sum_dforce(1) = sum_dforce(1) + (py * dvort(3) - pz * dvort(2))
            sum_dforce(2) = sum_dforce(2) + (pz * dvort(1) - px * dvort(3))
            sum_dforce(3) = sum_dforce(3) + (px * dvort(2) - py * dvort(1))

            res_enst     = res_enst + dvort(1)**2 &
                                  & + dvort(2)**2 &
                                  & + dvort(3)**2
            sum_enst_pen = sum_enst_pen + mesh_penal%vort(1,i,j,k)**2 &
                                      & + mesh_penal%vort(2,i,j,k)**2 & 
                                      & + mesh_penal%vort(3,i,j,k)**2

          END DO
        END DO
      END DO

! 4nd order FD
    ELSE IF (pmlib_fd_order .EQ. 4) THEN
      DO k = 1,ncell(3)
        pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
        DO j = 1,ncell(2)
          py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
          DO i = 1,ncell(1)
            px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)
! x-derivatives
            dvdx = (   c1 * mesh_penal%vel(2,i-2,j,k) &
                   & - c2 * mesh_penal%vel(2,i-1,j,k) & 
                   & + c2 * mesh_penal%vel(2,i+1,j,k) & 
                   & - c1 * mesh_penal%vel(2,i+2,j,k) & 
                   & ) * c_1_2dx2
            dwdx = (   c1 * mesh_penal%vel(3,i-2,j,k) &
                   & - c2 * mesh_penal%vel(3,i-1,j,k) & 
                   & + c2 * mesh_penal%vel(3,i+1,j,k) & 
                   & - c1 * mesh_penal%vel(3,i+2,j,k) & 
                   & ) * c_1_2dx2
! y-derivatives
            dudy = (   c1 * mesh_penal%vel(1,i,j-2,k) &
                   & - c2 * mesh_penal%vel(1,i,j-1,k) & 
                   & + c2 * mesh_penal%vel(1,i,j+1,k) & 
                   & - c1 * mesh_penal%vel(1,i,j+2,k) & 
                   & ) * c_1_2dy2
            dwdy = (   c1 * mesh_penal%vel(3,i,j-2,k) &
                   & - c2 * mesh_penal%vel(3,i,j-1,k) & 
                   & + c2 * mesh_penal%vel(3,i,j+1,k) & 
                   & - c1 * mesh_penal%vel(3,i,j+2,k) & 
                   & ) * c_1_2dy2
! y-derivatives
            dudz = (   c1 * mesh_penal%vel(1,i,j,k-2) &
                   & - c2 * mesh_penal%vel(1,i,j,k-1) & 
                   & + c2 * mesh_penal%vel(1,i,j,k+1) & 
                   & - c1 * mesh_penal%vel(1,i,j,k+2) & 
                   & ) * c_1_2dz2
            dvdz = (   c1 * mesh_penal%vel(2,i,j,k-2) &
                   & - c2 * mesh_penal%vel(2,i,j,k-1) & 
                   & + c2 * mesh_penal%vel(2,i,j,k+1) & 
                   & - c1 * mesh_penal%vel(2,i,j,k+2) & 
                   & ) * c_1_2dz2

! The penalisation vorticity
            dvort(1) = penal_relaxation * (dwdy - dvdz)
            dvort(2) = penal_relaxation * (dudz - dwdx)
            dvort(3) = penal_relaxation * (dvdx - dudy)

            mesh_penal%vort(1,i,j,k) = mesh_penal%vort(1,i,j,k) + dvort(1)
            mesh_penal%vort(2,i,j,k) = mesh_penal%vort(2,i,j,k) + dvort(2)
            mesh_penal%vort(3,i,j,k) = mesh_penal%vort(3,i,j,k) + dvort(3)

!---------------------------------------------------------------------------------!
! Calculate the added penalisation force
!---------------------------------------------------------------------------------!
            sum_dforce(1) = sum_dforce(1) + (py * dvort(3) - pz * dvort(2))
            sum_dforce(2) = sum_dforce(2) + (pz * dvort(1) - px * dvort(3))
            sum_dforce(3) = sum_dforce(3) + (px * dvort(2) - py * dvort(1))

            res_enst     = res_enst + dvort(1)**2 &
                                  & + dvort(2)**2 &
                                  & + dvort(3)**2
            sum_enst_pen = sum_enst_pen + mesh_penal%vort(1,i,j,k)**2 &
                                      & + mesh_penal%vort(2,i,j,k)**2 & 
                                      & + mesh_penal%vort(3,i,j,k)**2
          END DO
        END DO
      END DO
    ENDIF

!---------------------------------------------------------------------------------!
! Solve the Poisson equation
!---------------------------------------------------------------------------------!
    CALL poisson_solver_set_reprojection(2,.FALSE.)
    CALL poisson_solver_push(2,offset,mesh_penal%vort)
    CALL poisson_solver_solve3d(2)
    CALL poisson_solver_pull(2,offset,mesh_penal%vel)


! Map ghost cells
    DO i = 1,ndim
      CALL pmlib_mesh_map_ghost(topo_penal%cuboid,i,3,ierr,incl_edges = .FALSE.)
      CALL pmlib_comm_pack(mesh_penal%vel,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack(topo_penal%cuboid,mesh_penal%vel,0,ierr,clear=.FALSE.)
      CALL pmlib_comm_finalise(ierr)
    END DO

!    CALL pmlib_poisson_solve( patch_penal,topo_penal,mesh_penal,ierr, &
!                              & reg_vort = .FALSE.,  reproj = .FALSE. )
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(rank,caller,'Failed to set up Poisson solver.')
      GOTO 9999
    ENDIF

!---------------------------------------------------------------------------------!
! Calculate penalisation velocity and residual velocity
!---------------------------------------------------------------------------------!
    res_vel = 0.0_MK
    DO k = 1-nghost(5),ncell(3)+nghost(6)
      DO j = 1-nghost(3),ncell(2)+nghost(4)
        DO i = 1-nghost(1),ncell(1)+nghost(2)
          mesh_penal%vel(1,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * (mesh_penal%vel_pen(1,i,j,k) - mesh_penal%vel(1,i,j,k) )
          mesh_penal%vel(2,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * (mesh_penal%vel_pen(2,i,j,k) - mesh_penal%vel(2,i,j,k) )
          mesh_penal%vel(3,i,j,k) = mesh_penal%solid(1,i,j,k) &
                  & * (mesh_penal%vel_pen(3,i,j,k) - mesh_penal%vel(3,i,j,k) )

          res_vel = res_vel + SQRT( mesh_penal%vel(1,i,j,k)**2 & 
                                & + mesh_penal%vel(2,i,j,k)**2 & 
                                & + mesh_penal%vel(3,i,j,k)**2 )
        END DO
      END DO
    END DO

!---------------------------------------------------------------------------------!
! Reduce integrals
!---------------------------------------------------------------------------------!
    CALL MPI_ALLREDUCE( sum_vel_pen,vel_pen,1,mpi_prec_real,MPI_SUM,comm,ierr )
    CALL MPI_ALLREDUCE( sum_enst_pen,enst_pen,1,mpi_prec_real,MPI_SUM,comm,ierr )
    CALL MPI_ALLREDUCE( sum_dforce,dforce,ndim,mpi_prec_real,MPI_SUM,comm,ierr )

    CALL MPI_ALLREDUCE( res_enst,residual_enst,1,mpi_prec_real,MPI_SUM,comm,ierr )
    CALL MPI_ALLREDUCE( res_vel,residual_vel,1,mpi_prec_real,MPI_SUM,comm,ierr )

!---------------------------------------------------------------------------------!
! Calculate total force
!---------------------------------------------------------------------------------!
    dforce = -dforce * dx(1)*dx(2)*dx(3) * c_1_dt
    force  = force + dforce

!---------------------------------------------------------------------------------!
! Calculate residuals
!---------------------------------------------------------------------------------!
    residual_enst  = residual_enst/enst_pen
    residual_vel   = residual_vel/vel_pen
    residual_force = SQRT( dforce(1)**2 + dforce(2)**2 + dforce(3)**2 ) &
                 & / SQRT( force(1)**2 + force(2)**2 + force(3)**2 )

    IF(rank .EQ. 0)THEN
!      WRITE(*,*)'residual: ', residual_vel, residual_force, residual_enst
    END IF

  END DO

!---------------------------------------------------------------------------------!
! Map penalisation vorticity back to the main mesh (addition)
!---------------------------------------------------------------------------------!
  CALL pmlib_mesh_map(topo_penal%cuboid,topo_all%cuboid,nvort,ierr, map_ghost=.TRUE.)
  CALL pmlib_comm_pack(mesh_penal%vort,ierr)
  CALL pmlib_comm_send(ierr)
  CALL pmlib_comm_unpack(topo_all%cuboid,mesh%vort,1,ierr,clear=.FALSE.)
  CALL pmlib_comm_finalise(ierr)

!---------------------------------------------------------------------------------!
! Map ghost cells on the main mesh
!---------------------------------------------------------------------------------!
  DO i = 1,ndim
    CALL pmlib_mesh_map_ghost(topo_all%cuboid,i,nvort,ierr,incl_edges = .FALSE.)
    CALL pmlib_comm_pack(mesh%vort,ierr)
    CALL pmlib_comm_send(ierr)
    CALL pmlib_comm_unpack(topo_all%cuboid,mesh%vort,0,ierr,clear=.FALSE.)
    CALL pmlib_comm_finalise(ierr)
  END DO

!---------------------------------------------------------------------------------!
! Clear penalisation vorticity
!---------------------------------------------------------------------------------!
  mesh_penal%vort = 0.0_MK

!---------------------------------------------------------------------------------!
! Recalculate the velocity field on the full domain
!---------------------------------------------------------------------------------!
  IF(reproject )THEN
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
    CALL pmlib_write(rank,caller,'Failed to solve Poisson equation.')
    GOTO 9999
  ENDIF
! add free stream velocity
  ncell  = topo_all%cuboid(rank)%ncell
  nghost = topo_all%cuboid(rank)%nghost
  DO k = 1-nghost(5),ncell(3)+nghost(6)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        mesh%vel(1,i,j,k) = mesh%vel(1,i,j,k) + vel_inf(1)
        mesh%vel(2,i,j,k) = mesh%vel(2,i,j,k) + vel_inf(2)
        mesh%vel(2,i,j,k) = mesh%vel(2,i,j,k) + vel_inf(3)
      END DO
    END DO
  END DO

!---------------------------------------------------------------------------------!
! Output penalisation data
!---------------------------------------------------------------------------------!
!  IF(rank .EQ. 0)THEN
!    WRITE(*,*) npen, residual_vel, residual_force, residual_enst
!  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE penalisation_solve
