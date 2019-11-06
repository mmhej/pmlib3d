!---------------------------------------------------------------------------
! map_mesh3d.f
!---------------------------------------------------------------------------
!                  __    __
!    ___  _ _     / /__ / /
!   | _ \/ ` \ _ / // // _ \
!   | ,_/_,_,_\  \_/\_/\___/
!   |_|                      Mads M. Hejlesen* and Jens H. Walther*^
!                            * Technical University of Denmark
!                            ^ ETH-Zurich
!                            mmhej@mek.dtu.dk
!
!---------------------------------------------------------------------------
PROGRAM map_mesh3d

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patches
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_communication
USE pmlib_mod_output


IMPLICIT NONE
!---------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------------
  CHARACTER(LEN=10), PARAMETER                  :: caller = 'map_mesh3d'
  INTEGER                                       :: ierr

  INTEGER                                       :: ilevel,ipatch
  INTEGER                                       :: i,j,k
  REAL(MK)                                      :: px,py,pz

  REAL(MK)                                      :: r,r0,c

  CHARACTER(LEN=256)                            :: nodename
  INTEGER                                       :: name_len

  REAL(MK),DIMENSION(ndim)                      :: xmin, dx
  INTEGER,DIMENSION(ndim)                       :: ncell
  INTEGER,DIMENSION(2*ndim)                     :: nghost

  CHARACTER(LEN=256)                            :: pltfile

  INTEGER, PARAMETER                            :: nval = 3
  TYPE(class_mesh),DIMENSION(:,:),POINTER       :: mesh => NULL()

  TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_xpencil => NULL()
  TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_ypencil => NULL()
  TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_zpencil => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo_xpen => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo_ypen => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo_zpen => NULL()

  ierr = 0
!---------------------------------------------------------------------------
! Initialise pmlib parameters
!---------------------------------------------------------------------------
  CALL pmlib_parameters_initiate(ierr)


!-------------------------------------------------------------------------
! Start MPI
!-------------------------------------------------------------------------
  comm = MPI_COMM_WORLD
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(comm,rank,ierr)
  CALL MPI_COMM_SIZE(comm,nproc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(nodename, name_len, ierr)
  WRITE(*,*) 'rank: ',rank, ' starting on "',nodename(1:name_len),'"'

!---------------------------------------------------------------------------
! Input parameters
!---------------------------------------------------------------------------
  poisson_order          = 0

  nlevels  = 1
  ALLOCATE(npatches(nlevels))
  npatches = (/ 1 /)

  ALLOCATE( patch(nlevels,MAXVAL(npatches)) )

!---------------------------------------------------------------------------
! Setup input patch structure
!--------------------------------------------------------------------------
  patch(1,1)%ptype      = 2
  patch(1,1)%level      = 1
  patch(1,1)%parent     = (/ 0, 0 /)

  patch(1,1)%xmin       = (/ -1.0_MK, -1.0_MK, -1.0_MK /)
  patch(1,1)%xmax       = (/  1.0_MK,  1.0_MK,  1.0_MK /)
  patch(1,1)%nghost     = 4
  patch(1,1)%bound_cond = (/ 0, 0, 0 /)

  patch(1,1)%dx = 1.0_MK/REAL(32,MK)

!-------------------------------------------------------------------------
! Adjust patches
!-------------------------------------------------------------------------
  CALL pmlib_patches_adjust(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(0,rank,caller,'Failed to adjust patches.')
    GOTO 9999
  ENDIF

!-------------------------------------------------------------------------
! Topology create
!-------------------------------------------------------------------------
  CALL pmlib_topology_cuboid(topology,ierr)

  CALL pmlib_topology_pencil(topo_xpencil,1,(/ 0 , 0, 0 /),ierr)
  CALL pmlib_topology_pencil(topo_ypencil,2,(/ 1 , 1, 1 /),ierr)
  CALL pmlib_topology_pencil(topo_zpencil,3,(/ 0 , 0, 0 /),ierr)

  IF (ierr .NE. 0) THEN
    CALL pmlib_write(0,rank,caller,'Failed to create topology.')
    GOTO 9999
  ENDIF

!----------------------------------------------------------------------------
! Allocate field set
!----------------------------------------------------------------------------
  CALL pmlib_mesh_allocate(mesh,nval,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(0,rank,caller,'Failed to allocate mesh.')
    GOTO 9999
  ENDIF

!----------------------------------------------------------------------------
! Setup initial vorticity field
!----------------------------------------------------------------------------
  xmin   = topology(rank,1,1)%xmin
  dx     = topology(rank,1,1)%dx
  ncell  = topology(rank,1,1)%ncell
  nghost = topology(rank,1,1)%nghost

  DO k = 1-nghost(5),ncell(3)+nghost(6)
    pz = xmin(3) + (REAL(k - 1,MK) + 0.5_MK) *dx(3)
    DO j = 1-nghost(3),ncell(2)+nghost(4)
      py = xmin(2) + (REAL(j - 1,MK) + 0.5_MK) *dx(2)
      DO i = 1-nghost(1),ncell(1)+nghost(2)
        px = xmin(1) + (REAL(i - 1,MK) + 0.5_MK) *dx(1)
        
        mesh(1,1)%val(1,i,j,k) = REAL(rank,MK) + 1.0_MK
        mesh(1,1)%val(2,i,j,k) = 0.0_MK 
        mesh(1,1)%val(3,i,j,k) = 0.0_MK

      END DO !i
    END DO !j
  END DO !k

!----------------------------------------------------------------------------
! Setup topologies
!----------------------------------------------------------------------------
  ilevel = 1
  ipatch = 1

  ALLOCATE( topo(0:nproc-1), topo_xpen(0:nproc-1), topo_ypen(0:nproc-1), topo_zpen(0:nproc-1) )
  DO i = 0,nproc-1
    topo(i)      = topology(i,ilevel,ipatch)
    topo_xpen(i) = topo_xpencil(i,ilevel,ipatch)
    topo_ypen(i) = topo_ypencil(i,ilevel,ipatch)
    topo_zpen(i) = topo_zpencil(i,ilevel,ipatch)
  END DO

!----------------------------------------------------------------------------
! Map mesh interior
!----------------------------------------------------------------------------
!  CALL MPI_BARRIER(comm,ierr)

  IF(.TRUE.)THEN
    CALL pmlib_mesh_map(topo,topo_ypen,nval,ierr)
    CALL pmlib_comm_pack(mesh(ilevel,ipatch)%val,ierr)
    CALL pmlib_comm_send(ierr)
    CALL pmlib_comm_unpack(topo_ypen,mesh(ilevel,ipatch)%val,ierr)
    CALL pmlib_comm_finalise(ierr)

    WRITE(pltfile,'(A)') 'out'
    CALL pmlib_output_vtk_vec_mag(topo_ypencil(rank,ilevel,ipatch),mesh(ilevel,ipatch)%val,pltfile,ierr,ghost = .TRUE.)

  END IF
!----------------------------------------------------------------------------
! Map boundaries for interpolation
!----------------------------------------------------------------------------
  IF(.FALSE.)THEN
    DO i = 1,ndim
      CALL pmlib_mesh_map_interpolation(topo,i,nval,ierr)
      CALL pmlib_comm_pack(mesh(ilevel,ipatch)%val,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack(topo,mesh(ilevel,ipatch)%val,ierr)
      CALL pmlib_comm_finalise(ierr)
    END DO

    WRITE(pltfile,'(A)') 'out'
    CALL pmlib_output_vtk_vec_mag(topology(rank,ilevel,ipatch),mesh(ilevel,ipatch)%val,pltfile,ierr,ghost = .TRUE.)

  END IF

!----------------------------------------------------------------------------
! Map boundaries for ghost cells
!----------------------------------------------------------------------------
  IF(.FALSE.)THEN
    DO i = 1,ndim
      CALL pmlib_mesh_map_ghost(topo,i,nval,ierr)
      CALL pmlib_comm_pack(mesh(ilevel,ipatch)%val,ierr)
      CALL pmlib_comm_send(ierr)
      CALL pmlib_comm_unpack(topo,mesh(ilevel,ipatch)%val,ierr)
      CALL pmlib_comm_finalise(ierr)
    END DO

    WRITE(pltfile,'(A)') 'out'
    CALL pmlib_output_vtk_vec_mag(topology(rank,ilevel,ipatch),mesh(ilevel,ipatch)%val,pltfile,ierr,ghost = .TRUE.)
  END IF

!---------------------------------------------------------------------------
! Finalise MPI
!---------------------------------------------------------------------------
  CALL MPI_finalize(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(0,rank,caller,'Failed to finalize MPI.')
  END IF
  GOTO 1111

!---------------------------------------------------------------------------
! Return·
!---------------------------------------------------------------------------
 9999 CONTINUE
  CALL MPI_ABORT(comm,ierr)

 1111 CONTINUE
END PROGRAM map_mesh3d
