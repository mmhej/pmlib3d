# 1 "sources/topology3d.f"
# 1 "<command-line>"
# 1 "sources/topology3d.f"
!---------------------------------------------------------------------------
! fourier3d.f
!---------------------------------------------------------------------------
! __ __
! ___ _ _ / /__ / /
! | _ \/ ` \ _ /

! |_| Mads M. Hejlesen* and Jens H. Walther*^
! * Technical University of Denmark
! ^ ETH-Zurich
! mmhej@mek.dtu.dk
!
!---------------------------------------------------------------------------
PROGRAM fourier3d

USE pmlib_mod_parameters
USE pmlib_mod_write
USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_output

IMPLICIT NONE
!---------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------------
  CHARACTER(LEN=9), PARAMETER :: caller = 'topology3d'
  REAL(MK),PARAMETER :: PI = ACOS(-1.0_MK)

  CHARACTER(LEN=256) :: nodename
  INTEGER :: name_len
  INTEGER :: ierr



  TYPE(class_patch) :: patch
  TYPE(class_topology),DIMENSION(:),POINTER :: topology => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_xpencil => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_ypencil => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_zpencil => NULL()
  TYPE(class_topology),DIMENSION(:),POINTER :: topo_kernel => NULL()


! TYPE(class_patch),DIMENSION(:,:),POINTER :: patch => NULL()
! TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topology => NULL()

! TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_xpencil => NULL()
! TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_ypencil => NULL()
! TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_zpencil => NULL()
! TYPE(class_topology),DIMENSION(:,:,:),POINTER :: topo_kernel => NULL()

  ierr = 0
!---------------------------------------------------------------------------
! Initialise pmlib parameters
!---------------------------------------------------------------------------
  CALL pmlib_parameters_initiate(ierr)

!---------------------------------------------------------------------------
! Input parameters
!---------------------------------------------------------------------------
! nlevels = 1
! ALLOCATE(npatches(nlevels))
! npatches = (/ 1 /)
! ALLOCATE( patch(nlevels,MAXVAL(npatches)) )

!---------------------------------------------------------------------------
! Setup input patch structure
!--------------------------------------------------------------------------
  patch%ptype = 2
  patch%level = 1
  patch%parent = (/ 0, 0 /)
  patch%dx = 1.0_MK/REAL(64,MK)
  patch%xmin = (/ 0.0_MK, 0.0_MK, 0.0_MK /)
  patch%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
  patch%nghost = 4
  patch%bound_cond = (/ 0, 0, 0 /)

! patch(1,1)%ptype = 2
! patch(1,1)%level = 1
! patch(1,1)%parent = (/ 0, 0 /)
! patch(1,1)%dx = 1.0_MK/REAL(64,MK)
! patch(1,1)%xmin = (/ 0.0_MK, 0.0_MK, 0.0_MK /)
! patch(1,1)%xmax = (/ 1.0_MK, 1.0_MK, 1.0_MK /)
! patch(1,1)%nghost = 4
! patch(1,1)%bound_cond = (/ 0, 0, 0 /)

!-------------------------------------------------------------------------
! Start MPI
!-------------------------------------------------------------------------
  comm = MPI_COMM_WORLD
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(comm,rank,ierr)
  CALL MPI_COMM_SIZE(comm,nproc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(nodename, name_len, ierr)
  WRITE(*,*) 'rank: ',rank, ' starting on "',nodename(1:name_len),'"'

!-------------------------------------------------------------------------
! Adjust patches
!-------------------------------------------------------------------------
  CALL pmlib_patch_adjust(patch,ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to adjust patches.')
    GOTO 9999
  ENDIF

!-------------------------------------------------------------------------
! Topology create
!-------------------------------------------------------------------------
  CALL pmlib_topology_cuboid(patch,topology,ierr)

!----------------------------------------------------------------------------
! Set up pencil topologies
!----------------------------------------------------------------------------
  CALL pmlib_topology_pencil(patch,topo_xpencil,1,(/ 0 ,0, 0 /),ierr)
  CALL pmlib_topology_pencil(patch,topo_ypencil,2,(/ 0 ,0, 0 /),ierr)
  CALL pmlib_topology_pencil(patch,topo_zpencil,3,(/ 0 ,0, 0 /),ierr)

!----------------------------------------------------------------------------
! Set up kernel topology
!----------------------------------------------------------------------------
  CALL pmlib_topology_kernel(patch,topo_kernel,1,(/ 1 ,1, 1 /),ierr)

!----------------------------------------------------------------------------
! Plot
!----------------------------------------------------------------------------
  CALL pmlib_output_topology(patch,topology,'topology',ierr)


!---------------------------------------------------------------------------
! Finalise MPI
!---------------------------------------------------------------------------
  CALL MPI_finalize(ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(rank,caller,'Failed to finalize MPI.')
  END IF
  GOTO 1111

!---------------------------------------------------------------------------
! Return·
!---------------------------------------------------------------------------
 9999 CONTINUE
  CALL MPI_ABORT(comm,ierr)

 1111 CONTINUE
END PROGRAM fourier3d
