!---------------------------------------------------------------------------------!
! pmlib_comm_unpack.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
! Unpacks the communication buffer
!---------------------------------------------------------------------------------!




!=================================================================================!
! MESH OF TYPE LOGICAL
!=================================================================================!
SUBROUTINE pmlib_comm_unpack_mesh_logical(topo,mesh,operation,ierr,clear)

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN) :: topo
  LOGICAL,DIMENSION(:,:,:), POINTER                     :: mesh
  INTEGER, INTENT(IN)                                   :: operation
  INTEGER, INTENT(OUT)                                  :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                         :: clear

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=31)                  :: caller = 'pmlib_comm_unpack_mesh_logical'
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: iproc, jproc, ivar
  INTEGER,DIMENSION(pmlib_ndim)         :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)       :: nghost
  INTEGER, DIMENSION(pmlib_ndim)        :: imin

  INTEGER                               :: nvar
!  REAL(MK),DIMENSION(:),POINTER         :: recieved

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  nvar = 1
  ivar = communication%npacked - nvar
  communication%npacked = communication%npacked - nvar

!---------------------------------------------------------------------------------!
! Re-allocate the mesh variable
!---------------------------------------------------------------------------------!
  IF( PRESENT(clear) .AND. clear )THEN
    ncell  = topo(mpi_rank)%ncell
    nghost = topo(mpi_rank)%nghost

    IF( ASSOCIATED(mesh) ) THEN
      DEALLOCATE(mesh,stat=ierr)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF

    ALLOCATE(mesh( &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
      GOTO 9999
    ENDIF

    IF( operation .EQ. 0 .OR. operation .EQ. 1 )THEN
      mesh = .FALSE.
    ELSE IF( operation .EQ. -1 )THEN
      mesh = .TRUE.
    END IF

  END IF

!---------------------------------------------------------------------------------!
! Unpack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN

      ncell(1) = communication%ncell_recieve(1,m)  
      ncell(2) = communication%ncell_recieve(2,m) 
      ncell(3) = communication%ncell_recieve(3,m) 

      imin(1)  = communication%imin_recieve(1,m)  
      imin(2)  = communication%imin_recieve(2,m)
      imin(3)  = communication%imin_recieve(3,m)

      DO l = 1,ncell(3)
        DO k = 1,ncell(2)
          DO j = 1,ncell(1)

            IF( operation .EQ. 0 .OR. operation .EQ. 1 )THEN
! replace
              IF(NINT(comm_buffer(m)%mesh_recieve(1+ivar,j,k,l)) .EQ. 1)THEN
                mesh(imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = .TRUE.
              ELSE IF(NINT(comm_buffer(m)%mesh_recieve(1+ivar,j,k,l)) .EQ. 0)THEN
                mesh(imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = .FALSE.
              ELSE
                ierr = -1
                CALL pmlib_write( mpi_rank,caller, &
                                & 'Failed to read logical mesh value')
                GOTO 9999
              END IF

            ELSEIF( operation .EQ. -1 )THEN
! replace with inverse
              IF(NINT(comm_buffer(m)%mesh_recieve(1+ivar,j,k,l)) .EQ. 1)THEN
                mesh(imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = .FALSE.
              ELSE IF(NINT(comm_buffer(m)%mesh_recieve(1+ivar,j,k,l)) .EQ. 0)THEN
                mesh(imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = .TRUE.
              ELSE
                ierr = -1
                CALL pmlib_write( mpi_rank,caller, &
                                & 'Failed to read logical mesh value')
                GOTO 9999
              END IF

            ELSE
              ierr = -1
              CALL pmlib_write(mpi_rank,caller,'Operation unknown.')
              GOTO 9999
            END IF

          END DO
        END DO
      END DO

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_unpack_mesh_logical






!=================================================================================!
! MESH OF TYPE REAL
!=================================================================================!
SUBROUTINE pmlib_comm_unpack_mesh_real(topo,mesh,operation,ierr,clear,ndim)

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER, INTENT(IN) :: topo
  REAL(MK),DIMENSION(:,:,:,:), POINTER                  :: mesh
  INTEGER, INTENT(IN)                                   :: operation
  INTEGER, INTENT(OUT)                                  :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                         :: clear
  INTEGER, INTENT(IN), OPTIONAL                         :: ndim

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=28)           :: caller = 'pmlib_comm_unpack_mesh_real'
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: iproc, jproc, ivar
  INTEGER,DIMENSION(pmlib_ndim)         :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)       :: nghost
  INTEGER,DIMENSION(pmlib_ndim)         :: imin

  INTEGER                               :: nvar
!  REAL(MK),DIMENSION(:),POINTER         :: recieved

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh) )THEN
    nvar = SIZE(mesh,1)
  ELSEIF( PRESENT(ndim) )THEN
    nvar = ndim
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'nvar is not specified.')
    GOTO 9999
  END IF
  ivar = communication%npacked - nvar
  communication%npacked = communication%npacked - nvar

!---------------------------------------------------------------------------------!
! Re-allocate the mesh variable
!---------------------------------------------------------------------------------!
  IF( PRESENT(clear) .AND. clear )THEN
    ncell  = topo(mpi_rank)%ncell
    nghost = topo(mpi_rank)%nghost

    IF( ASSOCIATED(mesh) ) THEN
      DEALLOCATE(mesh,stat=ierr)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF

    ALLOCATE(mesh(nvar, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
      GOTO 9999
    ENDIF

    mesh = 0.0_MK

  END IF

!---------------------------------------------------------------------------------!
! Unpack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN

      ncell(1)   = communication%ncell_recieve(1,m)  
      ncell(2)   = communication%ncell_recieve(2,m)  
      ncell(3)   = communication%ncell_recieve(3,m)  

      imin(1)    = communication%imin_recieve(1,m)  
      imin(2)    = communication%imin_recieve(2,m)  
      imin(3)    = communication%imin_recieve(3,m) 

      DO l = 1,ncell(3)
        DO k = 1,ncell(2)
          DO j = 1,ncell(1)
            DO i = 1,nvar

              IF( operation .EQ. 0 )THEN
! replace
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = &
                  & comm_buffer(m)%mesh_recieve(i+ivar,j,k,l)

              ELSEIF( operation .EQ. 1 )THEN
! add
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = &
                  & mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1)  &
                  & + comm_buffer(m)%mesh_recieve(i+ivar,j,k,l)
              ELSEIF( operation .EQ. -1 )THEN
! subtract
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) = &
                  & mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1)  &
                  & - comm_buffer(m)%mesh_recieve(i+ivar,j,k,l)
              ELSE
                CALL pmlib_write(mpi_rank,caller,'Operation unknown.')
                GOTO 9999
              END IF

            END DO
          END DO
        END DO
      END DO

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_unpack_mesh_real




!=================================================================================!
! MESH OF TYPE COMPLEX
!=================================================================================!
SUBROUTINE pmlib_comm_unpack_mesh_complex(topo,mesh,operation,ierr,clear,ndim)
!---------------------------------------------------------------------------------!
! Unpacks the communication buffer
!---------------------------------------------------------------------------------!

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER,INTENT(IN)    :: topo
  COMPLEX(MKC),DIMENSION(:,:,:,:), POINTER,INTENT(INOUT)  :: mesh
  INTEGER, INTENT(IN)                                     :: operation
  INTEGER, INTENT(OUT)                                    :: ierr
  LOGICAL, INTENT(IN), OPTIONAL                           :: clear
  INTEGER, INTENT(IN), OPTIONAL                           :: ndim

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=31)                   :: caller = 'pmlib_comm_unpack_mesh_complex'
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: iproc, jproc, ivar
  INTEGER,DIMENSION(pmlib_ndim)         :: ncell
  INTEGER,DIMENSION(2*pmlib_ndim)       :: nghost
  INTEGER, DIMENSION(pmlib_ndim)        :: imin
  INTEGER                               :: nvar

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(mesh) )THEN
    nvar = SIZE(mesh,1)
  ELSEIF( PRESENT(ndim) )THEN
    nvar = ndim
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'nvar is not specified.')
    GOTO 9999
  END IF
  ivar = communication%npacked - 2*nvar
  communication%npacked = communication%npacked - 2*nvar

!---------------------------------------------------------------------------------!
! Re-allocate the mesh variable
!---------------------------------------------------------------------------------!
  IF( PRESENT(clear) .AND. clear )THEN
    ncell  = topo(mpi_rank)%ncell
    nghost = topo(mpi_rank)%nghost

    IF( ASSOCIATED(mesh) ) THEN
      DEALLOCATE(mesh,stat=ierr)
      IF (ierr .NE. 0) THEN
        CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
        GOTO 9999
      ENDIF
    END IF

    ALLOCATE(mesh(nvar, &
          & 1-nghost(1):ncell(1)+nghost(2), &
          & 1-nghost(3):ncell(2)+nghost(4), &
          & 1-nghost(5):ncell(3)+nghost(6)), &
          & stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
      GOTO 9999
    ENDIF

    mesh = CMPLX(0.0_MK,0.0_MK,MKC)

  END IF

!---------------------------------------------------------------------------------!
! Unpack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN

      ncell(1)    = communication%ncell_recieve(1,m)  
      ncell(2)    = communication%ncell_recieve(2,m)  
      ncell(3)    = communication%ncell_recieve(3,m)  

      imin(1)    = communication%imin_recieve(1,m)  
      imin(2)    = communication%imin_recieve(2,m)  
      imin(3)    = communication%imin_recieve(3,m) 

      DO l = 1,ncell(3)
        DO k = 1,ncell(2)
          DO j = 1,ncell(1)
            DO i = 1,nvar

              IF( operation .EQ. 0 )THEN
! replace
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) =       &
                  & CMPLX( comm_buffer(m)%mesh_recieve(2*i-1+ivar,j,k,l),         &
                  &        comm_buffer(m)%mesh_recieve(2*i+ivar,j,k,l) ,MKC)

              ELSEIF( operation .EQ. 1 )THEN
! add
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) =       &
                  & mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1)     &
                  & + CMPLX( comm_buffer(m)%mesh_recieve(2*i-1+ivar,j,k,l),       &
                  &          comm_buffer(m)%mesh_recieve(2*i+ivar,j,k,l) ,MKC)
              ELSEIF( operation .EQ. -1 )THEN
! subtract
                mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1) =       &
                  & mesh(i,imin(1) + j - 1, imin(2) + k - 1, imin(3) + l - 1)     &
                  & - CMPLX( comm_buffer(m)%mesh_recieve(2*i-1+ivar,j,k,l),       &
                  &          comm_buffer(m)%mesh_recieve(2*i+ivar,j,k,l) ,MKC)
              ELSE
                CALL pmlib_write(mpi_rank,caller,'Operation unknown.')
                GOTO 9999
              END IF

            END DO
          END DO
        END DO
      END DO

    END IF
  END DO


!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_unpack_mesh_complex




!=================================================================================!
! PARTICLES OF TYPE REAL
!=================================================================================!
SUBROUTINE pmlib_comm_unpack_part(topo,part,ierr,ndim)
!---------------------------------------------------------------------------------!
! Unpacks the communication buffer
!---------------------------------------------------------------------------------!

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER,INTENT(INOUT) :: topo
  REAL(MK),DIMENSION(:,:),POINTER,INTENT(OUT)             :: part
  INTEGER, INTENT(OUT)                                    :: ierr
  INTEGER, INTENT(IN), OPTIONAL                           :: ndim

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=23)                     :: caller = 'pmlib_comm_unpack_part'
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: iproc, jproc, ivar, ipart
  INTEGER                               :: npart

  INTEGER                               :: nvar
  INTEGER                               :: npart_tot

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(part) )THEN
    nvar = SIZE(part,1)
  ELSEIF( PRESENT(ndim) )THEN
    nvar = ndim
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'nvar is not specified.')
    GOTO 9999
  END IF
  ivar = communication%npacked - nvar
  communication%npacked = communication%npacked - nvar

  IF( ASSOCIATED(part) ) THEN
    DEALLOCATE(part,stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
      GOTO 9999
    ENDIF
  END IF

!---------------------------------------------------------------------------------!
! Find size of the output array and re-allocate particle
!---------------------------------------------------------------------------------!
  npart_tot = 0
  DO m = 1,communication%ncomm
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN
      npart_tot = npart_tot + communication%npart_recieve(m)
    END IF

  END DO

  topo(mpi_rank)%npart = npart_tot

  ALLOCATE(part(nvar,npart_tot), stat=ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
    GOTO 9999
  ENDIF
  part = 0.0_MK

!---------------------------------------------------------------------------------!
! Unpack the communication buffer
!---------------------------------------------------------------------------------!
  ipart = 0
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN

      npart = communication%npart_recieve(m)  

      DO j = 1,npart
        ipart = ipart + 1
        DO i = 1,nvar
          part(i,ipart) = comm_buffer(m)%part_recieve(i+ivar,j)
        END DO
      END DO

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_unpack_part




!=================================================================================!
! PARTICLES OF TYPE REAL (CHECKS IF MAPPING IS DONE)
!=================================================================================!
SUBROUTINE pmlib_comm_unpack_part_pos_check(topo,part,map_done,dir,ierr,ndim)
!---------------------------------------------------------------------------------!
! Unpacks the communication buffer
!---------------------------------------------------------------------------------!

USE pmlib_mod_topology

IMPLICIT NONE

!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER,INTENT(INOUT) :: topo
  REAL(MK),DIMENSION(:,:),POINTER,INTENT(OUT)             :: part
  LOGICAL,INTENT(OUT)                                     :: map_done
  INTEGER, INTENT(IN)                                     :: dir 
  INTEGER, INTENT(OUT)                                    :: ierr
  INTEGER, INTENT(IN), OPTIONAL                           :: ndim

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=33)           :: caller = 'pmlib_comm_unpack_part_pos_check'
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: iproc, jproc, ivar, ipart
  INTEGER                               :: npart
  REAL(MK), DIMENSION(pmlib_ndim)       :: xmin, xmax

  INTEGER                               :: nvar
  INTEGER                               :: npart_tot

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr     = 0
  map_done = .TRUE.

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  IF( ASSOCIATED(part) )THEN
    nvar = SIZE(part,1)
  ELSEIF( PRESENT(ndim) )THEN
    nvar = ndim
  ELSE
    ierr = -1
    CALL pmlib_write(mpi_rank,caller,'nvar is not specified.')
    GOTO 9999
  END IF
  ivar = communication%npacked - nvar
  communication%npacked = communication%npacked - nvar

  IF( ASSOCIATED(part) ) THEN
    DEALLOCATE(part,stat=ierr)
    IF (ierr .NE. 0) THEN
      CALL pmlib_write(mpi_rank,caller,'Failed to deallocate array.')
      GOTO 9999
    ENDIF
  END IF

!---------------------------------------------------------------------------------!
! Find size of the output array and re-allocate particle
!---------------------------------------------------------------------------------!
  npart_tot = 0
  DO m = 1,communication%ncomm
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN
      npart_tot = npart_tot + communication%npart_recieve(m)
    END IF

  END DO

  topo(mpi_rank)%npart = npart_tot
  xmin = topo(mpi_rank)%xmin
  xmax = topo(mpi_rank)%xmax 

  ALLOCATE(part(nvar,npart_tot), stat=ierr)
  IF (ierr .NE. 0) THEN
    CALL pmlib_write(mpi_rank,caller,'Failed to allocate array.')
    GOTO 9999
  ENDIF
  part = 0.0_MK

!---------------------------------------------------------------------------------!
! Unpack the communication buffer
!---------------------------------------------------------------------------------!
  ipart = 0
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. jproc)THEN

      npart = communication%npart_recieve(m)  

      DO j = 1,npart
        ipart = ipart + 1
        DO i = 1,nvar
          part(i,ipart) = comm_buffer(m)%part_recieve(i+ivar,j)
        END DO

!---------------------------------------------------------------------------------!
! Check if mapping is done
!---------------------------------------------------------------------------------!
        IF( part(dir,ipart) .LT. xmin(dir) .AND. &
          & topo(mpi_rank)%bound_cond(2*dir-1) .NE. 0 )THEN
          map_done = .FALSE.
        ELSEIF( part(dir,ipart) .GE. xmax(dir) .AND. &
              & topo(mpi_rank)%bound_cond(2*dir) .NE. 0 )THEN 
          map_done = .FALSE.
        END IF
      END DO

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_unpack_part_pos_check
