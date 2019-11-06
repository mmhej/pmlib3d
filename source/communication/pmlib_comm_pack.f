!---------------------------------------------------------------------------------!
! pmlib_comm_pack.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
! Packs the communication buffer
!---------------------------------------------------------------------------------!


!=================================================================================!
! MESH OF TYPE LOGICAL
!=================================================================================!
SUBROUTINE pmlib_comm_pack_mesh_logical(mesh,ierr)

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  LOGICAL,DIMENSION(:,:,:), POINTER             :: mesh
  INTEGER, INTENT(OUT)                          :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=26)                     :: caller = 'pmlib_comm_pack_mesh_real'
  INTEGER                               :: iproc, jproc, ivar
  INTEGER                               :: i,j,k,l,m
  INTEGER, DIMENSION(pmlib_ndim)        :: imin_send, imax_send
  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_send
  INTEGER                               :: nvar

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  nvar = 1
  ivar = communication%npacked
  communication%npacked = communication%npacked + nvar

!---------------------------------------------------------------------------------!
! Pack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. iproc)THEN

      ncell_send(1) = communication%ncell_send(1,m)  
      ncell_send(2) = communication%ncell_send(2,m)
      ncell_send(3) = communication%ncell_send(3,m)

      imin_send(1)  = communication%imin_send(1,m)  
      imin_send(2)  = communication%imin_send(2,m)
      imin_send(3)  = communication%imin_send(3,m)

! If send and recieve is the same then store directly in recieve buffer
      IF(jproc .EQ. iproc)THEN

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)

              IF( mesh( imin_send(1) + j - 1, imin_send(2) + k - 1, &
                      & imin_send(3) + l - 1 ))THEN
                comm_buffer(m)%mesh_recieve(1+ivar,j,k,l) = 1.0_MK
              ELSE
                comm_buffer(m)%mesh_recieve(1+ivar,j,k,l) = 0.0_MK
              END IF

            END DO
          END DO
        END DO

      ELSE

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)

              IF( mesh( imin_send(1) + j - 1, imin_send(2) + k - 1, &
                      & imin_send(3) + l - 1 ))THEN
                comm_buffer(m)%mesh_send(1+ivar,j,k,l) = 1.0_MK
              ELSE
                comm_buffer(m)%mesh_send(1+ivar,j,k,l) = 0.0_MK
              END IF

            END DO
          END DO
        END DO

      END IF

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_pack_mesh_logical




!=================================================================================!
! MESH OF TYPE REAL
!=================================================================================!
SUBROUTINE pmlib_comm_pack_mesh_real(mesh,ierr)

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  REAL(MK),DIMENSION(:,:,:,:), POINTER          :: mesh
  INTEGER, INTENT(OUT)                          :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=26)                     :: caller = 'pmlib_comm_pack_mesh_real'
  INTEGER                               :: iproc, jproc, ivar
  INTEGER                               :: i,j,k,l,m
  INTEGER, DIMENSION(pmlib_ndim)        :: imin_send, imax_send
  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_send
  INTEGER                               :: nvar

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  nvar = SIZE(mesh,1)
  ivar = communication%npacked
  communication%npacked = communication%npacked + nvar

!---------------------------------------------------------------------------------!
! Pack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. iproc)THEN

      ncell_send(1) = communication%ncell_send(1,m)  
      ncell_send(2) = communication%ncell_send(2,m)  
      ncell_send(3) = communication%ncell_send(3,m)  

      imin_send(1)  = communication%imin_send(1,m)  
      imin_send(2)  = communication%imin_send(2,m)  
      imin_send(3)  = communication%imin_send(3,m) 

! If send and recieve is the same then store directly in recieve buffer
      IF(jproc .EQ. iproc)THEN

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)
              DO i = 1,nvar
                comm_buffer(m)%mesh_recieve(i+ivar,j,k,l) = &
                            & mesh(i,imin_send(1) + j - 1, & 
                            &        imin_send(2) + k - 1, &
                            &        imin_send(3) + l - 1)
              END DO
            END DO
          END DO
        END DO

      ELSE

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)
              DO i = 1,nvar
                comm_buffer(m)%mesh_send(i+ivar,j,k,l) = & 
                            & mesh(i,imin_send(1) + j - 1, & 
                            &        imin_send(2) + k - 1, &
                            &        imin_send(3) + l - 1)
              END DO
            END DO
          END DO
        END DO
      END IF

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_pack_mesh_real




!=================================================================================!
! MESH OF TYPE COMPLEX
!=================================================================================!
SUBROUTINE pmlib_comm_pack_mesh_complex(mesh,ierr)
!---------------------------------------------------------------------------------!
! Packs the communication buffer
!---------------------------------------------------------------------------------!

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  COMPLEX(MKC),DIMENSION(:,:,:,:), POINTER      :: mesh
  INTEGER, INTENT(OUT)                          :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=29)           :: caller = 'pmlib_comm_pack_mesh_complex'
  INTEGER                               :: iproc, jproc, ivar
  INTEGER                               :: i,j,k,l,m
  INTEGER, DIMENSION(pmlib_ndim)        :: imin_send, imax_send
  INTEGER, DIMENSION(pmlib_ndim)        :: ncell_send
  INTEGER                               :: nvar

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  nvar = SIZE(mesh,1)
  ivar = communication%npacked
  communication%npacked = communication%npacked + 2*nvar

!---------------------------------------------------------------------------------!
! Pack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. iproc)THEN

      ncell_send(1) = communication%ncell_send(1,m)  
      ncell_send(2) = communication%ncell_send(2,m)  
      ncell_send(3) = communication%ncell_send(3,m)  

      imin_send(1)  = communication%imin_send(1,m)  
      imin_send(2)  = communication%imin_send(2,m)  
      imin_send(3)  = communication%imin_send(3,m) 

! If send and recieve is the same then store directly in recieve buffer
      IF(jproc .EQ. iproc)THEN

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)
              DO i = 1,nvar
                comm_buffer(m)%mesh_recieve(2*i-1+ivar,j,k,l) = &
                            & REAL(mesh(i,imin_send(1) + j - 1, & 
                            &             imin_send(2) + k - 1, &
                            &             imin_send(3) + l - 1) )

                comm_buffer(m)%mesh_recieve(2*i+ivar,j,k,l) = &
                            & AIMAG(mesh(i,imin_send(1) + j - 1, & 
                            &             imin_send(2) + k - 1, &
                            &             imin_send(3) + l - 1) )
              END DO
            END DO
          END DO
        END DO

      ELSE

        DO l = 1,ncell_send(3)
          DO k = 1,ncell_send(2)
            DO j = 1,ncell_send(1)
              DO i = 1,nvar
                comm_buffer(m)%mesh_send(2*i-1+ivar,j,k,l) = & 
                            & REAL(mesh(i,imin_send(1) + j - 1, & 
                            &             imin_send(2) + k - 1, &
                            &             imin_send(3) + l - 1) )
                comm_buffer(m)%mesh_send(2*i+ivar,j,k,l) = & 
                            & AIMAG(mesh(i,imin_send(1) + j - 1, & 
                            &             imin_send(2) + k - 1, &
                            &             imin_send(3) + l - 1) )
              END DO
            END DO
          END DO
        END DO
      END IF

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_pack_mesh_complex




!=================================================================================!
! PARTICLES OF TYPE REAL
!=================================================================================!
SUBROUTINE pmlib_comm_pack_part(part,ierr)
!---------------------------------------------------------------------------------!
! Packs the communication buffer
!---------------------------------------------------------------------------------!

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  REAL(MK),DIMENSION(:,:), POINTER              :: part
  INTEGER, INTENT(OUT)                          :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=21)                     :: caller = 'pmlib_comm_pack_part'
  INTEGER                               :: iproc, jproc, ivar
  INTEGER                               :: i,j,k,l,m
  INTEGER                               :: npart
  INTEGER                               :: nvar, npart_tot
  INTEGER, DIMENSION(:), POINTER        :: ipart

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Find the packing index
!---------------------------------------------------------------------------------!
  nvar      = SIZE(part,1)
  npart_tot = SIZE(part,2)
  ivar = communication%npacked
  communication%npacked = communication%npacked + nvar

!---------------------------------------------------------------------------------!
! Pack the communication buffer
!---------------------------------------------------------------------------------!
  DO m = 1,communication%ncomm

    iproc = communication%proc_send(m)
    jproc = communication%proc_recieve(m)

    IF(mpi_rank .EQ. iproc)THEN

      npart = communication%npart_send(m) 

      ALLOCATE(ipart(npart))

      IF(communication%comm_dir(m) .EQ. 0)THEN
        ipart = communication%ipart0
      ELSEIF(communication%comm_dir(m) .EQ. 1)THEN
        ipart = communication%ipart1
      ELSEIF(communication%comm_dir(m) .EQ. 2)THEN
        ipart = communication%ipart2
      END IF

      IF(jproc .EQ. iproc)THEN
        DO j = 1,npart
          DO i = 1,nvar
            comm_buffer(m)%part_recieve(i+ivar,j) = part(i,ipart(j))
          END DO
        END DO
      ELSE
        DO j = 1,npart
          DO i = 1,nvar
            comm_buffer(m)%part_send(i+ivar,j) = part(i,ipart(j))
          END DO
        END DO
      END IF

!---------------------------------------------------------------------------------!
! Deallocate particle list
!---------------------------------------------------------------------------------!
      DEALLOCATE(ipart,stat=ierr)
    END IF
  END DO

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN
END SUBROUTINE pmlib_comm_pack_part

