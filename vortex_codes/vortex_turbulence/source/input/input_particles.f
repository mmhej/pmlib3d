!---------------------------------------------------------------------------------!
! input_particles.f
!---------------------------------------------------------------------------------!
SUBROUTINE input_particles(infile,patch,topo_all,mesh,part,trunc,ierr)

USE pmlib_mod_patch
USE pmlib_mod_topology
USE pmlib_mod_mesh
USE pmlib_mod_particles
USE pmlib_mod_repatch

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! arguments
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN)                 :: infile
  TYPE(class_patch)                            :: patch
  TYPE(class_topology_all)                     :: topo_all
  TYPE(class_mesh)                             :: mesh
  TYPE(class_particles)                        :: part
  REAL(MK)                                     :: trunc
  INTEGER, INTENT(INOUT)                       :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=22), PARAMETER          :: caller = 'input_particles'
  INTEGER, PARAMETER                           :: IO_unit_part  = 10
  INTEGER                                      :: ltag
  INTEGER                                      :: i,j,k
  INTEGER                                      :: npart_global, npart
  INTEGER                                      :: nhead
  INTEGER, DIMENSION(nproc)                    :: npart_local
  REAL(MK),DIMENSION(ndim)                     :: dx
  REAL(MK)                                     :: dummy
  INTEGER                                      :: stat(MPI_STATUS_SIZE)

  REAL(MK)                                     :: mag_vort
  REAL(MK),DIMENSION(:),POINTER                :: input 
  REAL(MK),DIMENSION(:,:),POINTER              :: buffer

!---------------------------------------------------------------------------------!
! initiate subroutine
!---------------------------------------------------------------------------------!
  ierr  = 0
  nhead = 0

!---------------------------------------------------------------------------------!
! Deallocate the particles if they are already allocated
!---------------------------------------------------------------------------------!
  IF(ASSOCIATED( part%pos ))THEN
    DEALLOCATE(part%pos,stat=ierr)
  END IF

  IF(ASSOCIATED( part%vel ))THEN
    DEALLOCATE(part%vel,stat=ierr)
  END IF

  IF(ASSOCIATED( part%vort ))THEN
    DEALLOCATE(part%vort,stat=ierr)
  END IF

  IF(ASSOCIATED( part%dvort ))THEN
    DEALLOCATE(part%dvort,stat=ierr)
  END IF

  IF(ASSOCIATED( part%vel_rk ))THEN
    DEALLOCATE(part%vel_rk,stat=ierr)
  END IF

  IF(ASSOCIATED( part%dvort_rk ))THEN
    DEALLOCATE(part%dvort_rk,stat=ierr)
  END IF


  IF(rank .EQ. 0)THEN
!---------------------------------------------------------------------------------!
! Read restart parameters of the simulation 
!---------------------------------------------------------------------------------!
    OPEN( IO_unit_part,FILE=TRIM(infile), & 
        & FORM='UNFORMATTED',ACTION='READ')

    READ(IO_unit_part)itime
    nhead = nhead + 1

write(*,*)itime

    READ(IO_unit_part)time
    nhead = nhead + 1
    READ(IO_unit_part)npart_global
    nhead = nhead + 1

write(*,*)itime,time,npart_global

!---------------------------------------------------------------------------------!
! First read all and truncate particles
!---------------------------------------------------------------------------------!
    ALLOCATE( input(nvort + ndim) )
    j = 0
    DO i = 1,npart_global
      READ(IO_unit_part) input
      IF( SQRT(input(4)**2 + input(5)**2 + input(6)**2) .GE. trunc)THEN
        j = j + 1
      END IF
    END DO
    npart_global = j

    IF(npart_global .EQ. 0)THEN
      ierr = - 1
      CALL pmlib_write(rank,caller, &
                    & 'No particles above truncation value.')
      GO TO 9999
    END IF

    CLOSE(IO_unit_part)

!---------------------------------------------------------------------------------!
! Devide particles amongst the processors 
!---------------------------------------------------------------------------------!
    DO i = 1,nproc
      npart_local(i) = INT( REAL(npart_global,MK)/REAL(nproc,MK) )
      IF( i .LE. MOD(npart_global,npart_local(i)) )THEN
        npart_local(i) = npart_local(i) + 1
      END IF

    END DO

    WRITE(*,*)'Global number of particles:',npart_global

  END IF

!---------------------------------------------------------------------------------!
! Communicate number of particles
!---------------------------------------------------------------------------------!
  CALL MPI_SCATTER( npart_local, 1, mpi_prec_int, &
                  & npart, 1, mpi_prec_int, 0, comm,ierr)

!---------------------------------------------------------------------------------!
! Communicate simulation parameters
!---------------------------------------------------------------------------------!
  CALL MPI_BCAST(resolution, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(vort_max, 1, mpi_prec_real, 0, comm, ierr)

  CALL MPI_BCAST(itime, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(time, 1, mpi_prec_real, 0, comm, ierr)

  CALL MPI_BCAST(visc, 1, mpi_prec_real, 0, comm, ierr)

  CALL MPI_BCAST(reproject, 1, mpi_prec_log, 0, comm, ierr)
  CALL MPI_BCAST(regularise_vort, 1, mpi_prec_log, 0, comm, ierr)

  CALL MPI_BCAST(iremesh, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(irepatch, 1, mpi_prec_int, 0, comm, ierr)

  CALL MPI_BCAST(remesh_trunc, 1, mpi_prec_real, 0, comm, ierr)
  CALL MPI_BCAST(repatch_trunc, 1, mpi_prec_real, 0, comm, ierr)

  CALL MPI_BCAST(pmlib_poisson_order, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(pmlib_time_integration_order, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(pmlib_fd_order, 1, mpi_prec_int, 0, comm, ierr)
  CALL MPI_BCAST(pmlib_interpolation_order, 1, mpi_prec_int, 0, comm, ierr)

!---------------------------------------------------------------------------------!
! Allocate particles
!---------------------------------------------------------------------------------!
  topo_all%cuboid(rank)%npart = npart
  ALLOCATE( part%pos(ndim,npart), &
          & part%vort(nvort,npart), &
          & part%vel(ndim,npart), &
          & part%dvort(nvort,npart) )

  ALLOCATE( part%vel_rk(ndim,npart), &
          & part%dvort_rk(nvort,npart) )

!---------------------------------------------------------------------------------!
! Rewind and read particles
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0)THEN
    OPEN( IO_unit_part,FILE=TRIM(infile), & 
        & FORM='UNFORMATTED',ACTION='READ')

! Ignore header entries
    DO i = 1,nhead
      READ(IO_unit_part)
    END DO

! Particles on the root
    vort_max = 0.0_MK
    j = 0
    DO WHILE (j .LT. npart)
      READ(IO_unit_part) input

      mag_vort = SQRT(input(4)**2 + input(5)**2 + input(6)**2)
      IF( mag_vort .GE. trunc)THEN
        j = j + 1
        DO k = 1,ndim
          part%pos(k,j) = REAL(input(k),MK)
        END DO
        DO k = 1,nvort
          part%vort(k,j) = REAL(input(k+ndim),MK)
        END DO

        IF( mag_vort .GT. vort_max)THEN
          vort_max = mag_vort
        END IF

      END IF
    END DO
  END IF

! Send particles to the other processors
  DO i = 1,nproc-1

    IF(rank .EQ. 0)THEN
      ALLOCATE( buffer(nvort+ndim,npart_local(i+1)) )

      j = 0
      DO WHILE (j .LT. npart_local(i+1))
        READ(IO_unit_part) input

        mag_vort = SQRT(input(4)**2 + input(5)**2 + input(6)**2)
        IF( mag_vort .GE. trunc)THEN
          j = j + 1
          DO k = 1,nvort+ndim
            buffer(k,j) = REAL( input(k) ,MK)
          END DO

          IF( mag_vort .GT. vort_max)THEN
            vort_max = mag_vort
          END IF
        END IF
      END DO

      CALL MPI_SEND(buffer, (nvort+ndim)*npart_local(i+1), &
                     & mpi_prec_real,i,i,mpi_comm,stat,ierr)      

      DEALLOCATE( buffer )
    ELSEIF(rank .EQ. i)THEN

      ALLOCATE( buffer(nvort+ndim,npart) )

! recieve particles from root
      CALL MPI_RECV(buffer, (nvort+ndim)*npart, &
                     & mpi_prec_real,0,i,mpi_comm,stat,ierr)

! unpack buffer
      DO j = 1,npart
        DO k = 1,ndim
          part%pos(k,j) = buffer(k,j)
        END DO
        DO k = 1,nvort
          part%vort(k,j) = buffer(k+ndim,j)
        END DO  
      END DO

      DEALLOCATE( buffer )
    END IF
  END DO

  IF(rank .EQ. 0)THEN
    CLOSE(IO_unit_part)
  END IF

!---------------------------------------------------------------------------------!
! Broadcast maximum vorticity
!---------------------------------------------------------------------------------!
  CALL MPI_BCAST( vort_max, 1, mpi_prec_real, 0, mpi_comm, ierr)

!---------------------------------------------------------------------------------!
! Re-construct patches and topology
!---------------------------------------------------------------------------------!
  patch%dx = 1.0_MK/REAL(resolution,MK)
  CALL pmlib_repatch(patch,topo_all,part,mesh,repatch_trunc,10,ierr)

!---------------------------------------------------------------------------------!
! Output to terminal
!---------------------------------------------------------------------------------!
  CALL pmlib_write( rank,'pmlib_output_particle', &
                    & 'Particles loaded from file.',verb=.TRUE.)

!---------------------------------------------------------------------------------!
! De-allocate local pointers
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0)THEN
    DEALLOCATE( input )
  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE input_particles
