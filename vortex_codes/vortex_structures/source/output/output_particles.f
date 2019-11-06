!---------------------------------------------------------------------------------!
! output_particles.f
!---------------------------------------------------------------------------------!
SUBROUTINE output_particles(tag,topo,part,ierr)

USE pmlib_mod_particles
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! arguments
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN)                  :: tag
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo
  TYPE(class_particles)                         :: part
  INTEGER, INTENT(INOUT)                        :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER                         :: prec_output       = MK
  INTEGER, PARAMETER                         :: IO_unit_part      = 10
  CHARACTER(LEN=256)                         :: outfile
  INTEGER                                    :: ltag
  INTEGER                                    :: i,j,k
  INTEGER                                    :: npart
  REAL(MK),DIMENSION(ndim)                   :: dx
  REAL(prec_output),DIMENSION(:),POINTER     :: output 

!---------------------------------------------------------------------------------!
! initiate subroutine
!---------------------------------------------------------------------------------!
  ierr   = 0
  WRITE(outfile,'(2A)') TRIM(tag),'.dat'
  ALLOCATE( output(nvort+ndim) )

!---------------------------------------------------------------------------------!
! output setup 
!---------------------------------------------------------------------------------!
  CALL MPI_REDUCE( topo(rank)%npart,npart,1, &
                 & mpi_prec_int,MPI_SUM,0,comm,ierr )

!---------------------------------------------------------------------------------!
! Write particle file
!---------------------------------------------------------------------------------!
  DO i = 0,nproc-1
    IF(rank .EQ. i)THEN

      IF(rank .EQ. 0)THEN
        OPEN(IO_unit_part,FILE=TRIM(outfile),FORM='UNFORMATTED')

! Write simulation parameters
        WRITE(IO_unit_part)itime + 1
        WRITE(IO_unit_part)time + dtime
        WRITE(IO_unit_part)npart

      ELSE
        OPEN(IO_unit_part,FILE=TRIM(outfile),ACCESS = 'APPEND',FORM='UNFORMATTED')
      END IF

      DO j=1,topo(rank)%npart
        DO k = 1,ndim
          output(k) = REAL(part%pos(k,j),prec_output)
        END DO
        DO k = 1,nvort
          output(k+ndim) = REAL(part%vort(k,j),prec_output)
        END DO

        WRITE(IO_unit_part) output
      END DO

      CLOSE(IO_unit_part)

    END IF

    CALL MPI_BARRIER(comm,ierr)

  END DO

  CALL pmlib_write(rank, 'pmlib_output_particle', &
                  & 'Saved particles.', verb=.TRUE.)

!---------------------------------------------------------------------------------!
! De-allocate local pointers
!---------------------------------------------------------------------------------!
  DEALLOCATE( output )

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE output_particles
