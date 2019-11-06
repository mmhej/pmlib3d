!---------------------------------------------------------------------------------!
! diagnostics_energy_mesh.f
!---------------------------------------------------------------------------------!
! Calculates and outputs the diagnostics of the flowcase
!---------------------------------------------------------------------------------!
SUBROUTINE diagnostics_energy_mesh(topo,mesh,ierr)

USE pmlib_mod_topology
USE pmlib_mod_mesh

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER    :: topo
  TYPE(class_mesh)                             :: mesh
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  REAL(MK)                                     :: PI
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  CHARACTER(LEN=11)                            :: diag_format
  INTEGER,DIMENSION(ndim)                      :: ncell
  REAL(MK),DIMENSION(ndim)                     :: xmin, dx

  REAL(MK),DIMENSION(ndim)                     :: vort_1mom

  REAL(MK)                                     :: enst

  REAL(MK)                                     :: sum_kin, kin_energy
  REAL(MK)                                     :: sum_enst, enstrophy
  REAL(MK)                                     :: sum_hel, helicity
  REAL(MK)                                     :: sum_lamb, lamb

  REAL(MK)                                     :: max_vort

  REAL(MK)                                     :: eff_visc

  REAL(MK)                                     :: c_1_3

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

  sum_kin       = 0.0_MK
  sum_enst      = 0.0_MK
  sum_hel       = 0.0_MK
  sum_lamb      = 0.0_MK

  max_vort      = 0.0_MK

!---------------------------------------------------------------------------------!
! Calculate vorticity moments on the processor
!---------------------------------------------------------------------------------!
  xmin   = topo(rank)%xmin
  ncell  = topo(rank)%ncell
  dx     = topo(rank)%dx

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------!
! Mesh calculations
!---------------------------------------------------------------------------------!
  DO k = 1,ncell(3)
    DO j = 1,ncell(2)
      DO i = 1,ncell(1)

! The first vorticity moment
!        vort_1mom(1) = mesh%pos(2,i) * mesh%vort(3,i) &
!                   & - mesh%pos(3,i) * mesh%vort(2,i) 
!        vort_1mom(2) = mesh%pos(3,i) * mesh%vort(1,i) &
!                   & - mesh%pos(1,i) * mesh%vort(3,i)
!        vort_1mom(3) = mesh%pos(1,i) * mesh%vort(2,i) &
!                   & - mesh%pos(2,i) * mesh%vort(1,i) 

! the local enstrophy
        enst = mesh%vort(1,i,j,k) * mesh%vort(1,i,j,k) &
           & + mesh%vort(2,i,j,k) * mesh%vort(2,i,j,k) &
           & + mesh%vort(3,i,j,k) * mesh%vort(3,i,j,k)

! integrate the kinetic energy
        sum_kin = sum_kin + 0.5_MK * ( mesh%vel(1,i,j,k) * mesh%vel(1,i,j,k) &
                                   & + mesh%vel(2,i,j,k) * mesh%vel(2,i,j,k) &
                                   & + mesh%vel(3,i,j,k) * mesh%vel(3,i,j,k) &
                                   & ) *dx(1)*dx(2)*dx(3)

! integrate the kinetic energy Saffman:1992 eq.(3.11.3)
!        sum_kin = sum_kin + &
!              & ( mesh%vel(1,i) * vort_1mom(1) &
!              & + mesh%vel(2,i) * vort_1mom(2) &
!              & + mesh%vel(3,i) * vort_1mom(3) &
!              & ) *dx(1)*dx(2)*dx(3)

! integrate the enstropy
        sum_enst = sum_enst + enst *dx(1)*dx(2)*dx(3)

! integrate the helicity energy
        sum_hel = sum_hel + &
                & ( mesh%vort(1,i,j,k) * mesh%vel(1,i,j,k) &
                & + mesh%vort(2,i,j,k) * mesh%vel(2,i,j,k) &
                & + mesh%vort(3,i,j,k) * mesh%vel(3,i,j,k) &
                & )**2/enst *dx(1)*dx(2)*dx(3)

! integrate the lamb energy
        sum_lamb = sum_lamb + &
                 & ( ( mesh%vort(2,i,j,k) * mesh%vel(3,i,j,k)      &
                 &   - mesh%vort(3,i,j,k) * mesh%vel(2,i,j,k) )**2 &
                 & + ( mesh%vort(3,i,j,k) * mesh%vel(1,i,j,k)      &
                 &   - mesh%vort(1,i,j,k) * mesh%vel(3,i,j,k) )**2 &
                 & + ( mesh%vort(1,i,j,k) * mesh%vel(2,i,j,k)      &
                 &   - mesh%vort(2,i,j,k) * mesh%vel(1,i,j,k) )**2 &
                 & )/enst *dx(1)*dx(2)*dx(3)

! max vorticity
        IF( SQRT( mesh%vort(1,i,j,k)**2 &
              & + mesh%vort(2,i,j,k)**2 &
              & + mesh%vort(3,i,j,k)**2) .GT. max_vort ) THEN

          max_vort = SQRT( mesh%vort(1,i,j,k)**2 &
                       & + mesh%vort(2,i,j,k)**2 &
                       & + mesh%vort(3,i,j,k)**2 )

        END IF

      END DO
    END DO
  END DO

! Reduce and communicate
  CALL MPI_ALLREDUCE( max_vort,vort_max,1,mpi_prec_real,MPI_MAX,comm,ierr  )
  CALL MPI_ALLREDUCE( sum_enst,enstrophy,1,mpi_prec_real,MPI_SUM,comm,ierr  )

  CALL MPI_REDUCE( sum_kin,kin_energy,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_hel,helicity,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_lamb,lamb,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )

!---------------------------------------------------------------------------------!
! Calculate effective viscosity
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0 .AND. idiag .NE. 0)THEN
    eff_visc = - (kin_energy - diag(2,idiag))/dtime / enstrophy 
    IF(visc .GT. 0.0_MK )THEN
      WRITE(*,'(A,E10.3)') 'Effective viscosity error: ', ABS(eff_visc - visc)/visc
    ELSE
      WRITE(*,'(A,E10.3)') 'Effective viscosity error: ', eff_visc
    END IF
  END IF

!---------------------------------------------------------------------------------!
! Fill in diagnostics array
!---------------------------------------------------------------------------------!
  idiag = idiag + 1
  IF(rank .EQ. 0)THEN
    IF (ioutput_diag .NE. 0)THEN
!---------------------------------------------------------------------------------!
! Set up diagnostics array
!---------------------------------------------------------------------------------!
      ndiag = 5

      IF( .NOT. ALLOCATED(diag) )THEN
        ALLOCATE( diag(ndiag,ioutput_diag) )
      END IF

      diag(1,idiag) = time
      diag(2,idiag) = kin_energy
      diag(3,idiag) = enstrophy
      diag(4,idiag) = helicity
      diag(5,idiag) = lamb

!---------------------------------------------------------------------------------!
! Write diagnostics to file
!---------------------------------------------------------------------------------!
      IF( MOD(itime,ioutput_diag) .EQ. 0 &
        & .OR. abort .OR. itime .EQ. ntime )THEN
        IF( itime .EQ. 0 )THEN
          WRITE(diag_format,'(A,I2,A)') '(',ndiag, 'A)'
          OPEN(20,FILE = './output/diagnostics.dat')
          WRITE(20,TRIM(diag_format)) '#               Time', &
                                    & '      Kinetic energy', &
                                    & '           Enstrophy', &
                                    & '     Helicity energy', &
                                    & '         Lamb energy'
        ELSE
          OPEN(20,FILE = './output/diagnostics.dat',POSITION = 'APPEND')
        END IF

!      WRITE(diag_format,'(A,I2,A)') '(',ndiag, 'E11.3)'
        WRITE(diag_format,'(A,I2,A)') '(',ndiag, 'E20.12)'

        DO i = 1,idiag
          WRITE(20,TRIM(diag_format)) diag(:,i)
        END DO
        CLOSE(20)
        idiag = 0
      END IF
    END IF

  END IF

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE diagnostics_energy_mesh
