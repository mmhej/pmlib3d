!---------------------------------------------------------------------------------!
! diagnostics_helix.f
!---------------------------------------------------------------------------------!
! Calculates and outputs the diagnostics of the flowcase
!---------------------------------------------------------------------------------!
SUBROUTINE diagnostics_helix(topo,part,ierr)

USE pmlib_mod_topology
USE pmlib_mod_particles
!USE pmlib_mod_mesh

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  TYPE(class_topology),DIMENSION(:),POINTER    :: topo
  TYPE(class_particles)                        :: part
  INTEGER, INTENT(OUT)                         :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER                                      :: i,j,k
  REAL(MK)                                     :: px,py,pz

  CHARACTER(LEN=11)                            :: diag_format
  INTEGER,DIMENSION(ndim)                      :: ncell
  REAL(MK),DIMENSION(ndim)                     :: xmin, dx

  REAL(MK),DIMENSION(ndim)                     :: vort_1mom
  REAL(MK),DIMENSION(ndim)                     :: vort_2mom

  REAL(MK),DIMENSION(ndim)                     :: sum_imp, impulse
  REAL(MK),DIMENSION(ndim)                     :: sum_ang_imp, ang_impulse
  REAL(MK)                                     :: sum_kin, kin_energy
  REAL(MK)                                     :: sum_enst, enstrophy
  REAL(MK)                                     :: sum_hel, helicity
  REAL(MK),DIMENSION(ndim)                     :: sum_vort_cen
  REAL(MK)                                     :: sum_vort_rad, vort_radius
  REAL(MK)                                     :: max_vort

  REAL(MK)                                     :: eff_visc

  REAL(MK)                                     :: c_1_3

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

  sum_imp      = 0.0_MK
  sum_ang_imp  = 0.0_MK
  sum_kin      = 0.0_MK
  sum_enst     = 0.0_MK
  sum_hel      = 0.0_MK
  sum_vort_cen = 0.0_MK
  sum_vort_rad = 0.0_MK
  max_vort     = 0.0_MK

!---------------------------------------------------------------------------------!
! Calculate vorticity moments on the processor
!---------------------------------------------------------------------------------!
  xmin   = topo(rank)%xmin
  ncell  = topo(rank)%ncell
  dx     = topo(rank)%dx

!---------------------------------------------------------------------------------!
! Setup irrationel coefficients
!---------------------------------------------------------------------------------!
  c_1_3 = 1.0_MK/3.0_MK

!---------------------------------------------------------------------------------!
! Particle calculations
!---------------------------------------------------------------------------------!
  DO i = 1,topo(rank)%npart

! The first vorticity moment
    vort_1mom(1) = part%pos(2,i) * part%vort(3,i) &
               & - part%pos(3,i) * part%vort(2,i) 
    vort_1mom(2) = part%pos(3,i) * part%vort(1,i) &
               & - part%pos(1,i) * part%vort(3,i)
    vort_1mom(3) = part%pos(1,i) * part%vort(2,i) &
               & - part%pos(2,i) * part%vort(1,i) 

! The second vorticity moment
    vort_2mom(1) = part%pos(2,i) * vort_1mom(3) &
               & - part%pos(3,i) * vort_1mom(2) 
    vort_2mom(2) = part%pos(3,i) * vort_1mom(1) &
               & - part%pos(1,i) * vort_1mom(3)
    vort_2mom(3) = part%pos(1,i) * vort_1mom(2) &
               & - part%pos(2,i) * vort_1mom(1) 

! calculate the impulse Saffman:1992 eq.(3.2.8)
    sum_imp(1) = sum_imp(1) + 0.5_MK*vort_1mom(1)*dx(1)*dx(2)*dx(3)
    sum_imp(2) = sum_imp(2) + 0.5_MK*vort_1mom(2)*dx(1)*dx(2)*dx(3)
    sum_imp(3) = sum_imp(3) + 0.5_MK*vort_1mom(3)*dx(1)*dx(2)*dx(3)

! calculate the angular impulse Saffman:1992 eq.(3.5.5)
    sum_ang_imp(1) = sum_ang_imp(1) + c_1_3 * vort_2mom(1)*dx(1)*dx(2)*dx(3)
    sum_ang_imp(2) = sum_ang_imp(2) + c_1_3 * vort_2mom(2)*dx(1)*dx(2)*dx(3)
    sum_ang_imp(3) = sum_ang_imp(3) + c_1_3 * vort_2mom(3)*dx(1)*dx(2)*dx(3)

! calculate the kinetic energy Saffman:1992 eq.(3.11.3)
    sum_kin = sum_kin + &
          & ( part%vel(1,i) * vort_1mom(1) &
          & + part%vel(2,i) * vort_1mom(2) &
          & + part%vel(3,i) * vort_1mom(3) &
          & ) *dx(1)*dx(2)*dx(3)

! calculate the enstropy
    sum_enst = sum_enst + &
          & ( part%vort(1,i) * part%vort(1,i) &
          & + part%vort(2,i) * part%vort(2,i) &
          & + part%vort(3,i) * part%vort(3,i) &
          & ) *dx(1)*dx(2)*dx(3)

! calculate the helicity Saffman:1992 eq.(3.12.1)
    sum_hel = sum_hel + &
          & ( part%vel(1,i) * part%vort(1,i) &
          & + part%vel(2,i) * part%vort(2,i) &
          & + part%vel(3,i) * part%vort(3,i) &
          & ) *dx(1)*dx(2)*dx(3)

! Calculate the vortex centroid using the enstrophy
    sum_vort_cen(1) = sum_vort_cen(1) + & 
          & ( part%vort(1,i) * part%vort(1,i) &
          & + part%vort(2,i) * part%vort(2,i) &
          & + part%vort(3,i) * part%vort(3,i) &
          & ) * part%pos(1,i)*dx(1)*dx(2)*dx(3)
    sum_vort_cen(2) = sum_vort_cen(2) + & 
          & ( part%vort(1,i) * part%vort(1,i) &
          & + part%vort(2,i) * part%vort(2,i) &
          & + part%vort(3,i) * part%vort(3,i) &
          & ) * part%pos(2,i)*dx(1)*dx(2)*dx(3)
    sum_vort_cen(3) = sum_vort_cen(3) + & 
          & ( part%vort(1,i) * part%vort(1,i) &
          & + part%vort(2,i) * part%vort(2,i) &
          & + part%vort(3,i) * part%vort(3,i) &
          & ) * part%pos(3,i)*dx(1)*dx(2)*dx(3)

! Calculate the vortex ring radius using the enstrophy
    sum_vort_rad = sum_vort_rad + & 
          & ( part%vort(1,i) * part%vort(1,i) &
          & + part%vort(2,i) * part%vort(2,i) &
          & + part%vort(3,i) * part%vort(3,i) &
          & ) * SQRT(part%pos(2,i)**2 &
          & + part%pos(3,i)**2 )*dx(1)*dx(2)*dx(3)

! max vorticity
    IF( SQRT( part%vort(1,i)**2 &
            & + part%vort(2,i)**2 &
            & + part%vort(3,i)**2) .GT. max_vort ) THEN

      max_vort = SQRT( part%vort(1,i)**2 &
                   & + part%vort(2,i)**2 &
                   & + part%vort(3,i)**2 )

    END IF
  END DO

!---------------------------------------------------------------------------------!
! Reduce and communicate
!---------------------------------------------------------------------------------!
  CALL MPI_ALLREDUCE( max_vort,vort_max,1,mpi_prec_real,MPI_MAX,comm,ierr  )
  CALL MPI_ALLREDUCE( sum_vort_cen,vort_centroid,ndim,mpi_prec_real,MPI_SUM,comm,ierr  )
  CALL MPI_ALLREDUCE( sum_enst,enstrophy,1,mpi_prec_real,MPI_SUM,comm,ierr  )

  CALL MPI_REDUCE( sum_imp,impulse,ndim,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_ang_imp,ang_impulse,ndim,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_kin,kin_energy,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_hel,helicity,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )
  CALL MPI_REDUCE( sum_vort_rad,vort_radius,1,mpi_prec_real,MPI_SUM,0,comm,ierr  )

!---------------------------------------------------------------------------------!
! Normalise with enstrophy
!---------------------------------------------------------------------------------!
  DO i = 1,ndim
    vort_centroid(i) = vort_centroid(i)/enstrophy
  END DO
  vort_radius = vort_radius/enstrophy

!---------------------------------------------------------------------------------!
! Calculate effective viscosity
!---------------------------------------------------------------------------------!
  IF(rank .EQ. 0 .AND. idiag .NE. 0)THEN
    eff_visc = - (kin_energy - diag(8,idiag))/dtime / enstrophy 
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
!---------------------------------------------------------------------------------!
! Set up diagnostics array
!---------------------------------------------------------------------------------!
    ndiag = 12

    IF( .NOT. ALLOCATED(diag) )THEN
      ALLOCATE( diag(ndiag,ioutput_diag) )
    END IF

    diag( 1,idiag) = time
    diag( 2,idiag) = impulse(1)
    diag( 3,idiag) = impulse(2)
    diag( 4,idiag) = impulse(3)
    diag( 5,idiag) = ang_impulse(1)
    diag( 6,idiag) = ang_impulse(2)
    diag( 7,idiag) = ang_impulse(3)
    diag( 8,idiag) = kin_energy
    diag( 9,idiag) = enstrophy
    diag(10,idiag) = helicity
    diag(11,idiag) = vort_centroid(1)
    diag(12,idiag) = vort_radius

!---------------------------------------------------------------------------------!
! Write diagnostics to file
!---------------------------------------------------------------------------------!
    IF (ioutput_diag .NE. 0)THEN
      IF( MOD(itime,ioutput_diag) .EQ. 0 &
        & .OR. abort .OR. itime .EQ. ntime )THEN
        IF( itime .EQ. 0 )THEN
          WRITE(diag_format,'(A,I2,A)') '(',ndiag, 'A)'
          OPEN(20,FILE = 'diagnostics.dat')
          WRITE(20,TRIM(diag_format)) '#               Time', &
                                    & '           Impulse x', &
                                    & '           Impulse y', &
                                    & '           Impulse z', &
                                    & '      Ang. impulse x', &
                                    & '      Ang. impulse y', &
                                    & '      Ang. impulse z', &
                                    & '      Kinetic energy', &
                                    & '           Enstrophy', &
                                    & '            Helicity', &
                                    & '            Centroid', &
                                    & '         Ring radius'
        ELSE
          OPEN(20,FILE = 'diagnostics.dat',POSITION = 'APPEND')
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

END SUBROUTINE diagnostics_helix
