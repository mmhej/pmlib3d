!---------------------------------------------------------------------------------!
! input_setup.f
!---------------------------------------------------------------------------------!
SUBROUTINE input_setup(infile,patch,ierr)

USE pmlib_mod_patch

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN)                     :: infile
  TYPE(class_patch)                                :: patch
  INTEGER, INTENT(OUT)                             :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=18), PARAMETER             :: caller = 'input_setup'
  INTEGER                                          :: i,j,idx

  INTEGER                                          :: iline, len_line
  CHARACTER(LEN=256)                               :: msg, line, arg, val
  REAL(MK)                                         :: re

!---------------------------------------------------------------------------------!
! Initialise subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! Set up patch defaults ()
!---------------------------------------------------------------------------------!
  patch%level      = 1
  patch%parent     = (/ 0, 0 /)
  patch%nghost     = 4

!---------------------------------------------------------------------------------!
!  Open the file
!---------------------------------------------------------------------------------!
  OPEN(10, FILE=infile, IOSTAT=ierr, ACTION='READ')
  IF (ierr .NE. 0) THEN
    WRITE(msg,'(2A)')'Failed to open setup file: ',TRIM(infile)
    CALL pmlib_write(rank,caller,TRIM(msg))
    GOTO 9999
  ENDIF

!---------------------------------------------------------------------------------!
!  Scan file
!---------------------------------------------------------------------------------!
  iline    = 0
  DO
!---------------------------------------------------------------------------------!
!  Increment line counter
!---------------------------------------------------------------------------------!
    iline = iline + 1

!---------------------------------------------------------------------------------!
!  Read line
!---------------------------------------------------------------------------------!
    READ(10,'(A)',END=100,ERR=200) line
    len_line = LEN_TRIM(line)

!---------------------------------------------------------------------------------!
!  Skip comment or empty lines·
!---------------------------------------------------------------------------------!
    IF (len_line .GT. 0 .AND. line(1:1) .NE. '#') THEN
!---------------------------------------------------------------------------------!
!  Remove spaces in line
!---------------------------------------------------------------------------------!
      j = 0
      DO i = 1,len_line
        IF (line(i:i) .NE. ' ') THEN
          j = j + 1
          line(j:j) = line(i:i)
        ENDIF
      ENDDO

!---------------------------------------------------------------------------------!
!  Update length of string
!---------------------------------------------------------------------------------!
      len_line = j

!---------------------------------------------------------------------------------!
!  Find position of =
!---------------------------------------------------------------------------------!
      idx = INDEX(line,'=')
      IF (idx .LE. 0) GOTO 200

!---------------------------------------------------------------------------------!
!  Get argument and value
!---------------------------------------------------------------------------------!
      arg   = ADJUSTL(line(1:idx-1))
      val   = ADJUSTL(line(idx+1:len_line))

!---------------------------------------------------------------------------------!
!  Convert to upper case
!---------------------------------------------------------------------------------!
      CALL UpperCase(arg,ierr)

!---------------------------------------------------------------------------------!
! Update the appropriate input
!---------------------------------------------------------------------------------!
      IF (TRIM(arg) .EQ. 'FLOWCASE') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) flowcase

      ELSE IF (arg .EQ. 'RESTART') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) restart
      ELSE IF (arg .EQ. 'RESTART_FILE') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) restart_file

      ELSE IF (TRIM(arg) .EQ. 'FREE_STREAM_VELOCITY') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) vel_inf
      ELSE IF (TRIM(arg) .EQ. 'REYNOLDS_NUMBER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) re
        IF(re .EQ. 0.0)THEN
          visc = 0.0_MK
        ELSE
          visc = 1.0_MK/re
        END IF

      ELSE IF (arg .EQ. 'TIME_STEP_SIZE') THEN
        READ(val,*,IOSTAT=ierr,ERR=200)dtime

      ELSE IF (arg .EQ. 'TIME_STEP_PARAMETER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200)ctime

      ELSE IF (arg .EQ. 'NUMBER_TIME_STEPS') THEN
        READ(val,*,IOSTAT=ierr,ERR=200)ntime

      ELSE IF (TRIM(arg) .EQ. 'DOMAIN_MIN') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) patch%xmin
      ELSE IF (TRIM(arg) .EQ. 'DOMAIN_MAX') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) patch%xmax
      ELSE IF (TRIM(arg) .EQ. 'RESOLUTION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) resolution
        patch%dx = 1.0_MK/REAL(resolution,MK)

      ELSE IF (arg .EQ. 'DOMAIN_BOUNDARY_CONDITIONS') THEN
        READ(val,*,IOSTAT=ierr,ERR=200)patch%bound_cond

      ELSE IF (arg .EQ. 'RE-PATCHING_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) irepatch
      IF( irepatch .EQ. 0 )THEN
        patch%ptype = 1
      ELSE
        patch%ptype = 2
      END IF

      ELSE IF (arg .EQ. 'RE-PATCHING_TRUNCATION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) repatch_trunc

      ELSE IF (arg .EQ. 'RE-MESHING_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) iremesh
      ELSE IF (arg .EQ. 'RE-MESHING_TRUNCATION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) remesh_trunc

      ELSE IF (arg .EQ. 'REPROJECTION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) reproject

      ELSE IF (arg .EQ. 'REGULARISATION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) regularise_vort
      ELSE IF (arg .EQ. 'REGULARISATION_RADIUS') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_regularisation_radius

      ELSE IF (arg .EQ. 'PENALISATION_TOLERANCE') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) penal_tolerance
      ELSE IF (arg .EQ. 'PENALISATION_RELAXATION') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) penal_relaxation
      ELSE IF (arg .EQ. 'PENALISATION_MAX_ITERATIONS') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) penal_max

      ELSE IF (arg .EQ. 'POISSON_SOLVER_ORDER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_poisson_order
      ELSE IF (arg .EQ. 'POISSON_SOLVER_KERNEL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_poisson_kernel

      ELSE IF (arg .EQ. 'TIME_INTEGRATION_ORDER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_time_integration_order

      ELSE IF (arg .EQ. 'FINITE_DIFFERENCE_ORDER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_fd_order

      ELSE IF (arg .EQ. 'INTERPOLATION_ORDER') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) pmlib_interpolation_order

      ELSE IF (arg .EQ. 'OUTPUT_DIAGNOSTICS_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) ioutput_diag

      ELSE IF (arg .EQ. 'OUTPUT_PARTICLES_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) ioutput_part

      ELSE IF (arg .EQ. 'PLOT_MESH_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) iplot_mesh

      ELSE IF (arg .EQ. 'PLOT_PARTICLE_INTERVAL') THEN
        READ(val,*,IOSTAT=ierr,ERR=200) iplot_part

      ELSE
        WRITE(msg,'(A,I5,2A)')'Ignoring input line ',iline, &
            & ' : ', TRIM(arg)
        CALL pmlib_write(rank,caller,TRIM(msg))
      END IF
    END IF
  END DO

!---------------------------------------------------------------------------------!
! Reading error
!---------------------------------------------------------------------------------!
 200 CONTINUE
  ierr = -1
  WRITE(msg,'(A,I5,2A)') 'Error reading line: ',iline,     &
       & ' of file: ',TRIM(infile)
  CALL pmlib_write(rank,caller,TRIM(msg))
  GOTO 9999

!---------------------------------------------------------------------------------!
! End of file is reached
!---------------------------------------------------------------------------------!
 100 CONTINUE

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE

RETURN
END SUBROUTINE input_setup




SUBROUTINE UpperCase(string,ierr)

IMPLICIT NONE
!---------------------------------------------------------------------------------!
!  Arguments     
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(INOUT) :: string
  INTEGER         , INTENT(  OUT) :: ierr
!---------------------------------------------------------------------------------!
!  Local variables 
!---------------------------------------------------------------------------------!
  INTEGER          :: i,j
  INTEGER          :: i1,i2,i3,iadd
      
!---------------------------------------------------------------------------------!
!  Initialise subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
!  Convert to upper case
!---------------------------------------------------------------------------------!
  i1   = IACHAR('a') - 1
  i2   = IACHAR('z') + 1
  i3   = IACHAR('A')
  iadd = i3 - i1 - 1
  DO i=1,LEN_TRIM(string)
    j = IACHAR(string(i:i))
    IF (j .GT. i1 .AND. j .LT. i2) THEN
      string(i:i) = CHAR(j+iadd)
    ENDIF
  ENDDO

!---------------------------------------------------------------------------------!
!  Return 
!---------------------------------------------------------------------------------!
RETURN
END SUBROUTINE UpperCase


