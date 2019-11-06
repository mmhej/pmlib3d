!---------------------------------------------------------------------------------!
! pmlib_visualise_postscript.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_visualise_postscript(pix_rgb,tag,ierr)

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! arguments
!---------------------------------------------------------------------------------!
  REAL(MK), DIMENSION(:,:,:),POINTER              :: pix_rgb
  CHARACTER(LEN=*), INTENT(IN)                    :: tag
  INTEGER, INTENT(OUT)                            :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER               :: IO_unit_ps = 10
  INTEGER, PARAMETER               :: QK = SELECTED_REAL_KIND(18)
  CHARACTER(LEN=256)               :: outfile
  INTEGER                          :: ltag
  INTEGER                          :: i,j,k,l
  INTEGER,DIMENSION(4)             :: b
  INTEGER,DIMENSION(5)             :: c
  REAL(QK)                         :: base10
  CHARACTER(LEN=5)                 :: asc
  INTEGER                          :: lasc, lline

  INTEGER,DIMENSION(:),POINTER     :: rgb
  INTEGER, DIMENSION(2)            :: npixel

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0

!---------------------------------------------------------------------------------!
! find number of pixels
!---------------------------------------------------------------------------------!
  npixel(1) = SIZE(pix_rgb,2)
  npixel(2) = SIZE(pix_rgb,3)

!---------------------------------------------------------------------------------!
! Re-size to 1D array
!---------------------------------------------------------------------------------!
  ALLOCATE( rgb(3*npixel(1)*npixel(2)) )
  l = 0
  DO k = 1,npixel(2)
    DO j = 1,npixel(1)
      DO i = 1,3
        l = l + 1
        rgb(l) = FLOOR(pix_rgb(i,j,k)*pix_rgb(4,j,k))
      END DO
    END DO
  END DO

!---------------------------------------------------------------------------------!
! Write window to eps file
!---------------------------------------------------------------------------------!
  ltag = LEN_TRIM(tag)
  WRITE(outfile,'(A,A)') tag(1:ltag),'.ps'

!---------------------------------------------------------------------------------!
! Write postscript file
!---------------------------------------------------------------------------------!
  OPEN(IO_unit_ps,FILE=outfile,FORM='FORMATTED')

  WRITE(IO_unit_ps,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
  WRITE(IO_unit_ps,'(A,2I6)') '%%BoundingBox: 0 0 ',npixel(1),npixel(2)
  WRITE(IO_unit_ps,'(A)')''
  WRITE(IO_unit_ps,'(A)')'/DeviceRGB setcolorspace'
  WRITE(IO_unit_ps,'(A)')'gsave'
  WRITE(IO_unit_ps,'(A)')''
  WRITE(IO_unit_ps,'(A)') '%%%%BeginImage'
  WRITE(IO_unit_ps,'(A)') '0 0 translate'
  WRITE(IO_unit_ps,'(2I6,A)') npixel(1),npixel(2),' scale'
  WRITE(IO_unit_ps,'(A)') '<<'
  WRITE(IO_unit_ps,'(A)') '  /ImageType 1'
  WRITE(IO_unit_ps,'(A,I6)') '  /Width', npixel(1)
  WRITE(IO_unit_ps,'(A,I6)') '  /Height', npixel(2)
  WRITE(IO_unit_ps,'(A)') '  /BitsPerComponent 8'
  WRITE(IO_unit_ps,'(A,2(I6,A))') '  /ImageMatrix [ ', npixel(1), ' 0 0 ', npixel(2), ' 0 0 ]'
  WRITE(IO_unit_ps,'(A)') '  /Decode [ 0 1 0 1 0 1 ]'
  WRITE(IO_unit_ps,'(A)') '  /DataSource currentfile /ASCII85Decode filter'
  WRITE(IO_unit_ps,'(A)') '  /MultipleDataSources false'
  WRITE(IO_unit_ps,'(A)') '  /Interpolate true'
  WRITE(IO_unit_ps,'(A)') '>>'
  WRITE(IO_unit_ps,'(A)') 'image'

!---------------------------------------------------------------------------------!
! Encode to ASCII85
!---------------------------------------------------------------------------------!
  lline = 0
  b = 0
  c = 0
  DO i = 1,CEILING(REAL(3*npixel(1)*npixel(2),MK)/4.0_MK)

    DO j = 1,4
      IF( 4*i-4+j .LE. 3*npixel(1)*npixel(2) )THEN
        b(5-j) = rgb(4*i-4+j)
      ELSE
        b(5-j) = 0
      END IF
    END DO

    base10 = REAL(b(1),QK)
    base10 = base10 + REAL(b(2),QK)*256.0_QK
    base10 = base10 + REAL(b(3),QK)*(256.0_QK**2)
    base10 = base10 + REAL(b(4),QK)*(256.0_QK**3)

    DO j = 1,5
      c(6-j) = INT(MAX(0.0_QK, base10 - 85.0_QK*AINT(base10/85.0_QK) ))
      base10 = AINT(base10/85.0_QK)
    END DO

    IF( c(1) .EQ. 0 .AND. c(2) .EQ. 0 .AND. &
      & c(3) .EQ. 0 .AND. c(4) .EQ. 0 .AND. c(5) .EQ. 0)THEN
      WRITE(asc,'(A)') 'z'
      lasc = 1
    ELSE
      WRITE(asc,'(5A)') CHAR(c(1)+33),CHAR(c(2)+33), &
                      & CHAR(c(3)+33),CHAR(c(4)+33),CHAR(c(5)+33)
      lasc = 5
    END IF

    IF(lline .GT. 77)THEN
      WRITE(IO_unit_ps,'(A)',ADVANCE='YES') asc(1:lasc)
      lline = 0
    ELSE
      WRITE(IO_unit_ps,'(A)',ADVANCE='NO') asc(1:lasc)
      lline = lline + lasc
    END IF

  ENDDO

  WRITE(IO_unit_ps,'(A)',ADVANCE='YES')''
  WRITE(IO_unit_ps,'(A)')'~>'
  WRITE(IO_unit_ps,'(A)') '%%%%EndImage'
  WRITE(IO_unit_ps,'(A)') 'grestore' 
  WRITE(IO_unit_ps,'(A)') 'showpage'
  WRITE(IO_unit_ps,'(A)')'%%EOF'

  CLOSE(IO_unit_ps)


 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_visualise_postscript
