!---------------------------------------------------------------------------------!
! output_mesh.f
!---------------------------------------------------------------------------------!
! <<HEADER>>
!---------------------------------------------------------------------------------!
SUBROUTINE output_mesh(tag,patch,topo,mesh,ierr)

USE pmlib_mod_mesh
USE pmlib_mod_patch
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN)                  :: tag
  TYPE(class_patch)                             :: patch
  TYPE(class_topology),DIMENSION(:),POINTER     :: topo
  TYPE(class_mesh)                              :: mesh
  INTEGER, INTENT(INOUT)                        :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER, PARAMETER                            :: prec_output = KIND(1E0)
  CHARACTER(LEN=7), PARAMETER                   :: float_type  = 'Float32'
!  INTEGER, PARAMETER                            :: prec_output = KIND(1D0)
!  CHARACTER(LEN=7), PARAMETER                   :: float_type  = 'Float64'

  INTEGER                                       :: i,j,k
  CHARACTER(LEN=256)                            :: outfile
  CHARACTER(LEN=256)                            :: asciiline

  REAL(MK),DIMENSION(pmlib_ndim)                :: xmin, dx
  INTEGER,DIMENSION(pmlib_ndim)                 :: ncell, icell

  INTEGER                                       :: nbytes, ndata

!---------------------------------------------------------------------------------!
! initiate subroutine
!---------------------------------------------------------------------------------!
  ierr   = 0

!---------------------------------------------------------------------------------!
! Write the parallel .pvti file
!---------------------------------------------------------------------------------!
    IF (mpi_rank .EQ. 0) THEN
      WRITE(outfile,'(2A)') TRIM(tag),'.pvti'

      OPEN(17,FILE = TRIM(outfile))
      WRITE(17,'(A)') '<?xml version="1.0"?>'
      WRITE(17,'(2A)') '<VTKFile type="PImageData" version="0.1" ', &
                     & 'byte_order="LittleEndian">'
      WRITE(17,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)')'  <PImageData WholeExtent="', & 
        & 0,patch%ncell(1),0,patch%ncell(2),0,patch%ncell(3), &
        & '" Ghostlevel="0" Origin="',patch%xmin(1),&
        & patch%xmin(2),patch%xmin(3), &
        & '" Spacing="',patch%dx(1),patch%dx(2),patch%dx(3),'">'

      WRITE(17,'(A)') '    <PCellData Vectors="output">'
!---------------------------------------------------------------------------------!
! List field variables
!---------------------------------------------------------------------------------!
      WRITE(17,'(3A)') '      <PDataArray type="',float_type, &
       & '" Name="vort" NumberOfComponents="3" format="appended" offset="0"/>'
      WRITE(17,'(3A)') '      <PDataArray type="',float_type, &
       & '" Name="vel" NumberOfComponents="3" format="appended" offset="0"/>'

!---------------------------------------------------------------------------------!

      WRITE(17,'(A)') '    </PCellData>'

      DO i = 0,mpi_nproc-1
        ncell = topo(i)%ncell
        icell = topo(i)%icell
        WRITE(17,'(A,6(I8),3A,I5.5,A)')'    <Piece  Extent="', & 
        & icell(1) - 1, icell(1) + ncell(1) - 1, &
        & icell(2) - 1, icell(2) + ncell(2) - 1, &
        & icell(3) - 1, icell(3) + ncell(3) - 1, &
        & '" Source="',TRIM(tag),'_P',i,'.vti"/>'
      END DO

      WRITE(17,'(A)') '  </PImageData>'
      WRITE(17,'(A)') '</VTKFile>'
      CLOSE(17)
    END IF

!---------------------------------------------------------------------------------!
! Write individual .vti files
!---------------------------------------------------------------------------------!
    dx     = topo(mpi_rank)%dx
    xmin   = topo(mpi_rank)%xmin
    ncell  = topo(mpi_rank)%ncell
    icell  = topo(mpi_rank)%icell
    ndata  = ncell(1)*ncell(2)*ncell(3)
    nbytes = 0

    WRITE(outfile,'(2A,I5.5,A)') TRIM(tag),'_P',mpi_rank,'.vti'
    OPEN(10,FILE=outfile,FORM='UNFORMATTED',ACCESS='STREAM')

    WRITE(10)'<?xml version="1.0"?>'
    WRITE(10)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
    WRITE(asciiline,'(A,6(I8),A,3(E20.12),A,3(E20.12),A)')'  <ImageData WholeExtent="', & 
        & icell(1) - 1, icell(1) + ncell(1) - 1, &
        & icell(2) - 1, icell(2) + ncell(2) - 1, &
        & icell(3) - 1, icell(3) + ncell(3) - 1, &
        & '" Ghostlevel="0" Origin="',patch%xmin(1),&
        & patch%xmin(2),patch%xmin(3), &
        & '" Spacing="',patch%dx(1),patch%dx(2),patch%dx(3),'">'

    WRITE(10)TRIM(asciiline)
    WRITE(asciiline,'(A,6(I8),A)') '<Piece Extent="', &
        & icell(1) - 1, icell(1) + ncell(1) - 1, &
        & icell(2) - 1, icell(2) + ncell(2) - 1, &
        & icell(3) - 1, icell(3) + ncell(3) - 1,'">'

    WRITE(10)TRIM(asciiline)

    WRITE(10)'<PointData>'
    WRITE(10)'</PointData>'
    WRITE(10)'<CellData>'


!---------------------------------------------------------------------------------!
! List field variables and calculate the number of bytes
!---------------------------------------------------------------------------------!
    WRITE(asciiline,'(3A,I10,A)')'<DataArray type="',float_type, &
      & '" Name="vort" NumberOfComponents="3" format="appended" offset="', &
      & nbytes,'"/>'
    WRITE(10)TRIM(asciiline)
    nbytes = nbytes + 3*ndata*prec_output

    WRITE(asciiline,'(3A,I10,A)')'<DataArray type="',float_type, &
      & '" Name="vel" NumberOfComponents="3" format="appended" offset="', &
      & nbytes,'"/>'
    WRITE(10)TRIM(asciiline)
    nbytes = nbytes + 3*ndata*prec_output

    WRITE(10)'</CellData>'
    WRITE(10)'</Piece>' 
    WRITE(10)'</ImageData>'
    WRITE(10)'<AppendedData encoding="raw">'

    WRITE(10) '_'
    WRITE(10) nbytes

!---------------------------------------------------------------------------------!
! Append data
!---------------------------------------------------------------------------------!
! vorticity
    DO k=1, ncell(3)
      DO j=1, ncell(2)
        DO i=1, ncell(1)
          WRITE(10) REAL(mesh%vort(1,i,j,k),prec_output)
          WRITE(10) REAL(mesh%vort(2,i,j,k),prec_output)
          WRITE(10) REAL(mesh%vort(3,i,j,k),prec_output)
        END DO
      END DO
    END DO
! velocity
    DO k=1, ncell(3)
      DO j=1, ncell(2)
        DO i=1, ncell(1)
          WRITE(10) REAL(mesh%vel(1,i,j,k),prec_output)
          WRITE(10) REAL(mesh%vel(2,i,j,k),prec_output)
          WRITE(10) REAL(mesh%vel(3,i,j,k),prec_output)
        END DO
      END DO
    END DO

    WRITE(10)'</VTKFile>'


    CLOSE(10)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE output_mesh

