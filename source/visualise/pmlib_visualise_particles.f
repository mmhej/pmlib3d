!---------------------------------------------------------------------------------!
! pmlib_visualise_particles.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_visualise_particles( outfile,topo,part,max_vort, &
                                    & foc,dist,azi,ele,zoom,ierr,scl)

USE pmlib_mod_particles
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN)                    :: outfile
  TYPE(class_topology),DIMENSION(:),POINTER       :: topo
  TYPE(class_particles)                           :: part
  REAL(MK), INTENT(IN)                            :: max_vort
  REAL(MK),DIMENSION(pmlib_ndim), INTENT(IN)      :: foc
  REAL(MK), INTENT(IN)                            :: dist
  REAL(MK), INTENT(IN)                            :: azi
  REAL(MK), INTENT(IN)                            :: ele
  REAL(MK), INTENT(IN)                            :: zoom
  INTEGER, INTENT(OUT)                            :: ierr
  REAL(MK),INTENT(IN),OPTIONAL                    :: scl

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  CHARACTER(LEN=20)                      :: caller = 'pmlib_visualise_part'
  REAL(MK)                               :: PI
  INTEGER                                :: i,j,k,l,p

  REAL(MK),DIMENSION(3)                  :: eye_xpos
  REAL(MK)                               :: eye_ele, eye_azi, eye_dist

  REAL(MK)                               :: opac, gauss

  INTEGER                                :: npix
  REAL(MK)                               :: px, py
  REAL(MK),DIMENSION(3,2)                :: pix_bvec
  REAL(MK)                               :: pix_dist

  REAL(MK)                               :: eps,part_radius
  REAL(MK),DIMENSION(:),POINTER          :: part_dist => NULL()
  REAL(MK),DIMENSION(:),POINTER          :: part_vort => NULL()
  REAL(MK),DIMENSION(:,:),POINTER        :: part_pos  => NULL()

  INTEGER                                :: ix, iy, icol
  REAL(MK)                               :: epx, epy, epz, edotp
  REAL(MK)                               :: dx_mod, dy_mod

  INTEGER,DIMENSION(:),POINTER           :: plotorder
  REAL(MK),DIMENSION(:),POINTER          :: proc_dist

  INTEGER, PARAMETER                     :: ncolor = 256
  INTEGER, DIMENSION(2),PARAMETER        :: npixel   = (/ 512, 384 /)
  REAL(MK), DIMENSION(:,:,:),POINTER     :: pix_rgb
  REAL(MK), DIMENSION(:,:,:,:),POINTER   :: gl_pix_rgb

  REAL(MK), DIMENSION(4,ncolor)          :: palette
  REAL(MK), DIMENSION(2)                 :: pix_dx, pix_xmin


  REAL(MK),DIMENSION(:,:),POINTER                 :: palette_def
  INTEGER,DIMENSION(:),POINTER                    :: col_def
  INTEGER                                         :: ncol_def
  REAL(MK)                                        :: frac

  INTEGER                                         :: itemp
  REAL(MK)                                        :: temp
  LOGICAL                                         :: swapped

!---------------------------------------------------------------------------------!
! Initiate subroutine
!---------------------------------------------------------------------------------!
  ierr = 0
  PI = ACOS(-1.0_MK)

!---------------------------------------------------------------------------------!
! Define palette (value,red,green,blue,opacity)
!---------------------------------------------------------------------------------!
  ncol_def = 5
  ALLOCATE(palette_def(ncol_def,5))
  ALLOCATE(col_def(ncol_def))

!  palette_def(1,:) = (/ 0.00_MK,   0.0_MK,   0.0_MK,   0.0_MK, 0.00_MK /)
!  palette_def(2,:) = (/ 0.50_MK, 175.0_MK, 221.0_MK, 233.0_MK, 0.005_MK /)
!  palette_def(3,:) = (/ 0.75_MK, 175.0_MK, 221.0_MK, 233.0_MK, 0.05_MK /)
!  palette_def(4,:) = (/ 0.90_MK, 255.0_MK, 221.0_MK,  85.0_MK, 0.20_MK /)
!  palette_def(5,:) = (/ 1.00_MK, 255.0_MK,  42.0_MK,  42.0_MK, 0.30_MK /)

  palette_def(1,:) = (/ 0.00_MK,   0.0_MK,   0.0_MK,   0.0_MK, 0.00_MK /)
  palette_def(2,:) = (/ 0.25_MK, 175.0_MK, 221.0_MK, 233.0_MK, 0.05_MK /)
  palette_def(3,:) = (/ 0.50_MK, 175.0_MK, 221.0_MK, 233.0_MK, 0.10_MK /)
  palette_def(4,:) = (/ 0.75_MK, 255.0_MK, 221.0_MK,  85.0_MK, 0.15_MK /)
  palette_def(5,:) = (/ 1.00_MK, 255.0_MK,  42.0_MK,  42.0_MK, 0.20_MK /)


  col_def(1) = 0
  DO i = 2,ncol_def
    col_def(i) = INT(REAL(ncolor,MK) * palette_def(i,1))
  END DO

  icol = 0
  DO i = 1,ncol_def-1
    DO j = col_def(i)+1,col_def(i+1)
      icol = icol + 1

      frac  = REAL(j-col_def(i),MK)/REAL(col_def(i+1) - col_def(i),MK)

      palette(1,icol) = frac*palette_def(i+1,2) + (1-frac)*palette_def(i,2)
      palette(2,icol) = frac*palette_def(i+1,3) + (1-frac)*palette_def(i,3)
      palette(3,icol) = frac*palette_def(i+1,4) + (1-frac)*palette_def(i,4)
      palette(4,icol) = frac*palette_def(i+1,5) + (1-frac)*palette_def(i,5)

   END DO
  END DO

  DEALLOCATE(palette_def)
  DEALLOCATE(col_def)

!---------------------------------------------------------------------------------!
! Setup the picture
!---------------------------------------------------------------------------------!
  pix_dx(1) = 1.0/REAL(npixel(1),MK)
  pix_dx(2) = 0.75/REAL(npixel(2),MK)

  pix_xmin(1) = -AINT(REAL(npixel(1),MK)/2.0_MK)*pix_dx(1)
  pix_xmin(2) = -AINT(REAL(npixel(2),MK)/2.0_MK)*pix_dx(2)

!---------------------------------------------------------------------------------!
! setup eye and picture
!
!       | y
!       |
!       |
!       o---------- x
!        \
!         \
!          \ z
!---------------------------------------------------------------------------------!
! picture center coordinate in spherical coordinates with focus as origo
  eye_dist  = dist
  eye_ele   = (90.0_MK - ele)*2.0_MK*PI/360.0_MK
  eye_azi   = azi*2.0_MK*PI/360.0_MK

  eye_xpos(1) = eye_dist * SIN(eye_ele) * SIN(eye_azi)
  eye_xpos(2) = eye_dist * COS(eye_ele)
  eye_xpos(3) = eye_dist * SIN(eye_ele) * COS(eye_azi)

  pix_dist = zoom 

! the picture basis vectors are given by the azimuthal and elevation 
! unit vectors
  pix_bvec(:,1) = (/ -COS(eye_azi), 0.0_MK ,  SIN(eye_azi) /)
  pix_bvec(:,2) = (/ COS(eye_ele)*SIN(eye_azi), -SIN(eye_ele), &
                   & COS(eye_ele)*COS(eye_azi) /)

!---------------------------------------------------------------------------------!
! Calculate particle-to-eye vector and vorticity magnitude
!---------------------------------------------------------------------------------!
  ALLOCATE( part_dist( topo(mpi_rank)%npart ) )
  ALLOCATE( part_vort( topo(mpi_rank)%npart ) )
  ALLOCATE( part_pos(3, topo(mpi_rank)%npart ) )

  DO i = 1,topo(mpi_rank)%npart

    part_pos(1,i) = part%pos(1,i)
    part_pos(2,i) = part%pos(2,i)
    part_pos(3,i) = part%pos(3,i)

! eye to particle vector
    epx = (part_pos(1,i) - foc(1)) - eye_xpos(1)
    epy = (part_pos(2,i) - foc(2)) - eye_xpos(2)
    epz = (part_pos(3,i) - foc(3)) - eye_xpos(3)

    part_dist(i) = SQRT( epx**2 + epy**2 + epz**2 )

! Calculate particle magnitude
    IF(part_dist(i) .GE. pix_dist)THEN
      part_vort(i) = 0.0_MK
      DO j = 1,pmlib_nvort
        part_vort(i) = part_vort(i) + part%vort(j,i)**2
      END DO
      part_vort(i) = SQRT(part_vort(i))/max_vort

      IF(part_vort(i) .GT. 1.0_MK)THEN
        part_vort(i) = 1.0_MK
      END IF
    ELSE
      part_vort(i) = 0.0_MK
    END IF

  END DO

!---------------------------------------------------------------------------------!
! Sort particles according to distance in descending order
!---------------------------------------------------------------------------------!
  CALL pmlib_visualise_sort_particles(part_dist,part_vort,part_pos,ierr)
  IF(ierr .NE. 0)THEN 
    CALL pmlib_write(mpi_rank,caller,'Failed to sort particles.')
    GO TO 9999
  END IF

!---------------------------------------------------------------------------------!
! Setup the picture background
!---------------------------------------------------------------------------------!
  ALLOCATE( pix_rgb(4,npixel(1),npixel(2)) )
  ALLOCATE( gl_pix_rgb(4,npixel(1),npixel(2),mpi_nproc) )

  DO j = 1,npixel(2)
    DO i = 1,npixel(1)
      pix_rgb(1,i,j) = palette(1,1)
      pix_rgb(2,i,j) = palette(2,1)
      pix_rgb(3,i,j) = palette(3,1)
      pix_rgb(4,i,j) = palette(4,1)
    END DO
  END DO

!---------------------------------------------------------------------------------!
! Project particles onto picture
!---------------------------------------------------------------------------------!
  DO i = 1,topo(mpi_rank)%npart

    eps = pmlib_regularisation_radius*MAXVAL(topo(mpi_rank)%dx)

! eye to particle vector
    epx = eye_xpos(1) - (part_pos(1,i) - foc(1)) 
    epy = eye_xpos(2) - (part_pos(2,i) - foc(2))
    epz = eye_xpos(3) - (part_pos(3,i) - foc(3)) 

! calculate projection of vectors
    edotp   = ( eye_xpos(1)* epx + eye_xpos(2)* epy + eye_xpos(3)* epz )/eye_dist

! normalise eye to particle vector to intersection of pic
    epx = epx * pix_dist/edotp
    epy = epy * pix_dist/edotp
    epz = epz * pix_dist/edotp

    px = epx*pix_bvec(1,1) + epy*pix_bvec(2,1) + epz*pix_bvec(3,1)
    py = epx*pix_bvec(1,2) + epy*pix_bvec(2,2) + epz*pix_bvec(3,2)

    IF(PRESENT(scl))THEN
      part_radius = scl*0.2_MK*MAX(eps*pix_dist/edotp,MAXVAL(pix_dx))
    ELSE
      part_radius = 0.2_MK*MAX(eps*pix_dist/edotp,MAXVAL(pix_dx))
    END IF

    npix = CEILING(part_radius/MAXVAL(pix_dx))

! find index of south-west cell
    ix = 1 + FLOOR((px - (pix_xmin(1) + 0.5_MK*pix_dx(1)) )/pix_dx(1))
    iy = 1 + FLOOR((py - (pix_xmin(2) + 0.5_MK*pix_dx(2)) )/pix_dx(2))    

    IF( ix .GE. 1 .AND. ix .LE. npixel(1) .AND. &
      & iy .GE. 1 .AND. iy .LE. npixel(2) )THEN
       
      icol = FLOOR(part_vort(i)*255.0_MK) + 1

      DO j = -2*(npix+1),2*npix
        DO k = -2*(npix+1),2*npix

          IF( ix + j .GE. 1 .AND. ix + j .LE. npixel(1)  .AND. &
            & iy + k .GE. 1 .AND. iy + k .LE. npixel(2) )THEN

            dx_mod = (ix + j)*pix_dx(1) + (pix_xmin(1) - 0.5_MK*pix_dx(1)) - px
            dy_mod = (iy + k)*pix_dx(2) + (pix_xmin(2) - 0.5_MK*pix_dx(2)) - py

            gauss = EXP(-(dx_mod**2+dy_mod**2)/(2.0_MK*part_radius**2))

! Alpha compositing

              opac = 1 - ( 1 - gauss*palette(4,icol) )*( 1 - pix_rgb(4,ix+j,iy+k) )

              pix_rgb(1,ix+j,iy+k) = ( gauss*palette(4,icol)*palette(1,icol) + &
                           & (1 - gauss*palette(4,icol))* &
                           & pix_rgb(4,ix+j,iy+k)*pix_rgb(1,ix+j,iy+k) )/opac
              pix_rgb(2,ix+j,iy+k) = ( gauss*palette(4,icol)*palette(2,icol) + &
                           & (1 - gauss*palette(4,icol))* &
                           & pix_rgb(4,ix+j,iy+k)*pix_rgb(2,ix+j,iy+k) )/opac
              pix_rgb(3,ix+j,iy+k) = ( gauss*palette(4,icol)*palette(3,icol) + &
                           & (1 - gauss*palette(4,icol))* &
                           & pix_rgb(4,ix+j,iy+k)*pix_rgb(3,ix+j,iy+k) )/opac

              pix_rgb(4,ix+j,iy+k) = opac

          END IF

        END DO
      END DO
  
    END IF
  END DO ! particles

!---------------------------------------------------------------------------------!
! Deallocate particles
!---------------------------------------------------------------------------------!
  DEALLOCATE(part_dist)
  DEALLOCATE(part_vort)

!---------------------------------------------------------------------------------!
! Communicate pictures
!---------------------------------------------------------------------------------!
  CALL MPI_GATHER( pix_rgb,4*npixel(1)*npixel(2),mpi_prec_real,gl_pix_rgb, &
                 & 4*npixel(1)*npixel(2), mpi_prec_real,0,mpi_comm,ierr  )
  IF(ierr .NE. 0)THEN 
    CALL pmlib_write(mpi_rank,caller,'Failed to communicate pictures.')
    GO TO 9999
  END IF

!---------------------------------------------------------------------------------!
! Merge pictures
!---------------------------------------------------------------------------------!
  IF(mpi_rank .EQ. 0)THEN
!---------------------------------------------------------------------------------!
! Find plotting order
!---------------------------------------------------------------------------------!
    ALLOCATE( plotorder(mpi_nproc) )
    ALLOCATE( proc_dist(mpi_nproc) )
    DO i = 0,mpi_nproc-1
      plotorder(i+1) = i+1
      proc_dist(i+1) = SQRT( &
         & ( eye_xpos(1) - ( topo(i)%xmax(1) + topo(i)%xmin(1) )/2.0_MK )**2 + &
         & ( eye_xpos(2) - ( topo(i)%xmax(2) + topo(i)%xmin(2) )/2.0_MK )**2 + &
         & ( eye_xpos(3) - ( topo(i)%xmax(3) + topo(i)%xmin(3) )/2.0_MK )**2 )
    END DO

! Bubble sort
    DO j = mpi_nproc-1, 1, -1
      swapped = .FALSE.
      DO i = 1, j
        IF (proc_dist(i) < proc_dist(i+1)) THEN
          temp = proc_dist(i)
          proc_dist(i) = proc_dist(i+1)
          proc_dist(i+1) = temp

          itemp = plotorder(i)
          plotorder(i) = plotorder(i+1)
          plotorder(i+1) = itemp

          swapped = .TRUE.
        END IF
      END DO
      IF (.NOT. swapped) EXIT
    END DO

!---------------------------------------------------------------------------------!
! 
!---------------------------------------------------------------------------------!
    DO j = 1,npixel(2)
      DO i = 1,npixel(1)

! reset to background
        pix_rgb(1,i,j) = palette(1,1)
        pix_rgb(2,i,j) = palette(2,1)
        pix_rgb(3,i,j) = palette(3,1)
        pix_rgb(4,i,j) = palette(4,1)

        DO l = 1,mpi_nproc
          k = plotorder(l)

! Alpha compositing
          opac = 1 - ( 1 - gl_pix_rgb(4,i,j,k) )*( 1 - pix_rgb(4,i,j) )

          pix_rgb(1,i,j) = ( gl_pix_rgb(4,i,j,k)*gl_pix_rgb(1,i,j,k) + &
                         & (1 - gl_pix_rgb(4,i,j,k) )* &
                         & pix_rgb(4,i,j)*pix_rgb(1,i,j) )/opac
          pix_rgb(2,i,j) = ( gl_pix_rgb(4,i,j,k)*gl_pix_rgb(2,i,j,k) + &
                         & (1 - gl_pix_rgb(4,i,j,k) )* &
                         & pix_rgb(4,i,j)*pix_rgb(2,i,j) )/opac
          pix_rgb(3,i,j) = ( gl_pix_rgb(4,i,j,k)*gl_pix_rgb(3,i,j,k) + &
                         & (1 - gl_pix_rgb(4,i,j,k) )* &
                         & pix_rgb(4,i,j)*pix_rgb(3,i,j) )/opac
          pix_rgb(4,i,j) = opac
  
        END DO

      END DO
    END DO

!---------------------------------------------------------------------------------!
! Output post script file
!---------------------------------------------------------------------------------!
    CALL pmlib_visualise_postscript(pix_rgb,outfile,ierr)
    IF(ierr .NE. 0)THEN 
      CALL pmlib_write(mpi_rank,caller,'Failed to output postscript file.')
      GO TO 9999
    END IF

    DEALLOCATE(plotorder)
    DEALLOCATE(proc_dist)
  END IF

!---------------------------------------------------------------------------------!
! Deallocate local pictures
!---------------------------------------------------------------------------------!
  DEALLOCATE(pix_rgb)
  DEALLOCATE(gl_pix_rgb)

!---------------------------------------------------------------------------------!
! Return
!---------------------------------------------------------------------------------!
 9999 CONTINUE
RETURN

END SUBROUTINE pmlib_visualise_particles
