!---------------------------------------------------------------------------------!
! pmlib_visualise_sort_particles.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
SUBROUTINE pmlib_visualise_sort_particles(part_dist,part_vort,part_pos,ierr)

USE pmlib_mod_particles
USE pmlib_mod_topology

IMPLICIT NONE
!---------------------------------------------------------------------------------!
! Arguments
!---------------------------------------------------------------------------------!
  REAL(MK),DIMENSION(:),POINTER                   :: part_dist
  REAL(MK),DIMENSION(:),POINTER                   :: part_vort
  REAL(MK),DIMENSION(:,:),POINTER                 :: part_pos

  INTEGER, INTENT(OUT)                            :: ierr

!---------------------------------------------------------------------------------!
! Local variables
!---------------------------------------------------------------------------------!
  INTEGER,PARAMETER                  :: m      = 7
  INTEGER,PARAMETER                  :: nstack = 50
  REAL(MK),PARAMETER                 :: fm     = 7875.0_MK
  REAL(MK),PARAMETER                 :: fa     = 211.0_MK
  REAL(MK),PARAMETER                 :: fc     = 1663.0_MK
  REAL(MK),PARAMETER                 :: fmi    = 1.0_MK/fm

  INTEGER                            :: i, j, k, l
  INTEGER                            :: jstack, ir, iq

  INTEGER                            :: npart

  INTEGER, DIMENSION(nstack)         :: istack

  REAL(MK)                           :: fx
  REAL(MK), DIMENSION(5)              :: a

  ierr = 0

  npart = SIZE(part_dist)

  jstack = 0
  l      = 1
  ir     = npart
  fx     = 0.0_MK

 10 CONTINUE 
  IF (ir - l .LT. m) THEN
    DO j = l + 1, ir
      a(1) = part_dist(j)
      a(2) = part_pos(1,j)
      a(3) = part_pos(2,j)
      a(4) = part_pos(3,j)
      a(5) = part_vort(j)
      DO i = j-1,1,-1
        IF (part_dist(i) .GE. a(1)) GO TO 12
          part_dist(i+1)  = part_dist(i)
          part_pos(1,i+1) = part_pos(1,i)
          part_pos(2,i+1) = part_pos(2,i)
          part_pos(3,i+1) = part_pos(3,i)
          part_vort(i+1)  = part_vort(i)
      END DO
      i = 0

 12 CONTINUE
      part_dist(i+1)  = a(1)
      part_pos(1,i+1) = a(2)
      part_pos(2,i+1) = a(3)
      part_pos(3,i+1) = a(4)
      part_vort(i+1)  = a(5)

    END DO

    IF(jstack .EQ. 0) GO TO 9999
      ir = istack(jstack)
      l  = istack(jstack-1)
      jstack = jstack - 2
  ELSE

    i = l; j = ir
    fx = MOD(fx*fa+fc,fm) !Generate a random integer IQ between L and IR inclusive.
    iq = l + (ir - l + 1)*(fx*fmi)

    a(1) = part_dist(iq)
    a(2) = part_pos(1,iq)
    a(3) = part_pos(2,iq)
    a(4) = part_pos(3,iq)
    a(5) = part_vort(iq)

    part_dist(iq)  = part_dist(l)
    part_pos(1,iq) = part_pos(1,l)
    part_pos(2,iq) = part_pos(2,l)
    part_pos(3,iq) = part_pos(3,l)
    part_vort(iq)  = part_vort(l)

 20 CONTINUE

    IF (j .GT. 0) THEN
 21 CONTINUE
      IF(a(1) .GT. part_dist(j)) THEN
        j = j - 1
        GO TO 21
      END IF
    END IF

    IF (j .LE. i) THEN
      part_dist(i)  = a(1)
      part_pos(1,i) = a(2)
      part_pos(2,i) = a(3)
      part_pos(3,i) = a(4)
      part_vort(i)  = a(5)
      GO TO 30
    END IF

    part_dist(i)  = part_dist(j)
    part_pos(1,i) = part_pos(1,j)
    part_pos(2,i) = part_pos(2,j)
    part_pos(3,i) = part_pos(3,j)
    part_vort(i)  = part_vort(j)

    i = i + 1

 22 CONTINUE
    IF (i .LE. npart) THEN
      IF (a(1) .LT. part_dist(i)) THEN
        i = i + 1
        GO TO 22
      END IF
    END IF

    IF (j .LE. i) THEN
      part_dist(j)  = a(1)
      part_pos(1,j) = a(2)
      part_pos(2,j) = a(3)
      part_pos(3,j) = a(4)
      part_vort(j)  = a(5)

      i = j
      GO TO 30
    END IF
    part_dist(j)  = part_dist(i)
    part_pos(1,j) = part_pos(1,i)
    part_pos(2,j) = part_pos(2,i)
    part_pos(3,j) = part_pos(3,i)
    part_vort(j)  = part_vort(i)

    j = j - 1  
    GO TO 20

 30 CONTINUE
    jstack = jstack + 2

    IF (jstack .GT. nstack) THEN
      WRITE(*,*) 'Parameter nstack is too small!'
      ierr = -1
      GO TO 9999
    END IF

    IF (ir-i .GE. i-l) THEN
      istack(jstack)     = ir
      istack(jstack - 1) = i + 1
      ir = i - 1
    ELSE
      istack(jstack) = i - 1
      istack(jstack - 1) = l
      l = i + 1
    END IF
           
  END IF

  GO TO 10

 9999 CONTINUE

END SUBROUTINE pmlib_visualise_sort_particles
