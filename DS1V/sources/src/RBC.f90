!
!*****************************************************************************
!
SUBROUTINE RBC(XI,DX,DY,DZ,R,S)
!
!--calculates the trajectory fraction S from a point at radius XI with
!----displacements DX, DY, and DZ to a possible intersection with a
!----surface of radius R, IFX=1, 2 for cylindrical, spherical geometry
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
!
IMPLICIT NONE
!
REAL(KIND=8) :: A,B,C,XI,DX,DY,DZ,R,S,DD,S1,S2
!
DD=DX*DX+DY*DY
IF (IFX == 2) DD=DD+DZ*DZ
B=XI*DX/DD
C=(XI*XI-R*R)/DD
A=B*B-C
IF (A >= 0.D00) THEN
!--find the least positive solution to the quadratic
  A=DSQRT(A)
  S1=-B+A
  S2=-B-A
  IF (S2 < 0.D00) THEN
    IF (S1 > 0.D00) THEN
      S=S1
    ELSE
      S=2.D00
    END IF
  ELSE IF (S1 < S2) THEN
    S=S1
  ELSE
    S=S2
  END IF
ELSE
  S=2.D00
!--setting S to 2 indicates that there is no intersection
END IF
!
RETURN
!
END SUBROUTINE RBC
