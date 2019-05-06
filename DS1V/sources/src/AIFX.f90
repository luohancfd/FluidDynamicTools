!
!*****************************************************************************
!
SUBROUTINE AIFX(XI,DX,DY,DZ,X,U,V,W,IDT)
!
!--calculates the new radius and realigns the velocity components in
!----cylindrical and spherical flows
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: IDT !included IDT
REAL(KIND=8) :: A,B,C,XI,DX,DY,DZ,X,U,V,W,DR,VR,S,RANF !--isebasti: included RANF
!
IF (IFX == 1) THEN
  DR=DY
  VR=V
ELSE IF (IFX == 2) THEN
  DR=DSQRT(DY*DY+DZ*DZ)
  VR=DSQRT(V*V+W*W)
END IF
A=XI+DX
X=DSQRT(A*A+DR*DR)
S=DR/X
C=A/X
B=U
U=B*C+VR*S
V=-B*S+VR*C
IF (IFX == 2) THEN
  VR=V
  CALL ZGF(RANF,IDT)
  A=DPI*RANF
  V=VR*DSIN(A)
  W=VR*DCOS(A)
END IF
!
RETURN
!
END SUBROUTINE AIFX
