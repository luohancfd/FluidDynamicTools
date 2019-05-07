!
!*****************************************************************************
!
SUBROUTINE VHSS(LS,MS,BR,VR,VRC,VRCP,IDT)
!
!--calculate new scattering angles based on VHS/VSS model
! BR is (b/bmax)**2
!
USE GAS
USE CALC
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: LS,MS,IDT
REAL(8), INTENT(IN) :: BR
REAL(KIND=8) :: A,B,C,D,OC,SD,VR,VRC(3),VRCP(3),RANF
!
IF (ABS(SPM(8,LS,MS)-1.) < 1.d-3) THEN
!--use the VHS logic
  IF (BR < 0) THEN
    CALL ZGF(RANF,IDT)
  ELSE
    RANF = BR
  END IF
  B=2.D00*RANF-1.D00
!--B is the cosine of a random elevation angle
  A=DSQRT(1.D00-B*B)
  VRCP(1)=B*VR
  CALL ZGF(RANF,IDT)
  C=2.D00*PI*RANF
!--C is a random azimuth angle
  VRCP(2)=A*DCOS(C)*VR
  VRCP(3)=A*DSIN(C)*VR
ELSE
!--use the VSS logic
  IF (BR < 0) THEN
    CALL ZGF(RANF,IDT)
  ELSE
    RANF = BR
  END IF
  B=2.D00*(RANF**SPM(8,LS,MS))-1.D00  !isebasti: SP(4,1) was used instead of SPM(8,LS,MS)
!--B is the cosine of the deflection angle for the VSS model (eqn (11.8)
  A=DSQRT(1.D00-B*B)
  CALL ZGF(RANF,IDT)
  C=2.D00*PI*RANF
  OC=DCOS(C)
  SD=DSIN(C)
  D=SQRT(VRC(2)**2+VRC(3)**2)
  VRCP(1)=B*VRC(1)+A*SD*D
  VRCP(2)=B*VRC(2)+A*(VR*VRC(3)*OC-VRC(1)*VRC(2)*SD)/D
  VRCP(3)=B*VRC(3)-A*(VR*VRC(2)*OC+VRC(1)*VRC(3)*SD)/D
!--the post-collision rel. velocity components are based on eqn (2.22)
END IF
RETURN
END SUBROUTINE VHSS
