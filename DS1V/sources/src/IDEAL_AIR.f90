!
!***************************************************************************
!
SUBROUTINE IDEAL_AIR
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=1
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.07D-10
SP(2,1)=273.
SP(3,1)=0.77
SP(4,1)=1.
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.
SP(1,2)=4.17D-10
SP(2,2)=273.
SP(3,2)=0.74
SP(4,2)=1.
SP(5,2)=4.65D-26
ISPR(1,2)=2
ISPR(2,2)=0
SPR(1,2)=5.
RETURN
END SUBROUTINE IDEAL_AIR
