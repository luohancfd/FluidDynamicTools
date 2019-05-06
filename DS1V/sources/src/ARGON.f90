!
!***************************************************************************
!
SUBROUTINE ARGON
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=1
MMRM=0
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
SP(1,1)=4.11D-10   !--isebasti: use VSS; VHS value is 4.17D-10
SP(2,1)=273.15
SP(3,1)=0.81
SP(4,1)=1.D0/1.4D0 !--isebasti: use VSS
SP(5,1)=6.63D-26
ISPR(1,1)=0
ISPR(2,1)=0
!
RETURN
END SUBROUTINE ARGON
