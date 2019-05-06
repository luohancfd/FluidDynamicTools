!
!******************************************************************************************************************************************
!
SUBROUTINE HELIUM_XENON
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=0
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=2.33D-10
SP(2,1)=273.
SP(3,1)=0.66
SP(4,1)=1.
SP(5,1)=6.65D-27
ISPR(1,1)=0
ISPR(2,1)=0
SP(1,2)=5.74D-10
SP(2,2)=273.
SP(3,2)=0.85
SP(4,2)=1.
SP(5,2)=21.8D-26
ISPR(1,2)=0
ISPR(2,2)=0
RETURN
END SUBROUTINE HELIUM_XENON
