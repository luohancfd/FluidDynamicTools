!
!***************************************************************************
!******************************GAS DATABASE*********************************
!***************************************************************************
!
SUBROUTINE HARD_SPHERE
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
!
SP(1,1)=4.0D-10    !reference diameter
SP(2,1)=273.       !reference temperature
SP(3,1)=0.5        !viscosity-temperature index
SP(4,1)=1.         !reciprocal of VSS scattering parameter (1 for VHS)
SP(5,1)=5.D-26     !mass
ISPR(1,1)=0        !number of rotational degrees of freedom
!
RETURN
END SUBROUTINE HARD_SPHERE
