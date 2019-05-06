!
!***************************************************************************
!
SUBROUTINE REAL_OXYGEN
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=1
MMVM=1
MNRE=0
MTBP=0
MEX=0
MMEX=0
MNSR=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.07D-10
SP(2,1)=273.D00
SP(3,1)=0.77D00
SP(4,1)=1.D00
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0             ! 0,1 for constant,polynomial rotational relaxation collision number
SPR(1,1)=5.             ! the collision number or the coefficient of temperature in the polynomial (if a polynomial, the coeff. of T^2 is in spr_db(3  )
ISPV(1)=1               ! the number of vibrational modes
SPVM(1,1,1)=2256.D00          ! the characteristic vibrational temperature
SPVM(2,1,1)=90000.D00        ! a constant Zv, or the reference Zv
IF (IZV == 1) SPVM(2,1,1)=5.D00*SPVM(2,1,1)   !--to allow for the reduction in relaxation time when based on quantized collision temperature
SPVM(3,1,1)=2256.D00        ! -1 for a constant Zv, or the reference temperature
SPVM(4,1,1)=59500.D00       ! characteristic dissociation temperature
ISPVM(1,1,1)=2
ISPVM(2,1,1)=2
!
!--species 2 is atomic oxygen
SP(1,2)=3.D-10
SP(2,2)=273.D00
SP(3,2)=0.8D00
SP(4,2)=1.D00
SP(5,2)=2.656D-26
ISPR(1,2)=0
ISPV(2)=0     !--must be set!
!
!--set data needed for recombination
ISPRC=0
ISPRC(2,2)=1              !--O+O -> O2  recombined species code for an O+O recombination
SPRC(1,2,2,1)=0.04D00     !--UNADJUSTED
IF (ITCV == 1) SPRC(1,2,2,1)=SPRC(1,2,2,1)*0.23
SPRC(2,2,2,1)=-1.3D00
SPRC(1,2,2,2)=0.07D00     !--UNADJUSTED
IF (ITCV == 1) SPRC(1,2,2,2)=SPRC(1,2,2,2)*0.28
SPRC(2,2,2,2)=-1.2D00
!
NSPEX=0
SPEX=0.D00
ISPEX=0
!
RETURN
END SUBROUTINE REAL_OXYGEN
