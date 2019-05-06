!
!*****************************************************************************
!
SUBROUTINE INELASTIC_RT(KS, JS, EC, P, IDT)
!
!--calculate P=Erot(KS)/Ec for RT procedure
!
USE GAS,ONLY : INONVHS, SPM, ISPR
USE EXPCOL,ONLY : EXPCOL_RT
IMPLICIT NONE
INTEGER, INTENT(IN) :: KS, JS
INTEGER :: IDT
REAL(KIND=8), INTENT(IN) :: EC
REAL(KIND=8) :: P, RANF, PMAX

IF (INONVHS(KS, JS) == 2) THEN
  P =  EXPCOL_RT(KS,JS,EC,ISPR(1,KS),IDT)
ELSE
  ! use regular LB model
  IF (ISPR(1,KS) == 2) THEN
    CALL ZGF(RANF, IDT)
    P = 1.0d0 - RANF**(1.0d0 / (2.5d0 - SPM(3,KS,JS))) !EQN 5.46
  ELSE
    CALL LBS(DBLE(ISPR(1,KS))*0.5d0 - 1.0d0, 1.5d0 - SPM(3,KS,JS), P, IDT)
  END IF
END IF

END SUBROUTINE INELASTIC_RT
