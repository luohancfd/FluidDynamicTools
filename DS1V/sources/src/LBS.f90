!
!*****************************************************************************
!
SUBROUTINE LBS(XMA,XMB,ERM,IDT)
!
!--selects a Larsen-Borgnakke energy ratio using eqn (11.9)
!
IMPLICIT NONE
!
REAL(KIND=8) :: PROB,ERM,XMA,XMB,RANF
INTEGER :: I,N,IDT !--isebasti: included IDT
!
!--I is an indicator
!--PROB is a probability
!--ERM ratio of rotational to collision energy
!--XMA degrees of freedom under selection-1
!--XMB remaining degrees of freedom-1
!
I=0
DO WHILE (I == 0)
  CALL ZGF(RANF,IDT)
  ERM=RANF
  IF ((XMA < 1.D-6).OR.(XMB < 1.D-6)) THEN
!    IF (XMA < 1.E-6.AND.XMB < 1.E-6) RETURN
!--above can never occur if one mode is translational
    IF (XMA < 1.D-6) PROB=(1.D00-ERM)**XMB
    IF (XMB < 1.D-6) PROB=(1.D00-ERM)**XMA
  ELSE
    PROB=(((XMA+XMB)*ERM/XMA)**XMA)*(((XMA+XMB)*(1.D00-ERM)/XMB)**XMB)
  END IF
  CALL ZGF(RANF,IDT)
  IF (PROB > RANF) I=1
END DO
!
RETURN
!
END SUBROUTINE LBS
