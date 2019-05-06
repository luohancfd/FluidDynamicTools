!
!*****************************************************************************
!
FUNCTION ERF(S)
!
!--evaluates the error function of S
!
IMPLICIT NONE
REAL(8) :: S,B,C,T,D,ERF
B=DABS(S)
IF (B < 4.D0) THEN
  C=DEXP(-B*B)
  T=1.D0/(1.D0+0.3275911D0*B)
  D=1.D0-(0.254829592D0*T-0.284496736D0*T*T+1.421413741D0*T*T*T- &
    1.453152027D0*T*T*T*T+1.061405429D0*T*T*T*T*T)*C
ELSE
  D=1.D0
END IF
IF (S < 0.D0) D=-D
ERF=D
RETURN
END FUNCTION ERF
