!
!******************************************************************************
!
FUNCTION GAM(X)
!
!--calculates the Gamma function of X.
!
IMPLICIT NONE
REAL(8) ::X,A,Y,GAM
A=1.D0
Y=X
IF (Y < 1.D0) THEN
  A=A/Y
ELSE
  Y=Y-1.D0
  DO WHILE (Y >= 1.D0)
    A=A*Y
    Y=Y-1.D0
  END DO
END IF
GAM=A*(1.D0-0.5748646D0*Y+0.9512363D0*Y**2-0.6998588D0*Y**3+  &
       0.4245549D0*Y**4-0.1010678D0*Y**5)
!
RETURN
!
END FUNCTION GAM
