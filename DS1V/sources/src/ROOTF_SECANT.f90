!
!*****************************************************************************
!
SUBROUTINE ROOTF_SECANT(I,J,L,M,A,B,X2)
!
!--use the root finding secant method
!
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: NSTEP,I,J,L,M
REAL(KIND=8) :: A,B,X1,X2,ROOTF,F1,F2,XERROR,FERROR
!
X1=A                     !1st initial guess
X2=A*1.1d0               !2nd initial guess
F1=ROOTF(I,J,L,M,B,X1)
F2=ROOTF(I,J,L,M,B,X2)
XERROR=1.d0-X1/X2        !initializing error
FERROR=1.d0-F1/F2        !initializing error
NSTEP=0
!
DO WHILE (DABS(XERROR)>1.d-6 .AND. DABS(F2)/=0.d0 .AND. DABS(FERROR)/=0.d0)
  X1=X2                       !update 1st guess
  X2=X2*(1.d0-XERROR/FERROR)  !update 2nd guess; same as X3=X2-(X2-X1)*F2/(F2-F1)
  F1=F2
  F2=ROOTF(I,J,L,M,B,X2)
  XERROR=1.d0-X1/X2
  FERROR=1.d0-F1/F2
  NSTEP=NSTEP+1
  IF (NSTEP > 10000) THEN
    WRITE(*,*) 'ROOTF_SECANT is not converging:',DABS(XERROR),DABS(F2)
    XERROR=0.d0
    X2=A                 !same as initial guess
  END IF
  !WRITE(*,*) X1,X2,DABS(XERROR),DABS(F1),DABS(F2),DABS(FERROR)
END DO
!pause
!
RETURN
END SUBROUTINE ROOTF_SECANT
