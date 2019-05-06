!
!*****************************************************************************
!
SUBROUTINE SROT(L,TEMP,ROTE,IDT)
!
!--sets a typical rotational energy ROTE of species L
!
USE CALC
USE GAS
!
IMPLICIT NONE
!
INTEGER :: I,L,IDT !included IDT
REAL(KIND=8) :: A,B,ROTE,ERM,TEMP,RANF !--isebasti: included RANF
!
IF (ISPR(1,L).EQ.2) THEN
  CALL ZGF(RANF,IDT)
  ROTE=-DLOG(RANF)*BOLTZ*TEMP
ELSE
  A=0.5D00*ISPR(1,L)-1.D00
  I=0
  DO WHILE (I == 0)
    CALL ZGF(RANF,IDT)
    ERM=RANF*10.D00
!--there is an energy cut-off at 10 kT
    B=((ERM/A)**A)*DEXP(A-ERM)
    CALL ZGF(RANF,IDT)
    IF (B > RANF) I=1
  END DO
  ROTE=ERM*BOLTZ*TEMP
END IF
!
RETURN
!
END SUBROUTINE SROT
