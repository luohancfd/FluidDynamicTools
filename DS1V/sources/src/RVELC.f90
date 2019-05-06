!
!*****************************************************************************
!
SUBROUTINE RVELC(U,V,VMP,IDT)
!
USE CALC
!
IMPLICIT NONE
!
!--generates two random velocity components U and V in an equilibrium
!--gas with most probable speed VMP
INTEGER :: IDT !included IDT
REAL(KIND=8) :: U,V,VMP,A,B,RANF !--isebasti: included RANF
!
CALL ZGF(RANF,IDT)
A=DSQRT(-DLOG(RANF))
CALL ZGF(RANF,IDT)
B=DPI*RANF
U=A*DSIN(B)*VMP
V=A*DCOS(B)*VMP
RETURN
!
END SUBROUTINE RVELC
