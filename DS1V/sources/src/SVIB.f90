!
!*****************************************************************************
!
SUBROUTINE SVIB(L,TEMP,IVIB,K,IDT)
!
!--sets a typical vibrational state at temp. T of mode K of species L
!
USE GAS
USE CALC
USE MOLECS
!
IMPLICIT NONE
!
INTEGER :: II,K,L,IDT,IVIB,MAXLEVEL !included IDT
REAL(KIND=8) :: EVIB,TEMP,RANF,QVIB,PROB  !--isebasti: included RANF
!
IF (TEMP < 5.d0) THEN
  IVIB=0
  RETURN
END IF
!
!old approach (limited for SHO)
!CALL ZGF(RANF,IDT)
!IVIB=-DLOG(RANF)*TEMP/SPVM(1,K,L) !eqn(11.24)
!
!calculate the vibrational partition function
QVIB=0.d0
MAXLEVEL=IVMODEL(L,2)
DO IVIB=0,MAXLEVEL
  CALL VIB_ENERGY(EVIB,IVIB,K,L)
  QVIB=QVIB+DEXP(-EVIB/(BOLTZ*TEMP))
END DO
!
!sample a level from Boltzmann distribution
II=0
DO WHILE (II == 0)
  CALL ZGF(RANF,IDT)
  IVIB=RANF*(MAXLEVEL+0.99999999D00)
  CALL VIB_ENERGY(EVIB,IVIB,K,L)
  PROB=DEXP(-EVIB/(BOLTZ*TEMP))/QVIB
  CALL ZGF(RANF,IDT)
  IF (PROB > RANF) II=1
END DO
!
RETURN
END SUBROUTINE SVIB
