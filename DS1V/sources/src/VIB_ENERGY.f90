!
!*****************************************************************************
!
SUBROUTINE VIB_ENERGY(EVIB,IV,KV,KS)
!
!--for a given vibrational level IV and mode KV of species KS,
!--calculate the corresponding truncated vibrational energy EVIB
!
USE GAS
USE CALC
USE MOLECS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: KS,IV,KV
REAL(KIND=8) :: EVIB
!
IF (IVMODEL(KS,1) == 0) THEN
  EVIB=IV*BOLTZ*SPVM(1,KV,KS) !SHO model
ELSE
  EVIB=VIBEN(IV,KS) !QCT for N2 and O2 and NO
END IF
!
RETURN
END SUBROUTINE VIB_ENERGY
