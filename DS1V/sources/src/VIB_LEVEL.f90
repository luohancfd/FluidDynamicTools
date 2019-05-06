!
!*****************************************************************************
!
SUBROUTINE VIB_LEVEL(B,IV,KV,KS)
!
!--for a given energy EVIB and vibrational mode KV of species KS,
!--calculate the corresponding truncated vibrational level
!
USE GAS
USE CALC
USE MOLECS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: KS,IV,KV
REAL(KIND=8) :: EVIB,A(0:IVMODEL(KS,2)),B
!
EVIB=B*1.000001d0
!
IF (IVMODEL(KS,1) == 0) THEN
  IV=EVIB/(BOLTZ*SPVM(1,KV,KS)) !SHO model
ELSE
  A(:)=1.d0
  A(:)=EVIB-VIBEN(:,KS)
  !find the first, smallest positive value in A vector
  IV=MINLOC(A,DIM=1,MASK=A.GT.0.d0)-1 !subtract 1 because index starts at zero
  IF (IV < 0)  IV=0
  IF (IV > IVMODEL(KS,2)) IV=IVMODEL(KS,2)
END IF
!
RETURN
END SUBROUTINE VIB_LEVEL
