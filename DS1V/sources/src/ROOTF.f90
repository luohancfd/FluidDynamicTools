!
!*****************************************************************************
!
FUNCTION ROOTF(I,J,L,M,A,X)
!
!--evaluate the functions required in the root finding methods
!
USE CALC
USE GAS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: I,J,K,L,M,MAXLEV,IVIB
REAL(KIND=8) :: ROOTF,A,X,QSUM,QVIB,EVIB,AEVIB,TEMP
REAL(KIND=8) :: VDOF
!
IF (I == 1) THEN !to calculate TVIB(AEVIB)
  K=1            !assuming single vib mode
  AEVIB=A        !average vib energy
  TEMP=X         !vib temperature
  QSUM=0.d0
  QVIB=0.d0
  MAXLEV=IVMODEL(J,2)
  DO IVIB=0,MAXLEV
    CALL VIB_ENERGY(EVIB,IVIB,K,J)
    QSUM=QSUM+EVIB*DEXP(-EVIB/BOLTZ/TEMP)
    QVIB=QVIB+DEXP(-EVIB/BOLTZ/TEMP)
  END DO
  ROOTF=AEVIB-QSUM/QVIB
END IF
!
IF (I == 2) THEN !to calculate Zcorrection according to Gimelshein PoF 14 (2012)
  K=1            !assuming single vib mode
  TEMP=X         !equilibrium temperature
  VDOF=2.d0*(SPVM(1,K,J)/X)/(DEXP(SPVM(1,K,J)/X)-1.d0) !vibrational dof
  ROOTF=A-(VDOF+(5.d0-2.d0*SPM(3,L,M)))*X
END IF
!
RETURN
END FUNCTION ROOTF
