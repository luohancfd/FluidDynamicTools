!
!***************************************************************************
!
SUBROUTINE IDEAL_GAS
!
USE GAS
USE CALC
!
IMPLICIT NONE
INTEGER :: J !--isebasti: J included
!
IF(IRM == 1) THEN
  MSP=1; MMRM=2; MMVM=0
END IF
!
IF(IRM >= 2) THEN
  MSP=1; MMRM=2; MMVM=1
END IF
!
IF(IRM == 4) THEN
  MSP=1; MMRM=3; MMVM=3
END IF
!
IF(IRM == 5) THEN
  MSP=2; MMRM=3; MMVM=3
END IF
!
IF(IRM == 6) THEN
  MSP=2; MMRM=3; MMVM=3
END IF
!
MNRE=0
MTBP=0
MEX=0
MMEX=0
MNSR=0
!
CALL ALLOCATE_GAS
!
IF(IRM == 1) THEN !N2 (Trans+Rot; Zr const)
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
END IF
!
IF(IRM == 2) THEN !N2 (Trans+Rot+Vib; Zv const)
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=5.D00
SPVM(3,1,1)=-1.D00
SPVM(4,1,1)=113500.D00
END IF
!
IF(IRM == 3) THEN !N2 (Trans+Rot+Vib; Zv(T))
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=52600.D00
SPVM(3,1,1)=3371.D00
SPVM(4,1,1)=113500.D00
END IF
!
IF(IRM == 4) THEN !H2O (Trans+Rot+Vib; Zv(T))
J=1
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=0.75D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
END IF
!
IF(IRM == 5) THEN !H2O+H2O (Trans+Rot+Vib; Zv const)
J=1
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.0D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
!
J=2
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.0D00        !-estimate from SMILE
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
END IF
!
IF(IRM == 6) THEN !H2O+H (Trans+Rot+Vib)
J=1  !--species 1 is water vapor H2O
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.00D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=1.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=10. !20000.D00   !--estimate
SPVM(3,1,J)=-1. !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=10. !20000.D00   !--estimate
SPVM(3,2,J)=-1. !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=10. !20000.D00   !--estimate
SPVM(3,3,J)=-1. !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
!
J=2  !--species 2 is atomic hydrogen H
SP(1,J)=2.5D-10      !--estimate
SP(2,J)=273.D00
SP(3,J)=0.8D00
SP(4,J)=1.D00
SP(5,J)=1.67D-27
ISPR(1,J)=0
ISPV(J)=0
END IF
!
RETURN
END SUBROUTINE IDEAL_GAS
