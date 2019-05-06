!
!*****************************************************************************
!
SUBROUTINE ENERGY(M,I,TOTEN)
!
!--calculate the total energy (all molecules if I=0, otherwise molecule I)
!--used for diagnostic purposes only
!
USE MOLECS
USE GAS
USE CALC
!
IMPLICIT NONE
!
INTEGER :: K,L,N,I,II,M,IV,KV,J
REAL(KIND=8) :: A,TOTEN,TOTENI,HF(MSP)
!
TOTEN=0.
!
IF ((M == 6).OR.(M == 4)) THEN
  IF (I == 0) THEN
    DO N=1,NM
      IF (IPCELL(N) > 0) THEN
        L=IPSP(N)
        TOTENI=TOTEN
        IF (M == 6) THEN
          IF (L==3) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
          IF (L==4) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,2)
          IF (L==5) TOTEN=TOTEN+1.49D-19
        END IF
        IF (M == 4) THEN
          IF ((L==2).OR.(L == 3)) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
        END IF
        TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
        IF (ISPV(L) > 0) THEN
          DO KV=1,ISPV(L)
          J=IPVIB(KV,N)
          IF (J <0) THEN
            J=-J
            IF (J == 99999) J=0
          END IF
            TOTEN=TOTEN+DFLOAT(J)*BOLTZ*SPVM(1,KV,L)
          END DO
        END IF
      END IF
      IF ((TOTEN-TOTENI) > 1.D-16) WRITE (*,*) 'MOL',N,' ENERGY',TOTEN-TOTENI
    END DO
!
!    WRITE (9,*) 'Total Energy =',TOTEN,NM
!    WRITE (*,*) 'Total Energy =',TOTEN,NM
  ELSE
    N=I
    IF (IPCELL(N) > 0) THEN
      L=IPSP(N)
      IF (M == 6) THEN
        IF (L==3) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
        IF (L==4) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,2)
        IF (L==5) TOTEN=TOTEN+1.49D-19
      END IF
      IF ((M == 4).OR.(M > 8)) THEN
!        IF ((L==2).OR.(L == 3)) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
      END IF
      TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
      IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
      IF (ISPV(L) > 0) THEN
        DO KV=1,ISPV(L)
          J=IPVIB(KV,N)
          IF (J <0) THEN
            J=-J
            IF (J == 99999) J=0
          END IF
          TOTEN=TOTEN+DFLOAT(J)*BOLTZ*SPVM(1,KV,L)
        END DO
      END IF
    END IF
  END IF
END IF
IF (M == 8) THEN
  !HF(1)=0.d0         !H2 entalpy of formation [kcal/mol]; from O Conaire et al (2004)
  !HF(2)=52.098d0     !H
  HF(1)=0.d0         !N2 entalpy of formation [kcal/mol]; from O Conaire et al (2004)
  HF(2)=112.954d0    !N
  HF(3)=0.d0         !O2
  HF(4)=59.56d0      !O
  HF(5)=8.91d0       !OH
  HF(6)=-57.77d0     !H2O
  HF(7)=3.d0         !HO2
  HF(8)=-94.05d0     !CO2   !from wikipedia = 393509 J/mol
!  HF(8)=0.d0         !Ar
  HF(:)=HF(:)*4184.d3/AVOG  !converting to [J/molec]
  IF (I == 0) THEN
    DO N=1,NM
      IF (IPCELL(N) > 0) THEN  !energies are normalized by BOLTZ
        TOTENI=TOTEN
        L=IPSP(N)
        TOTEN=TOTEN+HF(L)/BOLTZ
        TOTEN=TOTEN+0.5D00*(SP(5,L)/BOLTZ)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)/BOLTZ
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (IPVIB(K,N) < 0) THEN
              WRITE (*,*) 'Dissociation marked molecule still in flow',N,IPVIB(K,N)
!              STOP
            END IF
            CALL VIB_ENERGY(A,IPVIB(K,N),K,L)
            TOTEN=TOTEN+A/BOLTZ !IPVIB(K,N)*SPVM(1,K,L)
         END DO
        END IF
        IF ((TOTEN-TOTENI) > 1.d-16/BOLTZ) WRITE (*,*) 'MOL',N,' SPECIES',L,' ENERGY/BOLTZ',TOTEN-TOTENI
      END IF
    END DO
!
!    WRITE (9,*) 'Total Energy =',TOTEN,NM
!    WRITE (*,*) 'Total Energy =',TOTEN,NM
  ELSE
    N=I
    IF (IPCELL(N) > 0) THEN
      L=IPSP(N)
      IF (IPCELL(N) > 0) THEN
        L=IPSP(N)
        IF (L==2) TOTEN=TOTEN+3.62D-19
        IF (L==4) TOTEN=TOTEN+4.14D-19
        IF (L==5) TOTEN=TOTEN+0.65D-19
        IF (L==6) TOTEN=TOTEN-4.02D-19
        IF (L==7) TOTEN=TOTEN+0.17D-19
        TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IV=IPVIB(K,N)
            IF (IV < 0) THEN
              IF (IV == -99999) THEN
                IV=0
              ELSE
                IV=-IV
              END IF
              IF (L == 1) TOTEN=TOTEN+7.24E-19
              IF (L == 3) TOTEN=TOTEN+8.28E-19
              IF (L == 5) TOTEN=TOTEN+7.76E-19
              IF (L == 6) TOTEN=TOTEN+8.29E-19
              IF (L == 7) TOTEN=TOTEN+3.45E-19
            END IF
            TOTEN=TOTEN+IV*BOLTZ*SPVM(1,K,L)
          END DO
        END IF
        IF (ABS(TOTEN) > 1.D-16) THEN
          CONTINUE
        END IF
      END IF
    END IF
  END IF
END IF
!
RETURN
!
END SUBROUTINE ENERGY
