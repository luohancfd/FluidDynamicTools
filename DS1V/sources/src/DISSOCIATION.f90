!
!*****************************************************************************
!
SUBROUTINE DISSOCIATION
!
!--dissociate molecules that have been marked with IPVIB(0,:) < 0
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
!
IMPLICIT NONE
!
INTEGER ::K,KK,L,M,LS,MS,KS,JS,KV,IKA,NC,NS,NSP,MAXLEV,II,IV,IDT=0,ITEST=0 !--isebasti: included IDT
REAL(KIND=8) :: A,B,C,ECT,ECC,EVIB,VRR,VR,RMM,RML,RANF,SVDOF(2),PROB,ERM !--isebasti: included RANF
REAL(KIND=8), DIMENSION(3) :: VRC,VCM,VRCP
!
L=0
DO WHILE (L < NM)
  L=L+1
  LS=IPSP(L)
  IF (ISPV(LS) > 0) THEN
    IF (IPVIB(0,L) < 0) THEN  !dissociate through reaction IKA
      IKA=-IPVIB(0,L)
      ECT=PROT(L)  !--ECT is the energy available for post-reaction states
      DO KV=1,ISPV(LS)
        CALL VIB_ENERGY(EVIB,IPVIB(KV,L),KV,LS)
        ECT=ECT+EVIB !IPVIB(KV,L)*BOLTZ*SPVM(1,KV,LS)
      END DO
      IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
      NM=NM+1
      M=NM
!--update species
      IPSP(L)=JREA(1,IKA,1)
      IPSP(M)=JREA(1,IKA,2)
      LS=IPSP(L)
      MS=IPSP(M)
!--set center of mass velocity as that of molecule
      VCM(1:3)=PV(1:3,L)
!--new born molecule
      PX(M)=PX(L)
      IPCELL(M)=IPCELL(L)
      PTIM(M)=PTIM(L)
!--set any internal mode to the ground state and reset IPVIB(0,L) and IPVIB(0,M) to 0
      PROT(L)=0.d0
      PROT(M)=0.d0
      IPVIB(:,L)=0
      IPVIB(:,M)=0
!--effective vibrational degrees of freedom
      NC=IPCELL(L)        !collision cell
      NS=ICCELL(3,NC)     !sampling cell
      A=VAR(10,NS)        !sampling cell vibrational temperature
      SVDOF(:)=0.d0
      IF (ISPV(LS) > 0) THEN
        DO KV=1,ISPV(LS)
          SVDOF(1)=SVDOF(1)+2.d0 !*(SPVM(1,KV,LS)/A)/(DEXP(SPVM(1,KV,LS)/A)-1.d0)
        END DO
      END IF
      IF (ISPV(MS) > 0) THEN
        DO KV=1,ISPV(MS)
          SVDOF(2)=SVDOF(2)+2.d0 !*(SPVM(1,KV,MS)/A)/(DEXP(SPVM(1,KV,MS)/A)-1.d0)
        END DO
      END IF
IF (ITEST == 1) THEN
!--Larsen-Borgnakke serial vibrational energy redistribution
      DO NSP=1,2
        IF (NSP == 1) THEN
          K=L ; KS=LS ; JS=MS
        ELSE
          K=M ; KS=MS ; JS=LS
        END IF
        IF (ISPV(KS) > 0) THEN
          DO KV=1,ISPV(KS)
            CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
            !EVIB=DFLOAT(IPVIB(KV,K))*BOLTZ*SPVM(1,KV,KS)
            ECC=ECT+EVIB
            CALL VIB_LEVEL(ECC,MAXLEV,KV,KS)
            !MAXLEV=ECC/(BOLTZ*SPVM(1,KV,KS))
            CALL ZGF(RANF,IDT)
            II=0
            DO WHILE (II == 0)
              CALL ZGF(RANF,IDT)
              IV=RANF*(MAXLEV+0.99999999D00)
              IPVIB(KV,K)=IV
              CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
              !EVIB=DFLOAT(IV)*BOLTZ*SPVM(1,KV,KS)
              IF (EVIB < ECC) THEN
                IF (SVDOF(NSP) > 0.d0) SVDOF(NSP)=SVDOF(NSP)-2.d0 !*(SPVM(1,KV,KS)/VAR(10,NN))/(DEXP(SPVM(1,KV,KS)/VAR(10,NN))-1.d0)
                IF (NSP == 1) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,KS)+ISPR(1,JS)+SVDOF(2)+SVDOF(1)
                IF (NSP == 2) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,KS)+ISPR(1,JS)+SVDOF(2)
                PROB=(1.D00-EVIB/ECC)**(A*0.5d0-1.d0)   !--see eqn (5.61) and (5.45)
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) II=1
              END IF
            END DO
            ECT=ECC-EVIB
          END DO
        END IF
      END DO
!-------------------------------------------------
!--Larsen-Borgnakke serial rotational energy redistribution
      DO NSP=1,2
        IF (NSP == 1) THEN
          K=L ; KS=LS ; JS=MS
        ELSE
          K=M ; KS=MS ; JS=LS
        END IF
        IF (ISPR(1,KS) > 0) THEN
          ECC=ECT+PROT(K)
          IF (NSP == 1) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,JS)
          IF (NSP == 2) A=(5.d0-2.d0*SPM(3,KS,JS))
          CALL LBS(ISPR(1,KS)*0.5d0-1.d0,A*0.5d0-1.d0,ERM,IDT)
          PROT(K)=ERM*ECC
          ECT=ECC-PROT(K)
        END IF
      END DO
END IF
!-------------------------------------------------
!--calculate new velocities
      VRR=2.d0*ECT/SPM(1,LS,MS)
      VR=DSQRT(VRR)
      VRC=0.d0                  !initial relative velocity
      RML=SPM(1,LS,MS)/SP(5,MS)
      RMM=SPM(1,LS,MS)/SP(5,LS)
      CALL VHSS(LS,MS,VR,VRC,VRCP,IDT)
      PV(1:3,L)=VCM(1:3)+RMM*VRCP(1:3)
      PV(1:3,M)=VCM(1:3)-RML*VRCP(1:3)
      IPCP(L)=M
      IPCP(M)=L
    END IF
  END IF
END DO
!
RETURN
!
END SUBROUTINE DISSOCIATION
