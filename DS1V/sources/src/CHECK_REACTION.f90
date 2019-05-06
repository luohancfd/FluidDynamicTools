!
!************************************************************************************
!
SUBROUTINE CHECK_REACTION(IKA,N,L,LM,M,LS,LMS,MS,VRR,VR,VRC,VRI,VCM,RML,RMM,IVDC,LOCAL_NPVIB, LOCAL_EVREM)
!
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
USE OMP_LIB
USE MFDSMC,only:IMF,IMFS,NMFETR,NMFERR,NMFEVR,NMFVTR
!
IMPLICIT NONE
!
!
INTEGER :: I,J,K,L,LM,M,N,LS,LMS,MS,IKA,JR,KV,ISTE(MNRE),NS,NPM,IVDC(MNRE)
REAL(8),EXTERNAL :: GAM
REAL(KIND=8) :: A,ECT,ECR,ECV,EC,VR,VRR,RML,RMM,&
                VDOF1(MMVM),VDOF2(MMVM),EV1(MMVM),EV2(MMVM),ECV1,ECV2,SVDOF1,SVDOF2,ECR1,ECR2,&
                EP,TEMP,ECT2,VRC(3),VCM(3),VRI
INTEGER :: II,IETDX,IERDX(2)

REAL(8) :: LOCAL_EVREM
INTEGER :: LOCAL_NPVIB(3, MMVM,0:100)
LOCAL_NPVIB = 0
LOCAL_EVREM = 0
              !
!--A,B,C working variables
!--J,K,KK,MK,KA,KAA,MKK,JR,KR,JS,KV,IA working integers
!--N the collision cell
!--L,S  the molecule numbers
!--LS,MS,KS species
!--IKA reaction indicator
!--ECR rotational energy
!--ECV vibrational energy
!--ECA available energy
!--STER steric factor
!--ISTE third body species
!--NRE number of reactions that may occur for this pair
!--STERT cumulative steric factor
!--THBCELL number of third body molecules in cell
!--WF weighting factor
!--PSTERT probability
!--VRC relative velocity components: can only be changed by recombination
!--VCM center-of-mass velocity components: can only be changed by recombination
!--VRI magnitude of VRC: can only be changed by recombination
!--VRR new value of u^2+v^2+w^2 containg all the energy
!--VR  new value of relative speed
!--RML,RMM molecule mass parameters
!
NS=ICCELL(3,N)     !sampling cell
TEMP=VAR(10,NS)    !sampling cell vibrational temperature
!
ECT=0.5D00*SPM(1,LS,MS)*VRR
ECR=0.D00; ECV=0.D00
EV1=0.D00; EV2=0.D00
ECR1=0.D00; ECR2=0.D00
ECV1=0.D00; ECV2=0.D00
VDOF1=0.D00; VDOF2=0.D00
SVDOF1=0.D00; SVDOF2=0.D00
IF (ISPR(1,LS) > 0) ECR1=PROT(L)
IF (ISPR(1,MS) > 0) ECR2=PROT(M)
ECR=ECR1+ECR2
IF (MMVM > 0) THEN
  IF (ISPV(LS) > 0) THEN
    DO KV=1,ISPV(LS)
      CALL VIB_ENERGY(EV1(KV),IPVIB(KV,L),KV,LS)
      VDOF1(KV)=2.d0*(SPVM(1,KV,LS)/TEMP)/(DEXP(SPVM(1,KV,LS)/TEMP)-1.d0)
      ECV1=ECV1+EV1(KV)
      SVDOF1=SVDOF1+VDOF1(KV)
    END DO
  END IF
  IF (ISPV(MS) > 0) THEN
    DO KV=1,ISPV(MS)
     CALL VIB_ENERGY(EV2(KV),IPVIB(KV,M),KV,MS)
      VDOF2(KV)=2.d0*(SPVM(1,KV,MS)/TEMP)/(DEXP(SPVM(1,KV,MS)/TEMP)-1.d0)
      ECV2=ECV2+EV2(KV)
      SVDOF2=SVDOF2+VDOF2(KV)
    END DO
  END IF
  ECV=ECV1+ECV2
END IF
!
EC=ECT+ECR+ECV
!
!--check if the molecules follow the order given by IREA
! IF (LS /= MS) THEN
  ! IF (LS == IREA(1,IKA)) IVDC=1  !corrected order
  ! IF (MS == IREA(1,IKA)) IVDC=2  !inverted order
! ELSE
  ! IF (ECV1 >= ECV2) IVDC=1  !choose the high energy molecule for dissociation
  ! IF (ECV2 >  ECV1) IVDC=2
! END IF
!-- Han changed this on Oct 15th 2018
! IVDC(IKA) = 1: LS = IREA(1,IKA), MS = IREA(2,IKA)
! IVDC(IKA) = 2: LS = IREA(2,IKA), MS = IREA(1,IKA)
! For dissociation reaction, IVDC(IKA) tell which molecule is the one to be dissociated
!
!
IF (IKA > 0) THEN
  NPM=NREA(1,IKA)+NREA(2,IKA)
!
!--sample pre-reaction vibrational levels
  DO KV=1,MMVM
    I=0; J=0; K=0
    IF (IVDC(IKA) == 1) THEN
      SELECT CASE(NPM)
        CASE(1) !recombination
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=IPVIB(KV,ISTE(JR)) !third body in recombination
        CASE(2) !exchange
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=0
        CASE(3) !dissociation
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=0
      END SELECT
    ELSE
      SELECT CASE(NPM)
        CASE(1) !recombination
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=IPVIB(KV,ISTE(JR)) !third body in recombination
        CASE(2) !exchange
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=0
        CASE(3) !dissociation
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=0
      END SELECT
    END IF
    IF (I > 100) I=100
    IF (J > 100) J=100
    IF (K > 100) K=100
    LOCAL_NPVIB(1,KV,I) = LOCAL_NPVIB(1,KV,I)+1
    LOCAL_NPVIB(2,KV,J) = LOCAL_NPVIB(2,KV,J)+1
    LOCAL_NPVIB(3,KV,K) = LOCAL_NPVIB(3,KV,K)+1
    ! !$omp atomic
    ! NPVIB(1,IKA,1,KV,I)=NPVIB(1,IKA,1,KV,I)+1  !bin counter
    ! !$omp atomic
    ! NPVIB(1,IKA,2,KV,J)=NPVIB(1,IKA,2,KV,J)+1
    ! !$omp atomic
    ! NPVIB(1,IKA,3,KV,K)=NPVIB(1,IKA,3,KV,K)+1

    !//TODO
    IF (IREAC > 0 .and. NPM == 3 .and. IMF .ne. 0 .and. KV == 1 .and. IMFS == 1 .and. MNRE > 0) THEN
    !$omp critical
      IETDX = FLOOR(ECT/BOLTZ/FTMP0/0.01D0)+1;
      IETDX=MIN(IETDX,1000)
      NMFETR(IETDX,IKA) = NMFETR(IETDX,IKA) + 1.0d0

      ! IREDX(1) is the one dissociated
      IF (IVDC(IKA) == 1) THEN
        IERDX(1) = FLOOR(ECR1/BOLTZ/FTMP0/0.01D0)+1
        IERDX(2) = FLOOR(ECR2/BOLTZ/FTMP0/0.01D0)+1
      ELSE
        IERDX(2) = FLOOR(ECR1/BOLTZ/FTMP0/0.01D0)+1
        IERDX(1) = FLOOR(ECR2/BOLTZ/FTMP0/0.01D0)+1
      END IF

      DO II = 1,2
        IERDX(II) = MIN(IERDX(II),1000)
        NMFERR(IERDX(II),II,IKA) = NMFERR(IERDX(II),II,IKA)+1.0d0
      END DO

      NMFEVR(I,1,IKA) = NMFEVR(I,1,IKA) + 1.0d0
      NMFEVR(J,2,IKA) = NMFEVR(J,2,IKA) + 1.0d0

      NMFVTR(I,IETDX,1,IKA) = NMFVTR(I,IETDX,1,IKA) + 1.0d0
      NMFVTR(J,IETDX,2,IKA) = NMFVTR(J,IETDX,2,IKA) + 1.0d0
    !$omp end critical
    ENDIF

    !remove!$omp critical
    IF (NPM == 3 .and. KV == 1 .and. NPM ==3) THEN
      IF (IVDC(IKA) == 1) THEN
        LOCAL_EVREM = LOCAL_EVREM + ECV1
!        EVREM(IKA) = EVREM(IKA) + ECV1  ! EVREM in Joule
      ELSE
        LOCAL_EVREM = LOCAL_EVREM + ECV2
 !       EVREM(IKA) = EVREM(IKA) + ECV2
      END IF
    END IF
    !remove!$omp end critical
  END DO
!
!--sample pre-reaction vibrational energies (normalized by 1000*BOLTZ)
  I=ECV/(1.d3*BOLTZ)         !bin index
  IF (I < 100) THEN
    !$omp atomic
    NEVIB(IKA,1,I)=NEVIB(IKA,1,I)+1 !bin accumulation
  END IF
!
!   -----------------------------------------------------------
!--modify molecule states according to post-reaction conditions
  IF (IREAC <= 1) THEN
    EP=EC+REA(5,IKA)         !post-reaction total energy (not accounting for third body yet)
!
    IF (NPM == 1) THEN
!--one post-collision molecule (recombination)
      IPSP(L)=JREA(1,IKA,1)  !molecule L becomes the recombined molecule
      LS=IPSP(L)
      PV(1:3,L)=VCM(1:3)
      IPCELL(M)=-IPCELL(M)   !molecule M is marked for removal
      M=ISTE(JR)             !calculate collision of molecule L with third body
      MS=IPSP(M)
      VRC(1:3)=PV(1:3,L)-PV(1:3,M)
      VRR=VRC(1)**2+VRC(2)**2+VRC(3)**2
      VRI=DSQRT(VRR)
      ECT2=0.5D00*SPM(1,LS,MS)*VRR
      ECR2=0.d0; ECV2=0.d0
      IF (ISPR(1,MS) > 0) ECR2=PROT(M)
      IF (ISPV(MS) > 0) THEN
          DO KV=1,ISPV(MS)
            CALL VIB_ENERGY(A,IPVIB(KV,M),KV,MS)
            ECV2=ECV2+A !DFLOAT(IPVIB(KV,M))*BOLTZ*SPVM(1,KV,MS)
          END DO
      END IF
      EP=EP+ECT2+ECR2+ECV2   !add third body Etra, Erot, and Evib
      IVDC(IKA)=1  !update order of the molecules
    END IF
!
    IF (NPM == 2) THEN
!--two post-collision molecules (exchange)
      IPSP(L)=JREA(1,IKA,1)
      IPSP(M)=JREA(2,IKA,1)
      LS=IPSP(L)
      MS=IPSP(M)
      IVDC(IKA)=1  !update order of the molecules
    END IF
!
    IF (NPM == 3) THEN
!--three post-collision molecules (dissociation)
!$omp critical
      IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
      NM=NM+1
      LM=NM                    !new born molecule
!$omp end critical
      IF (IVDC(IKA) == 1) THEN
        IPSP(L)=JREA(1,IKA,1)  !update species
        IPSP(LM)=JREA(1,IKA,2)
        LS=IPSP(L)
        LMS=IPSP(LM)
        PX(LM)=PX(L)
        IPCELL(LM)=IPCELL(L)
        PTIM(LM)=PTIM(L)
      ELSE
        IPSP(M)=JREA(1,IKA,1)
        IPSP(LM)=JREA(1,IKA,2)
        MS=IPSP(M)
        LMS=IPSP(LM)
        PX(LM)=PX(M)
        IPCELL(LM)=IPCELL(M)
        PTIM(LM)=PTIM(M)
      END IF
    PROT(LM)=0.d0
    IPVIB(0:MMVM,LM)=0
    END IF
!
!--set all internal energies to zero
    PROT(L)=0.d0
    PROT(M)=0.d0
    IPVIB(0:MMVM,L)=0
    IPVIB(0:MMVM,M)=0
!
!--update collision pair data
    VRR=2.D00*EP/SPM(1,LS,MS)
    VR=DSQRT(VRR)
    IF (NPM < 3) THEN
      RML=SPM(1,LS,MS)/SP(5,MS)
      RMM=SPM(1,LS,MS)/SP(5,LS)
    ELSE
      IF (IVDC(IKA) == 1) THEN
        RML=SPM(1,IREA(1,IKA),MS)/SP(5,MS)
        RMM=SPM(1,IREA(1,IKA),MS)/SP(5,IREA(1,IKA))
      ELSE
        RML=SPM(1,LS,IREA(1,IKA))/SP(5,IREA(1,IKA))
        RMM=SPM(1,LS,IREA(1,IKA))/SP(5,LS)
      END IF
    END IF
    IF (NPM == 1) THEN       !update VCM for recombination reactions
      VCM(1:3)=RML*PV(1:3,L)+RMM*PV(1:3,M)
    END IF
  END IF
!   -----------------------------------------------------------
END IF
!
RETURN
END SUBROUTINE CHECK_REACTION
