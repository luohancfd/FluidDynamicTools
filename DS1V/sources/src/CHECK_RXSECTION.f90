!
!************************************************************************************
!
SUBROUTINE CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,ECT,SXSECTION,BMAX,BR,IVDC,IDT)
!
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
USE OMP_LIB
USE MFDSMC,only : IMF, IMFdia, NMFANG,&
  MF_SAMPLE_PHASE,MF_CALC_COLL,MF_SAMPLE_ANGLE,MF_EVAL_F
USE EXPCOL,only : EXPCOL_SOLVERM
!
IMPLICIT NONE
!
!
REAL(8),INTENT(IN) :: ECT,SXSECTION
INTEGER,INTENT(IN) :: N,L,M,LS,MS

INTEGER :: J,K,NRE,KK,KS,MK,KA,KAA,MKK,KV,IA,ISTE(MNRE),&
           IDT,NS,NPM,I,IVDC(MNRE),II,MTYPE,IVIB,JROT,IMAX,JMAX
REAL(8),EXTERNAL :: GAM
REAL(KIND=8) :: A,B,ECR,ECV,EC,THBCELL,RANF,&
                VDOF1(MMVM),VDOF2(MMVM),EV1(MMVM),EV2(MMVM),EV,ECV1,ECV2,SVDOF1,SVDOF2,ECR1,ECR2,AL,&
                ATDOF,AIDOF,ANDOF,AVDOF,X,C,D,AA,BB,XI,PSI,PHI,CV,EE,SEV(2),TEMP,&
                RXSECTION(MNRE),STER(MNRE),ECA(MNRE),CF,&
                ALPHA1,ALPHA2,EA,CC,DD
REAL(KIND=8),ALLOCATABLE :: VEC_I(:),VEC_J(:),VALUES(:),ARRAY_TEMP(:,:)
REAL(KIND=8) :: MFANG(8),MFV1,MFV2,MFR1,MFR2,MFCtheta
REAL(KIND=8) :: MFF(2),MFcoll(2)
! variable for IMF=4
REAL(KIND=8) :: BMAX,BR,ECT_E,RM

INTEGER      :: NMFCALL, viblevel(2), IINPM
LOGICAL      :: ISAME
!
!--A,B,C working variables
!--J,K,KK,MK,KA,KAA,MKK,JR,KR,JS,KV,IA working integers
!--N the collision cell
!--L,S  the molecule numbers
!--LS,MS,KS species
!--ECR rotational energy
!--ECV vibrational energy
!--ECA available energy
!--ECT_E effective translational energy taking impact parameter into consideration
!--XSECTION cross-section
!--ISTE third body species
!--NRE number of reactions that may occur for this pair
!--TXSECTION total cross-section
!--THBCELL number of third body molecules in cell
!--WF weighting factor
!--VCM center-of-mass velocity components
!--VRR square of rel speed
!--RML,RMM molecule mass parameters
!--MFANG angles sampled by MF-MC model
!--MFF threshold function maximum 4 kinds of configuration
!--SXSECTION total collision cross sections
!
NMFCALL = 0           !count of MFmodel using
NS=ICCELL(3,N)        !sampling cell
TEMP=VAR(10,NS)       !sampling cell vibrational temperature
NRE=NRSP(LS,MS)       !number of reaction involving species LS and MS
ISAME = .false.       !check if two species are the same
IF (LS == MS) ISAME = .true.
RM = -1.0d0
!
IF (NRE > 0) THEN
!
  MFF = ECT*2.0d0  !initialization MF model
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
  EC=ECT+ECR+ECV  !total collision energy
  RXSECTION=0.d0  !initializing reaction cross-sections
  STER=0.d0       !initializing steric factors used in TCE model
!
  AL=SPM(3,LS,MS)-0.5d0
  ATDOF=(ISPR(1,LS)+ISPR(1,MS)+SVDOF1+SVDOF2+(4.d0-2.d0*AL))*0.5d0 !average number of dofs
  AIDOF=(ISPR(1,LS)+ISPR(1,MS)+SVDOF1+SVDOF2)*0.5d0                !average number of internal dofs
!
  DO J=1,NRE
    K=IRCD(J,LS,MS)
    NPM=NREA(1,K)+NREA(2,K)  !number of post-reaction species
    X=REA(4,K)-0.5d0+AL
    ECA(J)=EC
    CF=1.d0
!
!--collision energy is above reaction energy barrier
    IF ((EC >= REA(2,K)).AND.(EC+REA(5,K) > 0)) THEN !2nd condition prevents a negative post-reaction total energy
!
!--check if molecules follow the order given by IREA
      IVDC(K) = 1
      IF (ISAME) THEN
        IF (ISPV(LS) > 0)THEN
         ! choose the molecule with higher dissociation energy
         ! some model like MFDSMC may change this
          IF (ECV1 .ge. ECV2) THEN
            IVDC(K) = 1
          ELSE
            IVDC(K) = 2
          END IF
        END IF
      ELSE
        IF (LS == IREA(1,K)) IVDC(K) = 1  !corrected order
        IF (MS == IREA(1,K)) IVDC(K) = 2  !inverted order
      END IF
!
!--check reaction model
      MTYPE=0 !TCE/MF/VDC is the standard model
      IF (QCTMODEL==3.AND.GASCODE==8)THEN
        !Enable QCT reaction model
        IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) MTYPE=1 !N2-O collision
        IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) MTYPE=3 !O2-O collision
      END IF
!
!--------------------------------------------------------------------------------
      IF (MTYPE == 0) THEN !use TCE/MFDSMC/VDC model
!
!--possible recombination reaction
        IF (NPM == 1) THEN
          MKK=0
          CALL ZGF(RANF,IDT)
          IF (IRECOM == 1) THEN! .AND. RANF < 0.001d0) THEN  !only one in 1000 possible recombinations are considered and
                                                             !the probability is increased by a factor of 1000
            IF (JREA(2,K,1) < 0) THEN
!--third-body molecule may be any species in cell
              IF (ICCELL(2,N) >= 3) THEN
                MK=L
                DO WHILE ((MK == L).OR.(MK == M))
                  CALL ZGF(RANF,IDT)
                  KK=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                  MK=ICREF(KK)
                END DO
                KS=IPSP(MK)
                ISTE(J)=MK
                A=REA(6,K)
                IA=-JREA(2,K,1)
                A=A*THBP(IA,KS)
              ELSE
                A=0.D00
              END IF
              THBCELL=ICCELL(2,N)
            ELSE
              KS=JREA(2,K,1)
!--the third-body molecule must be species KS
              THBCELL=0.
              MKK=0
              DO KAA=1,ICCELL(2,N)
                KA=ICCELL(1,N)+KAA
                MK=ICREF(KA)
                IF ((IPSP(MK) == KS).AND.(IPCELL(MK) > 0).AND.(IPVIB(0,MK) >= 0)) THEN
                  THBCELL=THBCELL+1.
                  IF ((MK /= L).AND.(MK /= M)) MKK=MK
                END IF
              END DO
              IF (MKK > 0) THEN
                MK=MKK
                ISTE(J)=MK
                A=REA(6,K)
              ELSE
                A=0.
              END IF
            END IF
            B=THBCELL*FNUM/CCELL(1,N)  !number density of third body species in collision cell
            CC=3.5d0-AL !3.5d0-AIDOF-AL
            DD=3.5d0-AL+X !3.5d0-AIDOF-AL+X
            C=GAM(CC)
            D=GAM(DD)
            STER(J)=CF*B*A*(C/D)*ECA(J)**X  !*1000.d0  !note the factor
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          END IF
        END IF
!
!--possible exchange reaction
        IF (NPM == 2) THEN
          CC=ATDOF
          DD=AIDOF+1.5d0+REA(4,K)
          C=GAM(CC)
          D=GAM(DD)
          AA=X+ATDOF-1.d0
          BB=ATDOF-1.d0
          !--note measures to prevent underflows
          STER(J)=CF*(REA(6,K)*(C/D))*(((ECA(J)-REA(2,K))*1.d10)**AA/(ECA(J)*1.d10)**BB)*1.d10**(BB-AA)
          RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
        END IF
!
!--possible dissociation reaction
        IF (NPM == 3) THEN
          IF(REA(1,K) < 0.) THEN
            IF (IMF == 0 .or. (IMF .ne. 0 .and. NMFANG(K) == 0))THEN
              !TCE model (same as in exchange reactions)
              CC=ATDOF
              DD=AIDOF+1.5d0+REA(4,K)
              C=GAM(CC)
              D=GAM(DD)
              AA=X+ATDOF-1.d0
              BB=ATDOF-1.d0
              !--testing nonequilibrium reaction rate corrections
              IF (IVDC(K) == 1) KS=LS
              IF (IVDC(K) == 2) KS=MS
              !
              !--Kuznetsov (KSS)
              !E=(1.d0-DEXP(-SPVM(1,1,KS)/VAR(10,NS)))/(1.d0-DEXP(-SPVM(1,1,KS)/VAR(8,NS)))
              !A=1.d0/(BOLTZ*VAR(10,NS))
              !B=1.d0/(BOLTZ*VAR(8,NS))
              !CF=E*DEXP(-0.7d0*REA(2,K)*(A-B)) !note that the tecplot exact curves are considering a 0.78 factor
              !
              !--testing modified (0.78) Macharet-Fridam correction
              !IF (LS==3.AND.MS==3) THEN !TCE+MF only for O2-O2 collisions
              !  E=(1.d0-DEXP(-SPVM(1,1,KS)/VAR(10,NS)))/(1.d0-DEXP(-SPVM(1,1,KS)/VAR(8,NS)))
              !  A=1.d0/(BOLTZ*VAR(10,NS))
              !  B=1.d0/(BOLTZ*VAR(8,NS))
              !  ALP=(1.d0/2.d0)**2.d0                  !alpha; for homonuclear collisions
              !  TA=ALP*VAR(10,NS)+(1-ALP)*VAR(8,NS)    !T_a
              !  IF (ISPV(LS) >= 1 .AND. ISPV(MS) >= 1) THEN
              !    LF=2.d0*(1-ALP)/(PI**2.d0*ALP**0.75d0) * (BOLTZ*VAR(8,NS)/REA(2,K))**(1.5d0-REA(4,K)) &
              !                                           * (1+3.5*(1-ALP)*(1+DSQRT(ALP))*(BOLTZ*VAR(8,NS)/REA(2,K)))
              !  ELSE
              !    LF=(9.d0/64.d0)*DSQRT(PI*(1-ALP)) * (BOLTZ*VAR(8,NS)/REA(2,K))**(1.d0-REA(4,K)) &
              !                                      * (1+2.5*(1-ALP)*(BOLTZ*VAR(8,NS)/REA(2,K)))
              !  END IF
              !  PHI2=1.373274d6*VAR(8,NS)**(-1.309734d0)
              !  CF=E*(1-LF)*DEXP(-0.78d0*REA(2,K)*(A-B)) + PHI2*LF*DEXP(-0.78d0*REA(2,K)*(1.d0/(BOLTZ*TA)-B))
              !  IF (CF > 1.d0) CF=1.d0 !to correct only vibrationally cold conditions
              !END IF
              !
              !--note measures to prevent underflows
              STER(J)=CF*(REA(6,K)*(C/D))*(((ECA(J)-REA(2,K))*1.d10)**AA/(ECA(J)*1.d10)**BB)*1.d10**(BB-AA)
            ELSE
              IF (IMF == 4) THEN
                IF (RM < 0) THEN
                  IF (INONVHS(LS,MS) == 2) THEN
                    RM = EXPCOL_SOLVERM(LS,MS, BMAX*DSQRT(BR),ECT/EVOLT)
                  ELSE
                    RM = BMAX
                  END IF
                END IF
                ECT_E = ECT * (1-BR*BMAX**2/RM**2)
              ELSE
                ECT_E = ECT
              END IF

              ! Macheret-Fridman-MC mode
              !
              !---get energy of each one particle
              ! here I assume that if the molecule doesn't have vibrational energy
              ! IPVIB will be 0
              IF (IVDC(K) == 1)THEN
                MFV1 = ECV1; MFR1 = ECR1
                MFV2 = ECV2; MFR2 = ECR2
                viblevel(1) = IPVIB(1,L); viblevel(2) = IPVIB(1,M)
              ELSE
                MFV1 = ECV2; MFR1 = ECR2
                MFV2 = ECV1; MFR2 = ECR1
                viblevel(1) = IPVIB(1,M); viblevel(2) = IPVIB(1,L)
              END IF
              !
              NMFCALL = NMFCALL + 1
              IF (NMFCALL > 2) THEN
                WRITE (*,"(A,I3,'+',I3)") "There are more than 2 MF-diss for ",LS,MS
                STOP
              END IF
              !---set mass of mb (MFcol(1)) and mc (MFcol(2))
              CALL MF_CALC_COLL(K,NMFCALL,MFcoll,IDT)
              !
              !---generate geometries of collision
              !
              ! atom-diatom: gamma1, gamma2, theta, phi
              ! diatom-diatom: gamma1, gamma2, theta, phi0, phi1, beta1, beta2, delta, beta
              ! MFANG:
              !  1. gamma_1
              !  2. gamma_2
              !  3. gamma_
              !  4. phi_0
              !  5. phi_1
              !  6. beta_1
              !  7. beta_2
              !  8. delta
              CALL MF_SAMPLE_ANGLE(K,NMFCALL,MFANG,MFCtheta,viblevel,MFV1,MFV2,IDT)
                ! log angles for binning
                !IF (IREAC == 2)THEN
                  !DO II = 1,NMFANG(K)
                    !JJ = FLOOR(MFANG(II)*180.0d0)+1
                    !JJ = MIN(JJ,180)
                    !!              NCANGLE(II,JJ) = NCANGLE(II,JJ)+1.d0
                  !END DO
                !END IF
              CALL MF_EVAL_F(K,NMFCALL,MFV1,MFV2,MFR1,MFR2,MFcoll,MFANG,MFCtheta,MFF(NMFCALL),IDT)
              STER(J) = 0.0d0
              IF (ECT_E >= MFF(NMFCALL)) STER(J) = 1.0D0

              IF (ISAME .and. IMFdia == 1) THEN
                ! check the possibility of dissociate the molecule with lower vibrational energy
                NMFCALL = NMFCALL + 1
                ! exchange of viblevel
                ! this is actually not needed, but we have it here to keep consistency
                II = viblevel(1)
                viblevel(1) = viblevel(2)
                viblevel(2) = II
                ! exchange MFcoll
                CALL MF_CALC_COLL(K,NMFCALL,MFcoll,IDT)
                ! exchange some angles
                CALL MF_SAMPLE_ANGLE(K,NMFCALL,MFANG,MFCtheta,viblevel,MFV2,MFV1,IDT)
                ! exchange internal energy
                CALL MF_EVAL_F(K,NMFCALL,MFV2,MFV1,MFR2,MFR1,MFcoll,MFANG,MFCtheta,MFF(NMFCALL),IDT)

                IF (MFF(2) < MFF(1) .and. MFF(2) <= ECT_E) THEN
                  STER(J) = 1.0d00
                  IVDC(K) = 3-IVDC(K)
                  ! the configuration IVDC(K) is prefered
                END IF
              ELSE
                IF ((NMFCALL == 2) .and. (STER(J) > 0.9d0) )THEN
                  IF (MFF(2) < MFF(1))THEN
                    DO II=1,J-1
                      IINPM = NREA(1,IRCD(II,LS,MS))+NREA(2,IRCD(II,LS,MS))
                      IF (IINPM == 3 .and. NMFANG(II) > 0) then
                        STER(II) =0.0d0 ! set previous dissociation to 0
                        RXSECTION(II) = 0.0d0
                      END IF
                    END DO
                  ELSE
                    STER(J) = 0.0d0
                  END IF
                END IF
              !ELSE  ! MF probability version
                !STER(J)=0.0D0
                !IF (MFV1 .LE. 1.0D-5) MFV1 = 0.5d0*BOLTZ*SPVM(1,1,IREA(1,K))
                !C = MFALPHA(K)*MFDSTAR
                !CC = DSQRT(MFV1 / MFDSTAR)
                !DD = DSQRT(MFALPHA(K))
                !IF (MFV1 .LT. C) THEN
                  !MFF = MFDSTAR/(1.0d0-MFALPHA(K))*(1.0d0-DD*CC)**2
                  !IF ( ECT > MFF) THEN
                    !STER(J) = 4.0D0*(ECT-MFF)**1.5d0 /(3.0d0*PI*PI*MFF*DSQRT(MFDSTAR/(1.0d0-MFALPHA(K))))
                    !STER(J) = STER(J) / DSQRT((1.0d0-CC*DD)*(1.0d0-(2.0d0-DD)*CC))
                  !END IF
                !ELSE
                  !MFF = MFDSTAR - MFV1
                  !IF ( ECT > MFF) THEN
                    !STER(J) = (ECT-MFF)**2/(2.0d0*PI*PI*MFF*MFDSTAR*DD/DSQRT((1.0d0+DD)/(1.0d0-DD)))
                    !STER(J) = STER(J) / DSQRT((1.0d0+MFV1/MFDSTAR)*(MFV1/MFDSTAR/MFALPHA(K)-1.0D0))
                  !END IF
                !END IF
                !IF ( DABS(MFV1-MFDSTAR) < 1.0D0 ) STER(J) = 1.0D0
                !IF (STER(J) > 0.999D0)  STER(J) = 1.0D0
              !END IF
              END IF
            END IF
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          ELSE
            !VDC model
            IF (IVDC(K) == 1) AVDOF=SVDOF1*0.5d0    !internal dofs contributing to vibrational favoring
            IF (IVDC(K) == 2) AVDOF=SVDOF2*0.5d0
            ANDOF=ATDOF-AVDOF                    !internal dofs not contributing to vibrational favoring
            XI=ANDOF-1.d0
            PSI=REA(4,K)+ATDOF+0.5d0-(4.d0-2.d0*AL)*0.5d0
            PHI=REA(1,K)
            CC=ANDOF
            DD=PSI+1.d0
            C=GAM(CC)
            D=GAM(DD)
            CV=1.d0
            BB=1.d0
            EE=1.d0
            EV=0.d0
            IF (MMVM > 0) THEN
              IF (IVDC(K) == 1) THEN
                DO KV=1,ISPV(LS)
                  AA=(TEMP/SPVM(1,KV,LS))**(VDOF1(KV)*0.5d0)
                  SEV=0.d0
                  DO I=1,IDL(KV,LS)
                    SEV(1)=SEV(1)+(DFLOAT(I)*SPVM(1,KV,LS))**PHI
                    SEV(2)=SEV(2)+DEXP(-DFLOAT(I)*SPVM(1,KV,LS)/TEMP)
                  END DO
                  CV=CV*AA*SEV(1)/SEV(2)
                  BB=BB*SPVM(1,KV,LS)**(VDOF1(KV)*0.5d0)
                  EE=EE*(EV1(KV)/BOLTZ)**PHI
                END DO
                EV=ECV1
              ELSE
                DO KV=1,ISPV(MS)
                  AA=(TEMP/SPVM(1,KV,MS))**(VDOF2(KV)*0.5d0)
                  SEV=0.d0
                  DO I=1,IDL(KV,MS)
                    SEV(1)=SEV(1)+(DFLOAT(I)*SPVM(1,KV,MS))**PHI
                    SEV(2)=SEV(2)+DEXP(-DFLOAT(I)*SPVM(1,KV,MS)/TEMP)
                  END DO
                  CV=CV*AA*SEV(1)/SEV(2)
                  BB=BB*SPVM(1,KV,MS)**(VDOF2(KV)*0.5d0)
                  EE=EE*(EV2(KV)/BOLTZ)**PHI
                END DO
                EV=ECV2
              END IF
            END IF
            STER(J)=CF*(REA(6,K)*(C/D)/(CV*BB))*(((ECA(J)-REA(2,K))/BOLTZ)**PSI/((ECA(J)-EV)/BOLTZ)**XI)*EE
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          END IF
        END IF
!
      ELSE
!--------------------------------------------------------------------------------
!
!--possible exchange reaction with QCT model
        IF (NPM == 2) THEN
          IF (MTYPE==1) THEN
            ! data only for N2+O=>NO+N reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 50) IVIB=50
!
            IMAX=60
            ALLOCATE (VEC_I(0:IMAX))
            VEC_I=(/2.454482d-4,2.427769d-4,2.401402d-4,2.375291d-4,2.348600d-4,2.322024d-4,&
                    2.296469d-4,2.270147d-4,2.243892d-4,2.217688d-4,2.191515d-4,2.166443d-4,&
                    2.140356d-4,2.114257d-4,2.088131d-4,2.061956d-4,2.036905d-4,2.010655d-4,&
                    1.984311d-4,1.957851d-4,1.932510d-4,1.905840d-4,1.878997d-4,1.851937d-4,&
                    1.826079d-4,1.798666d-4,1.772361d-4,1.745758d-4,1.717679d-4,1.689237d-4,&
                    1.661840d-4,1.632698d-4,1.604614d-4,1.576103d-4,1.547141d-4,1.516307d-4,&
                    1.486383d-4,1.455892d-4,1.424788d-4,1.393016d-4,1.360511d-4,1.327185d-4,&
                    1.294563d-4,1.259450d-4,1.224778d-4,1.187449d-4,1.150340d-4,1.111769d-4,&
                    1.071570d-4,1.029527d-4,9.853521d-5,9.402503d-5,8.922437d-5,8.408293d-5,&
                    7.852777d-5,7.259131d-5,6.598628d-5,5.858944d-5,4.994199d-5,3.901624d-5,&
                    1.683039d-5/)
            B=VEC_I(IVIB)  !ro-vibrational energy ladder alpha_1 fitting parameter for N2
            JROT=-0.5d0+DSQRT(0.25d0+(ECR/EVOLT)/B)
            B=2.88d0    !characteristic rot temperature for N2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 0) JROT=0
            IF (JROT > 150) JROT=150
            DEALLOCATE (VEC_I)
!
            IMAX=5
            JMAX=10
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,20,50,100,150/)         !rot levels
            VEC_J=(/0,1,5,7,10,15,20,25,30,50/) !vib levels
!
            VALUES=(/ 106.6040, 104.3597, 108.0698, 351.2132,  0.    ,&
                       87.8864,  70.9021, 197.9063, 359.9827, 46.6272,&
                      266.9474, 259.8274, 206.4693, 463.3182, 53.9792,&
                      417.9959, 385.7538, 305.3864, 409.3383, 48.7572,&
                      800.0000, 672.3581, 240.2589, 196.8796, 42.5553,&
                      800.0000, 622.1744, 260.5081, 160.8333, 37.3285,&
                      321.1994, 269.5515,  58.4387,  84.0821, 36.1900,&
                      213.5262, 100.3110,  44.1216,  73.1788, 21.1884,&
                      107.6799,  76.4448,  93.3640,  88.2306, 17.3448,&
                      80.4866,   94.9898,  71.1551,   0.    ,  0.     /)
            VALUES=VALUES/3.d0
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),A)
!
            VALUES=(/-2.1495, -2.1276, -2.0766, -2.5219, -2.3018,&
                     -2.0134, -1.8961, -2.3555, -2.5117, -2.3504,&
                     -2.4319, -2.4155, -2.2638, -2.7787, -2.5659,&
                     -2.6048, -2.5582, -2.4071, -2.8969, -2.5814,&
                     -2.8641, -2.7715, -2.2370, -2.6994, -2.6589,&
                     -3.1547, -3.0446, -2.6986, -2.9661, -2.8068,&
                     -3.0753, -2.9859, -1.9479, -2.7587, -3.1888,&
                     -3.2309, -2.5819, -1.9158, -3.0487, -2.9818,&
                     -3.0438, -2.6741, -3.1186, -3.9082, -2.9116,&
                     -6.6570, -7.7534, -9.3717,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA1)
!
            VALUES=(/ 3.2628, 3.2630, 3.2628, 3.2629, 6.4555,&
                      3.2628, 3.2637, 3.2628, 3.2632, 6.5782,&
                      3.2629, 3.2628, 3.2629, 3.7551, 6.9618,&
                      3.2628, 3.2635, 3.2629, 4.2386, 7.2907,&
                      3.2628, 3.2629, 3.3366, 4.9327, 7.8277,&
                      3.9871, 4.0736, 4.5092, 6.0078, 8.5064,&
                      5.0910, 5.1721, 5.5801, 6.9790, 9.1292,&
                      6.0918, 6.1672, 6.5467, 7.8421, 9.9821,&
                      6.9871, 7.0568, 7.4069, 8.5925,10.2970,&
                      9.4269, 9.4681, 9.6689, 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),EA)
!
            D=9.8216d0                                                              !dissociation energy in eV
            ALPHA2=ALPHA1*(1.d0-D/EA)                                               !2nd exponential coefficient
            IF ((EC/EVOLT) <= EA) A=0.d0                                            !to satisfy the fitting constraint
            RXSECTION(J)=A*(((EC/EVOLT)/EA)**ALPHA1)*((1.d0-EA/(EC/EVOLT))**ALPHA2) !SIGMA_EXCHANGE in Angstrons^2
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
        END IF
!
!--possible dissociation reaction with QCT model
        IF (NPM == 3) THEN
          IF (MTYPE==1) THEN
            ! data only for N2+O=>2N+O reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 55) IVIB=55
!
            IMAX=60
            ALLOCATE (VEC_I(0:IMAX))
            VEC_I=(/2.454482d-4,2.427769d-4,2.401402d-4,2.375291d-4,2.348600d-4,2.322024d-4,&
                    2.296469d-4,2.270147d-4,2.243892d-4,2.217688d-4,2.191515d-4,2.166443d-4,&
                    2.140356d-4,2.114257d-4,2.088131d-4,2.061956d-4,2.036905d-4,2.010655d-4,&
                    1.984311d-4,1.957851d-4,1.932510d-4,1.905840d-4,1.878997d-4,1.851937d-4,&
                    1.826079d-4,1.798666d-4,1.772361d-4,1.745758d-4,1.717679d-4,1.689237d-4,&
                    1.661840d-4,1.632698d-4,1.604614d-4,1.576103d-4,1.547141d-4,1.516307d-4,&
                    1.486383d-4,1.455892d-4,1.424788d-4,1.393016d-4,1.360511d-4,1.327185d-4,&
                    1.294563d-4,1.259450d-4,1.224778d-4,1.187449d-4,1.150340d-4,1.111769d-4,&
                    1.071570d-4,1.029527d-4,9.853521d-5,9.402503d-5,8.922437d-5,8.408293d-5,&
                    7.852777d-5,7.259131d-5,6.598628d-5,5.858944d-5,4.994199d-5,3.901624d-5,&
                    1.683039d-5/)
            B=VEC_I(IVIB)  !ro-vibrational energy ladder alpha_1 fitting parameter for N2
            JROT=-0.5d0+DSQRT(0.25d0+(ECR/EVOLT)/B)
            B=2.88d0    !characteristic rot temperature for N2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 0) JROT=0
            IF (JROT > 150) JROT=150
            DEALLOCATE (VEC_I)
!
            IMAX=5
            JMAX=11
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,20,50,100,150/)              !rot levels
            VEC_J=(/0,1,5,7,10,15,20,25,30,50,55/) !vib levels
!
            VALUES=(/ 13.6954,  3.8612,  8.8992, 32.0071, 22.7469,&
                      15.9766, 16.6928, 10.1433, 18.1889, 20.5395,&
                      21.7207, 28.8407, 14.1536, 43.3566, 33.4931,&
                      17.5769, 15.7801, 27.7189, 29.9199, 45.3146,&
                      21.7649, 23.5739, 29.0426, 33.5678, 50.3760,&
                      46.4568, 38.8535, 47.1537, 44.6373, 66.0602,&
                      62.2184, 53.2935, 55.2718, 71.2303, 67.1038,&
                      58.4835, 67.2176, 59.4007, 65.7182, 76.9170,&
                      76.7714, 70.4267, 46.9692, 55.4398, 62.7420,&
                      52.5293, 52.5709, 56.9886,  0.    ,  0.    ,&
                      34.5883, 33.3195, 40.4456,  0.    ,  0.     /)
            VALUES=VALUES/3.d0
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),A)
!
            VALUES=(/-0.8562,  0.1535, -0.4418, -1.1983, -0.7820,&
                     -0.9275, -0.9514, -0.5159, -0.7844, -0.6867,&
                     -0.9779, -1.1816, -0.6282, -1.3260, -0.9121,&
                     -0.7563, -0.6698, -1.0632, -1.0196, -1.0701,&
                     -0.8428, -0.9005, -1.0460, -1.0323, -1.0748,&
                     -1.2716, -1.1245, -1.2475, -1.0968, -1.1749,&
                     -1.3609, -1.2548, -1.2502, -1.3362, -1.1218,&
                     -1.2437, -1.3355, -1.2301, -1.2168, -1.1619,&
                     -1.3797, -1.3113, -1.0134, -1.0469, -0.9061,&
                     -0.9389, -0.9274, -0.9389,  0.    ,  0.    ,&
                     -0.4158, -0.3554, -0.2485,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA1)
!
            VALUES=(/ 2.7421, 1.9663, 2.4130, 3.1189, 2.4395,&
                      2.8119, 2.8250, 2.4530, 2.6669, 2.3048,&
                      2.8314, 3.0095, 2.4938, 3.0070, 2.4202,&
                      2.6106, 2.5168, 2.8200, 2.5957, 2.5188,&
                      2.5600, 2.5917, 2.6326, 2.4992, 2.4386,&
                      2.7926, 2.6603, 2.7038, 2.4336, 2.3362,&
                      2.7358, 2.6093, 2.5622, 2.4924, 2.0266,&
                      2.4224, 2.5002, 2.3485, 2.1483, 1.7466,&
                      2.3265, 2.2610, 1.9071, 1.7467, 1.3238,&
                      0.9529, 0.9338, 0.8428, 0.    , 0.    ,&
                      0.3998, 0.3449, 0.1476, 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA2)
!
            D=9.8216d0                                                            !dissociation energy in eV
            RXSECTION(J)=A*(((EC/EVOLT)/D)**ALPHA1)*((1.d0-D/(EC/EVOLT))**ALPHA2) !SIGMA_DISS in Angstrons^2
!    IF (TEMP < 6.d3)  RXSECTION(J)=RXSECTION(J)*1.d3                      !trick to sample low T rates
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
!
          IF (MTYPE==3) THEN
            ! data only for O2+O=>3O reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 40) IVIB=40
!
            B=2.07d0    !characteristic rot temperature for O2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 11) JROT=11
            IF (JROT > 231) JROT=231
!
            IMAX=5
            JMAX=6
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,10,20,30,40/)         !vib levels
            VEC_J=(/11,51,101,151,201,231/) !rot levels
!
            VALUES=(/-0.2146,-1.4838,-0.9152,-0.6273, 0.3143,&
                     -0.4933,-1.2553,-0.6653,-0.8075, 0.    ,&
                     -0.7825,-1.1472,-0.6812, 0.    , 0.    ,&
                      0.3910, 0.1276, 0.    , 0.    , 0.    ,&
                      0.4847, 0.    , 0.    , 0.    , 0.    ,&
                      0.6097, 0.    , 0.    , 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),ALPHA1)
!
            VALUES=(/ 2.7549, 3.0977, 2.2660, 1.4933, 0.2179,&
                      2.8395, 2.8527, 1.9294, 1.4879, 0.    ,&
                      2.6707, 2.4818, 1.7044, 0.    , 0.    ,&
                      1.9083, 1.8382, 0.    , 0.    , 0.    ,&
                      2.0407, 0.    , 0.    , 0.    , 0.    ,&
                      0.6837, 0.    , 0.    , 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),ALPHA2)
!
            VALUES=(/ 13.3864, 182.0998,115.1628,215.0688, 93.7781,&
                      20.2472, 165.1482, 80.5651,300.0000,  0.    ,&
                      53.4346, 160.0991,252.3550,  0.    ,  0.    ,&
                      28.6315, 198.8150,  0.    ,  0.    ,  0.    ,&
                     106.1007,   0.    ,  0.    ,  0.    ,  0.    ,&
                     128.6496,   0.    ,  0.    ,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),A)
!
            D=5.21275d0                                                           !dissociation energy in eV
            RXSECTION(J)=A*(((EC/EVOLT)/D)**ALPHA1)*((1.d0-D/(EC/EVOLT))**ALPHA2) !SIGMA_DISS in Angstrons^2
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
        END IF
      END IF
!--------------------------------------------------------------------------------
    END IF
!
  END DO
END IF
!
RETURN
END SUBROUTINE CHECK_RXSECTION
