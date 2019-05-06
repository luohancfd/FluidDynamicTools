!
!************************************************************************************
!
SUBROUTINE COLLISIONS
!
!--calculate the collisions
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
USE MFDSMC,only: NMFET0,NMFER0,NMFEV0,NMFVT0,NMFEV,NMFET,NMFER,NMFVT, &
  IMF, IMFS, IMFpair
!
IMPLICIT NONE
!
INTEGER :: N,NS,NN,M,MM,L,LL,K,KK,KT,J,I,II,III,NSP,MAXLEV,IV,IVP,NSEL,KV,LS,MS,KS,JS,IKA,LZ,KL,IS,IREC,&
           NLOOP,IA,IDISS,IE,IEX,JJ,NPRI,LIMLEV,KVV,KW,ISP1,ISP2,ISP3,INIL,INIM,JI,LV,IVM,NMC,NVM,LSI,&
           JX,MOLA,KR,JKV,NSC,KKV,NUMEXR,IAX,NSTEP,ISWITCH,NPM,LM,LMS,IVDC(MNRE),ISHUF(3),IT,IDT=0,IVDC0
REAL(KIND=8) :: A,AA,AAA,AB,B,BB,BBB,ASEL,DTC,SEP,VRI,VR,VRR,ECT,ET0,EVIB,ECC,ZV,COLT,ERM,C,OC,SD,D,CVR,PROB,&
                RML,RMM,ECTOT,ETI,EREC,ET2,XMIN,XMAX,WFC,CENI,CENF,VRRT,TCOLT,EA,DEN,SVDOF(3,3),RANF,DT(ITMAX),& !--isebasti: included RANF,DT
                RXSECTION(MNRE),SXSECTION,TXSECTION,&
                QNU,A1,A2,B1,B2,C1,C2,EROT,EV,SIGMA_REF,EV_POST,SUMF,E,F,EL,ED,EF,S !ME-QCT variables
REAL(KIND=8),DIMENSION(0:100) :: VTXSECTION
REAL(KIND=8),DIMENSION(3) :: VRC,VCM,VRCP,VRCT
REAL(KIND=8) :: ECR(2),EVIBEV,BMAX,REST_DOF
INTEGER :: IETDX,IERDX(2),IEVDX(2),IVPS(2)
logical :: IREACSP
! variable used for nonVHS == 3
LOGICAL :: IVHS
REAL(8) :: CVR2, SXSECTION2
REAL(8),EXTERNAL :: GAM

INTEGER :: LOCAL_NPVIB(3,MMVM,0:100)
REAL(8) :: LOCAL_EVREM

!
!--N,M,K working integer
!--LS,MS,KS,JS molecular species
!--VRC components of the relative velocity
!--RML,RMM molecule mass parameters
!--VCM components of the center of mass velocity
!--VRCP post-collision components of the relative velocity
!--SEP the collision partner separation
!--VRR the square of the relative speed
!--VR the relative speed
!--ECT relative translational energy (update after different process)
!--ET0 relative intitial translational energy in eV
!--EVIB vibrational energy
!--ECC collision energy (rel trans +vib)
!--MAXLEV maximum vibrational level
!--ZV vibration collision number
!--COLT collision temperature (rel trans + vib)
!--TCOLT collision temperature based rel trans
!--SDF the number of degrees of freedom associated with the collision
!--ERM rotational energy
!--NSEL integer number of selections
!--CVR product of collision cross-section and relative speed
!--PROB a probability
!--IKA reaction indicator  ( 0 for no reaction, >0 for reaction )
!--KT third body molecule code
!--ECTOT energy added at recmbination
!--IREC initially 0, becomes 1 if a recombination occurs
!--NPRI the priority collision molecule (the last one rejected as a previous collision molecule)
!--WFC weighting factor in the cell
!--IEX is the reaction that occurs (1 if only one is possible)
!--EA activation energy
!--IVDC order of the two colliding particles
!
! REAL(8) :: MVELC(3), MSVELC(3), VARVELC(3)
!--MVELC mean velocity of the cell, mass averaged
!--MSVELC mean square of the velocity
!--VARIANCE

!IDT: ALWAYS =0 OUTSIDE PARALLEL REGIONS

!$OMP PARALLEL &
!$OMP& PRIVATE(IDT) &
!$OMP& PRIVATE(N,NS,DTC,NN,WFC,AAA,BBB,C,XMIN,XMAX,ASEL,NSEL,I,IE,KL,NPRI,III,K,L,LM,M,RANF,LL) &
!$OMP& PRIVATE(SEP,J,MM,A,KK,VRC,VRR,VRI,VR,CVR,ECT,NSP,KV,EVIB,ECC,MAXLEV,COLT,B,ZV,II,IV,IVP,PROB,SVDOF) &
!$OMP& PRIVATE(ERM,VCM,VRCP,OC,SD,D,LS,LMS,MS,IVDC,IVDC0,ISHUF,RML,RMM,ISWITCH,IDISS,KS,JS,LIMLEV,LZ,IEX,IREC) &
!$OMP& PRIVATE(TCOLT,KT,AA,AB,BB,KVV,JI,LSI,IVM,NMC,NVM,ECTOT,PSF,JJ,EA,DEN,IAX,JX,IKA,NPM,NSTEP,IT,DT) &
!$OMP& PRIVATE(SXSECTION,RXSECTION,VTXSECTION,TXSECTION) &
!$OMP& PRIVATE(QNU,A1,A2,B1,B2,C1,C2,E,F,EL,ED,ET0,EROT,EV,SIGMA_REF,EV_POST,SUMF,EF,S) &
!$OMP& PRIVATE(ECR,EVIBEV,IVPS,BMAX,REST_DOF)&
!$OMP& PRIVATE(CVR2,IVHS, SXSECTION2)&
!$OMP& PRIVATE(IETDX,IEVDX,IERDX,IREACSP,LOCAL_NPVIB,LOCAL_EVREM) &
!$OMP& PRIVATE(ELACK) &
!$OMP& REDUCTION(+:NDISSOC,NDISSL,TRECOMB,NRECOMB,TREACL,TREACG,TNEX,TFOREX,TREVEX) & !Q-K
!$OMP& REDUCTION(+:TOTDUP,TOTCOL,PCOLLS,TCOL,CSCR,COLLS,WCOLLS,CLSEP,REAC,NPVIB, EVREM) &
!$OMP& REDUCTION(+:NMFET0,NMFER0,NMFEV0,NMFVT0,NMFEV,NMFET,NMFER,NMFVT)
!$    IDT=OMP_GET_THREAD_NUM()  !THREAD ID
!$OMP DO SCHEDULE (DYNAMIC)

! OPEN(19, FILE="VelocityVAR.txt", ACCESS="APPEND")
DO N=1,NCCELLS
!
  IF (FTIME-CCELL(5,N) > CCELL(3,N)) THEN
!
!--calculate collisions appropriate to  time DTC
    DTC=2.D00*CCELL(3,N)
    IF (ICCELL(2,N) > 1) THEN  !--no collisions calculated if there are less than two molecules in collision cell
      NN=ICCELL(3,N)
      WFC=1.D00
      IF ((IWF == 1).AND.(IVB == 0)) WFC=1.D00+WFM*CELL(1,NN)**IFX
      CCELL(5,N)=CCELL(5,N)+DTC
      IF (IVB == 0) AAA=CCELL(1,N)
      IF (IVB == 1) THEN
        C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV*NCIS)
        XMIN=XB(1)+DFLOAT(N-1)*C
        XMAX=XMIN+C
        WFC=1.D00+WFM*(0.5D00*(XMIN+XMAX))**IFX
        IF (IFX == 0) AAA=XMAX-XMIN
        IF (IFX == 1) AAA=PI*(XMAX**2-XMIN**2)  !--assumes unit length of full cylinder
        IF (IFX == 2) AAA=1.33333333333333333333D00*PI*(XMAX**3-XMIN**3)    !--flow is in the full sphere
      END IF
!
!--these statements implement the N(N-1) scheme
      ASEL=0.5D00*ICCELL(2,N)*(ICCELL(2,N)-1)*WFC*FNUM*CCELL(4,N)*DTC/AAA+CCELL(2,N)
      NSEL=ASEL
      CCELL(2,N)=ASEL-DFLOAT(NSEL)
!
      IF (NSEL > 0) THEN
!
        ! MVELC = 0.0d0
        ! MSVELC = 0.0d0
        ! DO K = 1, ICCELL(2,N)
          ! M = ICREF(ICCELL(1,N) + K)
          ! LS = ABS(IPSP(M))
          ! DO J = 1,3
            ! MVELC(J) = MVELC(J) + PV(J,M)
            ! MSVELC(J) = MSVELC(J) + PV(J,M)**2
          ! END DO
        ! END DO
        ! DO J=1,3
          ! MVELC(J) = MVELC(J) / ICCELL(2,N)
          ! MSVELC(J) = MSVELC(J) / ICCELL(2,N)
          ! VARVELC(J) = MSVELC(J) - MVELC(J)**2
        ! END DO
        ! WRITE(19,'(I3,1X,20G13.5)') N, (MVELC(J),J=1,3), (MSVELC(J),J=1,3), (VARVELC(J),J=1,3)
!
        NS=0   !--counts the number of selections
        KL=0   !--becomes 1 if it is the last selection
        NPRI=0
!
        DO KL=1,NSEL
          NS=NS+1
          IE=0  !--becomes >0 if an inelastic collision occurs
!
          IF (ICCELL(2,N) == 2) THEN
            K=1+ICCELL(1,N)
            L=ICREF(K)
            K=2+ICCELL(1,N)
            M=ICREF(K)
            IF (M == IPCP(L)) THEN
              III = 1
              CCELL(5,N)=CCELL(5,N)-DTC
            END IF
          ELSE
            K=NPRI
            IF (K > 0) THEN
              L=K
              NPRI=0
            ELSE
              CALL ZGF(RANF,IDT)
              K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
              L=ICREF(K)
            END IF
!--one molecule has been selected at random
            IF (NNC == 0) THEN
!--select the collision partner at random
              M=L
              DO WHILE (M == L)
                CALL ZGF(RANF,IDT)
                K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                M=ICREF(K)
              END DO
            ELSE
!--select the nearest from the total number (< 30) or a random 30
              IF (ICCELL(2,N) < 30) THEN
                LL=ICCELL(2,N)
              ELSE
                LL=30
              END IF
              SEP=1.D10
              M=0
              DO J=1,LL
                IF (LL < 30) THEN
                  K=J+ICCELL(1,N)
                ELSE
                  CALL ZGF(RANF,IDT)
                  K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                END IF
                MM=ICREF(K)
                IF (MM /= L) THEN   !exclude the already selected molecule
                  IF (MM /= IPCP(L)) THEN   !exclude the previous collision partner
                    A=DABS(PX(L)-PX(MM))
                    IF ((A < SEP).AND.(A >1.D-8*DDIV)) THEN
                      M=MM
                      SEP=A
                    END IF
                  ELSE
                    NPRI=MM    !MM is made the priority molecule
                  END IF
                END IF
              END DO
            END IF
          END IF
!
!--Calculate the relative velocity
          DO KK=1,3
            VRC(KK)=PV(KK,L)-PV(KK,M)
          END DO
          VRR=VRC(1)*VRC(1)+VRC(2)*VRC(2)+VRC(3)*VRC(3)
          VR=DSQRT(VRR)
          VRI=VR
!
!--Prepare variables
          LS=ABS(IPSP(L))
          MS=ABS(IPSP(M))
          RML=SPM(1,LS,MS)/SP(5,MS)
          RMM=SPM(1,LS,MS)/SP(5,LS)
          DO KK=1,3
            VCM(KK)=RML*PV(KK,L)+RMM*PV(KK,M)
          END DO
!
          IF (IREAC .EQ. 2) THEN
            DO I = 1,MNRE
              IF ( (LS == IREA(1,I)) .and. (MS == IREA(2,I)) ) THEN
                IREACSP = .true.
                EXIT
              ELSE IF ( (MS == IREA(1,I)) .and. (LS == IREA(2,I))) THEN
                IREACSP = .true.
                EXIT
              ELSE
                IREACSP = .false.
              END IF
            END DO
          ELSE
            IREACSP = .true.
          END IF

!--Simple gas
          IF (MSP == 1) THEN
            CVR=VR*CXSS*((2.D00*BOLTZ*SP(2,1)/(RMAS*VRR))**(SP(3,1)-0.5D00))*RGFS
            IF (CVR > CCELL(4,N)) CCELL(4,N)=CVR
            CALL ZGF(RANF,IDT)
            IF (RANF < CVR/CCELL(4,N)) THEN
!--the collision occurs
              IF ((M == IPCP(L)).AND.(L == IPCP(M))) TOTDUP=TOTDUP+1.D00
              TOTCOL=TOTCOL+1.D00
              TCOL(1,1)=TCOL(1,1)+2.D00
              CSCR(1,1)=CSCR(1,1)+2.d0*CVR
              COLLS(NN)=COLLS(NN)+1.D000
              WCOLLS(NN)=WCOLLS(NN)+WFC
              IF(NN == NSPDF) PCOLLS=PCOLLS+WFC
              SEP=DABS(PX(L)-PX(M))
              CLSEP(NN)=CLSEP(NN)+SEP
              IF ((ISPR(1,1) > 0)) THEN
!--Larsen-Borgnakke serial redistribution
                ECT=0.5D00*RMAS*VRR
                DO NSP=1,2
!--consider the molecules in turn
                  IF (NSP == 1) THEN
                    K=L
                  ELSE
                    K=M
                  END IF
                  IF (MMVM > 0) THEN
                    IF (ISPV(1) > 0) THEN
                      DO KV=1,ISPV(1)
                        CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,1)
                        ECC=ECT+EVIB
                        CALL VIB_LEVEL(ECC,MAXLEV,KV,1)
                        IF (SPVM(3,KV,1) > 0.) THEN
!--note quantizing of COLT in the following statements
                          IF (IZV == 1) THEN
                            COLT=(DFLOAT(MAXLEV)*SPVM(1,KV,1))/(3.5D00-SPM(3,1,1))        !--quantized collision temperature
                          ELSE
                            COLT=VAR(8,NN)        !--isebasti: Ttrans, according to DSMC.f90 (2013)
                          END IF
                          B=SPVM(4,KV,1)/SPVM(3,KV,1)    !Tdiss/Tref
                          A=SPVM(4,KV,1)/COLT            !Tdiss/T
       ZV=(A**SPM(3,1,1))*(SPVM(2,KV,1)*(B**(-SPM(3,1,1))))**(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                        ELSE
                          ZV=SPVM(2,KV,1)
                        END IF
                        CALL ZGF(RANF,IDT)
                        IF (1.D00/ZV > RANF) THEN
                          II=0
                          DO WHILE (II == 0)
                            CALL ZGF(RANF,IDT)
                            IV=RANF*(MAXLEV+0.99999999D00)
                            IPVIB(KV,K)=IV
                            CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,1)
                            IF (EVIB < ECC) THEN
                              CALL INELASTIC_VT(1,1,EVIB,ECC,5.0d0-2.0d0*SPM(3,1,1),PROB)
!--PROB is the probability ratio of eqn (5.61)
                              CALL ZGF(RANF,IDT)
                              IF (PROB > RANF) II=1
                            END IF
                          END DO
                          ECT=ECC-EVIB
                        END IF
                      END DO
                    END IF
                  END IF
!--now rotation of this molecule
                  IF (ISPR(1,1) > 0) THEN
                    IF (ISPR(2,1) == 0) THEN
                      B=1.D00/SPR(1,1)
                    ELSE    !--use molecule rather than mean value
                      COLT=ECC/((2.5D00-SP(3,1))*BOLTZ)
                      B=1.D00/(SPR(1,1)+SPR(2,1)*COLT+SPR(3,1)*COLT*COLT)
                    END IF
                    CALL ZGF(RANF,IDT)
                    IF (B > RANF) THEN
                      ECC=ECT+PROT(K)
                      CALL INELASTIC_RT(1,1,ECC,ERM,IDT)
                      PROT(K)=ERM*ECC
                      ECT=ECC-PROT(K)
                    END IF
                  END IF
                END DO
!--adjust VR for the change in energy
                VR=SQRT(2.D00*ECT/SPM(1,1,1))
              END IF
!--end of L-B redistribution
              DO KK=1,3
                VCM(KK)=0.5d0*(PV(KK,L)+PV(KK,M))
              END DO
!--calculate new velocities
              CALL VHSS(1,1,VR,VRC,VRCP,IDT)
              PV(1:3,L)=VCM(1:3)+0.5d0*VRCP(1:3)
              PV(1:3,M)=VCM(1:3)-0.5d0*VRCP(1:3)
              IPCP(L)=M
              IPCP(M)=L
            END IF   !--collision occurrence
!
!--Gas Mixture
          ELSE IF (IPCELL(L)>0.AND.IPCELL(M)>0.AND.IPVIB(0,L)>=0.AND.IPVIB(0,M)>=0.AND. IREACSP) THEN !trick
            !IPCELL<0 for recombined molecule marked for removal
            !IPVIB(0,L)<0 for molecule marked for dissociation
!
!--Calculate the total scattering (VHS) cross-section
            ECT=0.5D00*SPM(1,LS,MS)*VRR         !collision translational energy in J
            ET0=ECT/EVOLT                        !convert to eV
            CALL CALC_TOTXSEC(LS, MS, VR, VRR, ET0, -1.0d0, SXSECTION, BMAX, CVR)
            IF (nonVHS == 3) THEN
              CVR2 = CVR   !save true vhs cross sections
              SXSECTION2 = SXSECTION

              ! MF-DSMC-Correct
              ! Han added the following to correct reaction rate for MF-DSMC model
              IF ((LS == 3 .and. MS == 4) .or. (LS == 4 .and. MS == 3)) THEN
                ! For O2+O with MF-DSMC model, the total cross sections are corrected to match
                ! O2+O dissociation rates calculated by Marat Kulakhmetov
                IF (IMF == 3) THEN
                  ! MF-DSMC-AHO, Morse version
                  CVR  = CVR*5.8863204392427d0
                  ! original dref of O2+O VHS: 3.4420
                  ! dref_correct= 8.3509
                  ! the multiplier is dref_correct^2/dref^2
                ELSE IF (IMF == 1) THEN
                  ! MF-DSMC-SHO, 2018 Han Luo's PRF
                  CVR  = CVR*6.681561661633024d0
                  ! original dref of O2+O VHS: 3.4420
                  ! dref_correct= 9.6196
                ELSE
                  WRITE(*,*) "Something wrong in collision Code:1"
                  STOP
                END IF
                SXSECTION = CVR/CVR2*SXSECTION2
                IVHS = .false.
              ELSE IF (LS == 3 .and. MS == 3) THEN
                IF (IMF == 1) THEN
                  ! MF-DSMC-SHO model: the following one try to match Park's rate, or Byron's rate
                  ! to use this one, you first need to change dref of O2 to 17.474245e-10,
                  ! and omega to 0.99531
                  !
                  ! the original value of dref should be 4.1515e-10, omega should be 0.7318
                  !
                  ! What I do here is that I first set CVR and SXSECTION as the corrected one with parameter
                  ! dref = 17.474245e-10 and omega = 0.99531
                  !
                  ! And then I set CVR2 and SXSECTION2 as the real VHS one
                  CVR = CVR2*(1.0d0 - DEXP(-2238.0d0/VAR(8,NN)))
                  SXSECTION = SXSECTION2 * CVR / CVR2

                  CVR2 = VR*PI*(4.1515d-10)**2*((2.D00*BOLTZ*SPM(5,LS,MS)/(SPM(1,LS,MS)*VRR))**(0.7318D0-0.5D00))
                  CVR2 = CVR2 / GAM(2.5D0 - 0.7318D0)
                  SXSECTION2 = CVR2/VR*1e20

                  IVHS = .false.
                ELSE IF (IMF == 3) THEN
                  ! MF-DSMC-AHO Morse version
                  ! don't do anything because the model matches well with Ross Chaudry's QCT rate
                  ! dref = 4.1515, omega = 0.7318
                  IVHS = .true.
                ELSE
                  WRITE(*,*) "Something wrong in collision Code:2"
                  STOP
                END IF
              END IF
            ELSE
              IVHS = .true.
            END IF


! get rotational energy and vibrational energy, ECR1 is always the one with lower sp number
! 1 is always the one with lower sp number  or the one with higher Ev
            ECR = 0.0d0; IVPS = 0
            IF (LS < MS ) THEN
              IF (ISPV(LS) > 0) THEN
                ECR(1) = PROT(L); IVPS(1) = IPVIB(1,L)
              END IF
              IF (ISPV(MS) > 0) THEN
                ECR(2) = PROT(M); IVPS(2) = IPVIB(1,M)
              END IF
            ELSE IF (LS > MS) THEN
              IF (ISPV(LS) > 0) THEN
                ECR(2) = PROT(L); IVPS(2) = IPVIB(1,L)
              END IF
              IF (ISPV(MS) > 0) THEN
                ECR(1) = PROT(M); IVPS(1) = IPVIB(1,M)
              END IF
            ELSE ! same specie
              IF (ISPV(LS) > 0) THEN
                IF (IPVIB(1,L) >= IPVIB(1,M)) THEN ! select one with higher Ev
                  ECR(1) = PROT(L); IVPS(1) = IPVIB(1,L)
                  ECR(2) = PROT(M); IVPS(2) = IPVIB(1,M)
                ELSE
                  ECR(1) = PROT(M); IVPS(1) = IPVIB(1,M)
                  ECR(2) = PROT(L); IVPS(2) = IPVIB(1,L)
                END IF
              END IF
            END IF
            ! Sample rotational, translational and vibrational energy
            IF (IREAC .eq. 2 .and. IMF .ne. 0 .and. MNRE > 0 .and. IMFS == 1) THEN
              IF (MNRE > 2) THEN
                WRITE(*,*) " IMFS shouldn't be enabled if MNRE > 2"
              END IF
              IETDX = FLOOR(ECT/BOLTZ/FTMP0/0.010D0)+1
              IETDX = MIN(IETDX,1000)
              NMFET0(IETDX,IMFpair(LS,MS)) = NMFET0(IETDX,IMFpair(LS,MS)) + 1.0D0
              ! KK=1 the one with lower smaller sp number
              DO KK = 1,2
                K = KK
                ! for same specie, X(:,2,:) is always zero
                ! rotational energy
                IERDX(KK) = FLOOR(ECR(KK)/BOLTZ/FTMP0/0.01D00)+1
                IERDX(KK) = MIN(IERDX(KK),1000)
                NMFER0(IERDX(KK),K,IMFpair(LS,MS)) = NMFER0(IERDX(KK),K,IMFpair(LS,MS)) + 1.0D0
                ! vibrational energy
                IEVDX(KK) = IVPS(KK)
                IF (IEVDX(KK) > 100)    IEVDX(KK) = 100
                NMFEV0(IEVDX(KK),K,IMFpair(LS,MS)) = NMFEV0(IEVDX(KK),K,IMFpair(LS,MS)) + 1.0d0
                NMFVT0(IEVDX(KK),IETDX,K,IMFpair(LS,MS)) = NMFVT0(IEVDX(KK),IETDX,K,IMFpair(LS,MS))+1.0D0
              END DO
            END IF

!
!--Calculate the VT (ME-QCT) cross-sections (only for O2+O and N2+O collisions)
            J=0
            IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) J=1 !N2-O collision
            IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) J=3 !O2-O collision
            ! J is only larger than 0 for the above two pairs

            ! Calcuate pre-collision vibrational energy
            ! IVP is the vibrational level of the molecule which is more
            ! likely to be dissociated
            ! For J>0, predefined QCT model *MIGHT* be used later
            ! For J<0, QCT model will never be used
            IVP = -1
            EVIB = 0.0d0
            IF (LS == J) THEN
              IVP=IPVIB(1,L);
            ELSE IF  (MS == J) THEN
              IVP=IPVIB(1,M)
            ELSE IF ((ISPV(LS) == 1) .and. (ISPV(MS) == 0))THEN
              IVP = IPVIB(1,L)
              J = -LS
            ELSE IF ((ISPV(MS) == 1) .and. (ISPV(LS) == 0))THEN
              IVP = IPVIB(1,M)
              J = -MS
            ELSE ! both are diatom, select the one with higher vibrational energy to be dissociated
              IF (IPVIB(1,L) >= IPVIB(1,M)) THEN
                J = -LS
                IVP = IPVIB(1,L)
              ELSE
                J = -MS
                IVP = IPVIB(1,M)
              ENDIF
            END IF
            IF (IVP >= 0)THEN
              CALL VIB_ENERGY(EVIB,IVP,1,ABS(J))
            END IF
            ! the following value of ECC is only used in ME-QCT-VT model
            ECC=(ECT+EVIB)/EVOLT                !available collision energy in eV

            ! Han: add nonVHS cross section model for N2+O
            ! This model has dependency on vibrational level
            IF (nonVHS == 1) THEN
              IF ((LS == 1 .and. MS == 4) .or. (LS ==4 .and. MS == 1)) THEN
                CALL CALC_TOTXSEC(LS, MS, VR, VRR, ET0, EVIB, SXSECTION,BMAX,CVR)
              END IF
            END IF

!--Calculate the reaction cross-sections
            IVDC=0                 !to track the order of reacting species
            RXSECTION=0.d0
            IF (MNRE>0 .and. IMF .ne.  0 .and. QCTMODEL == 3) THEN
              ! precalculate Reaction cross section before collision occur
              ! This is for the case when QCT-SSD/SSE is used
              ! Reaction cross sections might become larger than others
              CALL CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,ECT,SXSECTION,IVDC,IDT)
            END IF
!


            ! QCT VT cross sections
            VTXSECTION=0.d0
            IF (QCTMODEL>=2.AND.GASCODE==8.AND.J>0) THEN
!
              IF(J==1) THEN  !N2+O collisions
                SIGMA_REF=14.5351d0/3.d0 !24.2788d0/3.d0
                QNU=0.7074d0 !0.3994d0
                A1=1.0570d0  !0.1253d0
                A2=-0.0503d0 !0.0221d0
                B1=8.2246d0  !5.8506d0
                B2=0.1639d0  !0.1398d0
                C1=0.4170d0  !0.2493d0
                C2=-0.3942d0 !-0.2075d0
                D=0.6974d0   !0.6988d0
                E=-7.3212d0  !-5.0131d0
                F=-0.0273d0  !-0.0868d0
                ED=9.8216d0  !E_dissociation in eV
              END IF
!
              IF(J==3) THEN  !O2+O collisions
                SIGMA_REF=6.0108d0/27.d0
                QNU=0.6390d0
                A1=0.0331d0
                A2=0.0016d0
                B1=2.0921d0
                B2=0.2459d0
                C1=-0.1357d0
                C2=0.4272d0
                D=0.5859d0
                E=-3.1009d0
                F=23.3814d0
                ED=5.21275d0  !E_dissociation in eV
              END IF
!
!
              A=A1+A2*ECC
              B=B1+B2*ECC
              C=C1+C2*ECC
!
              SUMF=0.d0
              EL=DMIN1(ECC,ED)                    !in eV
!
              IF(J==1) THEN !N2+O collisions
                AB=ECT+EVIB
                CALL VIB_LEVEL(AB,MAXLEV,1,J)              !max level based on the collision energy
                DO IV=0,MAXLEV
                  CALL VIB_ENERGY(EV_POST,IV,1,J)
                  AB=DABS(EVIB-EV_POST)/EVOLT              !energy difference in eV
                  BB=DABS(AB/C)
                  EF=(EV_POST/EVOLT)/ECC                   !energy ratio
                  IF(EF > 1.d0) EF=1.d0                    !to avoid round off errors when EV_POST=ECC
                  AAA=EL-(EVIB/EVOLT)
                  BBB=EL-(EV_POST/EVOLT)
                  S=A*AB+B*DEXP(-BB**D)+E*(AAA*BBB)**F     !surprisal function
                  VTXSECTION(IV)=((1.d0-EF)**QNU)*DEXP(S)
                  SUMF=SUMF+(1.d0-EF)**QNU                 !normalization constant
                END DO
              END IF
!
              IF(J==3) THEN !O2+O collisions
                MAXLEV=IVMODEL(J,2)
                DO IV=0,MAXLEV
                  CALL VIB_ENERGY(EV_POST,IV,1,J)
                  AB=DABS(EVIB-EV_POST)/EVOLT             !energy difference in eV
                  BB=DABS(AB/C)
                  EF=(EV_POST/EVOLT)/ECC                  !energy ratio
                  IF(EF > 1.d0) EF=1.d0                   !to avoid round off errors when EV_POST=ECC
                  AAA=(EV_POST/EVOLT)/ED-1.d0
                  BBB=(EVIB/EVOLT)/ED-1.d0
                  S=A*AB+B*DEXP(-BB**D)+E*DEXP(-F*AAA*BBB) !surprisal function
                  VTXSECTION(IV)=((1.d0-EF)**QNU)*DEXP(S)
                  SUMF=SUMF+(1.d0-EF)**QNU                 !normalization constant
                END DO
              END IF
!
              VTXSECTION(:)=VTXSECTION(:)/SUMF                     !normalized
              VTXSECTION(:)=SIGMA_REF*ET0**(QNU-1.d0)*VTXSECTION(:) !SIGMA_ME-QCT-VT in Angstrons^2
            END IF
!
!--Define the total cross-section
            TXSECTION=DMAX1(SXSECTION,SUM(VTXSECTION(:)))      !trick
            IF (MNRE>0 .and. IMF .ne.  0 .and. QCTMODEL == 3) THEN
              TXSECTION=DMAX1(TXSECTION,SUM(RXSECTION(:)))
            END IF
            !TXSECTION=DMAX1(TXSECTION,SUM(VTXSECTION(:)))     !total cross-section in Angstrons^2

            IF (TXSECTION*VRI/1.d20 > CCELL(4,N)) THEN
              !IF (((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) .and. nonVHS ==1) THEN
                !WRITE(*,*) TXSECTION*VRI/1.d20/CCELL(4,N)
              !ENDIF
              ! WRITE(*,*) TXSECTION*VRI/1.d20/CCELL(4,N)
              CCELL(4,N)=TXSECTION*VRI/1.d20 !update maximum value in (m2)*(m/s)
            END IF
!
!--NTC approach based on total cross-section
            PROB=(TXSECTION*VRI/1.d20)/CCELL(4,N)

            CALL ZGF(RANF,IDT)
            IF (PROB > RANF) THEN !collision occurs
!
!--Sample and initialize collision data
!-------------------------------------------------
              IF ((M==IPCP(L)).AND.(L==IPCP(M))) TOTDUP=TOTDUP+1.D00
              TOTCOL=TOTCOL+1.D00
              TCOL(LS,MS)=TCOL(LS,MS)+1.D00
              TCOL(MS,LS)=TCOL(MS,LS)+1.D00
              CSCR(LS,MS)=CSCR(LS,MS)+TXSECTION*VRI/1.d20 !accumulate species-dependent SigmaTotal*Cr (in m3/s)
              CSCR(MS,LS)=CSCR(MS,LS)+TXSECTION*VRI/1.d20
              COLLS(NN)=COLLS(NN)+1.D00
              WCOLLS(NN)=WCOLLS(NN)+WFC
              IF(NN == NSPDF) PCOLLS=PCOLLS+WFC
              SEP=DABS(PX(L)-PX(M))
              CLSEP(NN)=CLSEP(NN)+SEP

              !
              !--Sample internal states and energy
              !-------------------------------------------------
              ! sample t,r,v energy for particles that collide
              IF ((IREAC .eq. 2) .and. (IMF .ne. 0) .and. (MNRE > 0) .and. (IMFS == 1)) THEN
                IF (MNRE > 2) THEN
                  WRITE(*,*) 'MNRE should be less than 3 when IMFS is enabled'
                END IF
                NMFET(IETDX,IMFpair(LS,MS)) = NMFET(IETDX,IMFpair(LS,MS)) + 1.0D0
                ! KK=1 the one with lower smaller sp number
                DO KK = 1,2
                  K = KK
                  ! for same specie, X(:,2,:) is always zero
                  ! rotational energy
                  NMFER(IERDX(KK),K,IMFpair(LS,MS)) = NMFER(IERDX(KK),K,IMFpair(LS,MS)) + 1.0D0
                  ! vibrational energy
                  NMFEV(IEVDX(KK),K,IMFpair(LS,MS)) = NMFEV(IEVDX(KK),K,IMFpair(LS,MS)) + 1.0d0
                  NMFVT(IEVDX(KK),IETDX,K,IMFpair(LS,MS)) = NMFVT(IEVDX(KK),IETDX,K,IMFpair(LS,MS))+1.0D0
                END DO
              END IF
!
!--Check chemical reaction
!-------------------------------------------------
              IKA=0                  !becomes >0 if reaction occurs
              JX=2                   !standard number of post reaction species
              LMS=1                  !initialize
              IF (MNRE>0) THEN
                IF (QCTMODEL .ne. 3) THEN
                  RXSECTION = 0.0d0
                  ! For non QCT model, the function is called here
                  CALL CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,ECT,SXSECTION,IVDC,IDT)
                  ! for a chemical reaction K
                  ! IVDC(K) = 1:  LS = IREA(1,K), MS = IREA(2,K)
                  ! IVDC(K) = 2:  LS = IREA(2,K), MS = IREA(1,K)
                  ! For the case IREA(1,K) = IREA(2,K)
                  ! IVDC(K) = 1:  Molecule L is the one to be dissociated
                  ! IVDC(K) = 2:  Molecule M is the one to be dissociated
                END IF
                SUMF=SUM(RXSECTION(:)) !total reaction cross-section in Angstron^2

                PROB=SUMF/TXSECTION
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) THEN  !reaction occurs
                  IF (PROB > 1.d0) write(*,*) 'Prob_RXN above unity:',PROB,SUMF,TXSECTION
                  J=0
                  A=0.d0
                  CALL ZGF(RANF,IDT)
                  DO WHILE (A < RANF)
                    J=J+1
                    A=A+RXSECTION(J)/SUMF
                  END DO
                  IKA=IRCD(J,LS,MS) !reaction IKA occurs
                  CALL CHECK_REACTION(IKA,N,L,LM,M,LS,LMS,MS,VRR,VR,VRC,VRI,VCM,RML,RMM,IVDC,LOCAL_NPVIB,LOCAL_EVREM)
                  IF (IKA > 0) THEN
                    DO I=1,MMVM
                      NPVIB(1,IKA,1,I,:) = NPVIB(1,IKA,1,I,:) + LOCAL_NPVIB(1,I,:)
                      NPVIB(1,IKA,2,I,:) = NPVIB(1,IKA,2,I,:) + LOCAL_NPVIB(2,I,:)
                      NPVIB(1,IKA,3,I,:) = NPVIB(1,IKA,3,I,:) + LOCAL_NPVIB(3,I,:)
                    END DO
                    EVREM(IKA) = EVREM(IKA) + LOCAL_EVREM
                  END IF

                  ! VRR, VR become the velocites make 0.5*mr*VRR = Etot
                  ! For recombination, VRC is changed from A-B relative velocity to AB-C, the same applies to VCM
                END IF
              END IF
!
!--Update variables for VRT energy redistribution
!-------------------------------------------------
              IF (IKA /= 0) THEN
                NPM=NREA(1,IKA)+NREA(2,IKA)
                IF (NPM == 3) JX=3          !dissociation
                IF (IREAC == 0) THEN
                  REAC(IKA)=REAC(IKA)+1.D00 !count reactions in entire domain
                ELSE
                  IF (NN == NSPDF) REAC(IKA)=REAC(IKA)+1.D00 !count reactions only in NSPDF cell
                END IF
                IF (IPRS == 1) THEN
                  DO IT=1,ITMAX
                    DT(IT)=DABS(VAR(8,NN)-FPTEMP(IT))
                  END DO
                  IT=MINLOC(DT,1)
                END IF
                IF (IREAC == 2) THEN
                  IKA=0
                  JX=2
                  IVDC = 0
                  IVDC0 = 0
                ELSE
                  IVDC0 = IVDC(IKA)
                END IF
              ELSE
                IVDC0 = 0
              END IF

!             !ECT is updated here
              ! :cry:
              ! TLDR, ECT is precollision translational energy for nonreaction case
              !       ECT is total energy available to distribute for reactive case
              !
              ECT=(0.5D00*SPM(1,LS,MS)*VRR)
              IF (ECT < 0.d0) ECT = 1.d-26

              !! ======== Israel include this ELACK due to recombination reaction, but this ELACK =====
              !!         can cause data race when openmp is enabled
              ! ECT=ECT-ELACK
              ! IF (ECT > 0.d0) THEN
              !   ELACK=0.d0
              ! ELSE
              !   ELACK=-ECT+1.d-26
              !   ECT=1.d-26
              ! END IF
              ! ========================================================================================
!
              CALL ZGF(RANF,IDT)
              IF (RANF > 0.5d0) THEN  !shuffle the sequence of molecules in energy redistribution
                IAX=1;  KK=JX; IVM=1
              ELSE
                IAX=JX; KK=1;  IVM=-1
              END IF
!
!--Check VT energy redistribution
!-------------------------------------------------
              J=0
              IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) J=1 !N2-O collision
              IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) J=3 !O2-O collision
!
              IF (IREAC<=1.AND.QCTMODEL>=2.AND.IKA==0.AND.GASCODE==8.AND.J>0) THEN
                ! The following code is run under the condition that
                !  1. Reaction can occur (IREA <=1)
                !  2. QCT-VT Modell is on
                !  3. Reaction doesn't occur (IKA = 0)
                !  4. GASCODE == 8
                !  5. the species is either N2-O or O2-O
              !--Vibrational relaxation based on ME-QCT model
!
                SUMF=SUM(VTXSECTION(:))      !total VT cross-section (SIGMA_ME-QCT-VT,T)
                PROB=SUMF/TXSECTION          !SIGMA_ME-QCT-VT,T to SIGMA_TOTAL ratio
                IF (PROB > 1.d0) write(*,*) 'Prob_VT above unity:',PROB
!
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) THEN        !vibrational relaxation occurs
                  II=0
                  MAXLEV=IVMODEL(J,2)
                  DO WHILE (II == 0)
                    CALL ZGF(RANF,IDT)
                    IV=RANF*(MAXLEV+0.99999999D00)
                    PROB=VTXSECTION(IV)/SUMF !SIGMA_ME-QCT-VT to SIGMA_ME-QCT-VT,T ratio
                    CALL ZGF(RANF,IDT)
                    IF (PROB > RANF) THEN    !a vibrational level is accepted
                      CALL VIB_ENERGY(EV_POST,IV,1,J)
                      II=1
                    END IF
                  END DO
!
                  IF (LS == J) IPVIB(1,L)=IV
                  IF (MS == J) IPVIB(1,M)=IV
!
                  ! for diatom-atom ME-QCT-VT, ECC is the collisional energy + vibrational energy of diatom molecule
                  ! ECC is in eV
                  ! Note that, ME-QCT-VT don't work when reaction occurs
                  IF (IRELAX /= 0) ECT=(ECC*EVOLT)-EV_POST  !remaining energy available for redistribution; comment this line for isothermal relaxations
                  ! ECT don't have precollision rotational energy
                END IF
!
              ELSE
!-------------------------------------------------
                !--Larsen-Borgnakke serial vibrational energy redistribution
                REST_DOF = 0.0d0
                IF ((IKA /= 0).AND.(IPRS == 2)) THEN
                      A=VAR(8,NN)  !sampling cell translational temperature
                      SVDOF=0.d0
                      IF (ISPV(LS) > 0) THEN
                        DO KV=1,ISPV(LS)
                          SVDOF(1,KV)=2.d0*(SPVM(1,KV,LS)/A)/(DEXP(SPVM(1,KV,LS)/A)-1.d0)
                        END DO
                      END IF
                      IF (ISPV(MS) > 0) THEN
                        DO KV=1,ISPV(MS)
                          SVDOF(2,KV)=2.d0*(SPVM(1,KV,MS)/A)/(DEXP(SPVM(1,KV,MS)/A)-1.d0)
                        END DO
                      END IF
                      IF (JX == 3) THEN
                        IF (ISPV(LMS) > 0) THEN
                          DO KV=1,ISPV(LMS)
                            SVDOF(3,KV)=2.d0*(SPVM(1,KV,LMS)/A)/(DEXP(SPVM(1,KV,LMS)/A)-1.d0)
                          END DO
                        END IF
                      END IF
                      REST_DOF=5.d0-2.d0*SPM(3,LS,MS)   !translational dofs
                      REST_DOF=REST_DOF+ISPR(1,LS)+ISPR(1,MS) !rotational dofs
                      REST_DOF=REST_DOF+SUM(SVDOF)            !vibrational dofs
                END IF
!
                DO NSP=IAX,KK,IVM
                      IF (IVDC0 < 2) THEN
                        IF (NSP == 1) THEN
                          K=L ; KS=LS ; JS=MS
                        ELSE IF (NSP == 2) THEN
                          K=M ; KS=MS ; JS=LS
                        ELSE
                          K=LM ; KS=LMS ; JS=MS  !new born molecule from dissociation reaction
                        END IF
                      ELSE
                        IF (NSP == 1) THEN
                          K=M ; KS=MS ; JS=LS
                        ELSE IF (NSP == 2) THEN
                          K=L ; KS=LS ; JS=MS
                        ELSE
                          K=LM ; KS=LMS ; JS=LS  !new born molecule from dissociation reaction
                        END IF
                      END IF
!
                      IF (MMVM > 0) THEN
                        IF (ISPV(KS) > 0) THEN
                          IF (ISPV(KS) == 3) THEN  !shuffle the sequence of vibrational modes over serial LB
                            CALL ZGF(RANF,IDT)
                            JI=RANF*(6.d0+0.9999999d0)
                            SELECT CASE (JI)
                              CASE DEFAULT; ISHUF(1)=1; ISHUF(2)=2; ISHUF(3)=3;
                              CASE (1); ISHUF(1)=1; ISHUF(2)=2; ISHUF(3)=3;
                              CASE (2); ISHUF(1)=1; ISHUF(2)=3; ISHUF(3)=2;
                              CASE (3); ISHUF(1)=2; ISHUF(2)=3; ISHUF(3)=1;
                              CASE (4); ISHUF(1)=2; ISHUF(2)=1; ISHUF(3)=3;
                              CASE (5); ISHUF(1)=3; ISHUF(2)=1; ISHUF(3)=2;
                              CASE (6); ISHUF(1)=3; ISHUF(2)=2; ISHUF(3)=1;
                            END SELECT
                          END IF
                          DO I=1,ISPV(KS)
                            KV=I
                            IF (ISPV(KS) == 3) KV=ISHUF(I)
                            IF (IPVIB(KV,K) >= 0) THEN
!--do not redistribute to a dissociating molecule marked for removal or for a recombined molecule
                              IF (IKA .eq. 0) THEN
                                ! IF IKA == 0, ECT is the precollision translational energy
                                ! This part of the code is called when no reaction occurs
                                ! and LB model is used to redistribute energy
                                CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
                                ECC=ECT+EVIB
                              ELSE
                                ECC=ECT
                                ! This part of the code is called when reaction
                                ! occurs and LB model distributes internal energy of each molecule
                                ! from the total energy available. ECC was updated around line 11651
                                ! In Israel's old version, he didn't differenciate these two cases
                                ! because IPVIB~=0 => EVIB~=0. But it's not true
                              END IF
                              CALL VIB_LEVEL(ECC,MAXLEV,KV,KS)
                              IF (SPVM(3,KV,KS) > 0.) THEN
                                IF (IZV == 1) THEN
                                  COLT=(DFLOAT(MAXLEV)*SPVM(1,KV,KS))/(3.5D00-SPM(3,KS,JS))  !--quantized collision temperature
                                ELSE
                                  COLT=VAR(8,NN)   !--isebasti: Ttrans, according to DSMC.f90 (2013)
                                END IF
                                B=SPVM(4,KV,KS)/SPVM(3,KV,KS)    !Tdiss/Tref
                                A=SPVM(4,KV,KS)/COLT             !Tdiss/T
       ZV=(A**SPM(3,KS,JS))*(SPVM(2,KV,KS)*(B**(-SPM(3,KS,JS))))**(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                              ELSE
                                IF (SPVM(3,1,KS) == -1.) THEN
                                  ZV=SPVM(2,KV,KS)
                                ELSE IF (SPVM(3,1,KS) == -2.) THEN
                                  IF (KS == JS) THEN
                                    ZV=SPVM(2,1,KS)
                                  ELSE
                                    ZV=SPVM(2,2,KS)
                                  END IF
                                END IF
                              END IF
                              IF (IREAC<=1) THEN
!VAR(8,NN)=0.154*SPVM(1,KV,KS)
!do II=0,20
!VAR(10,NN)=3.5d0*SPVM(1,KV,KS)*DFLOAT(ii)/20.d0
                                IF (VAR(8,NN) >= VAR(11,NN)) THEN
                                  A=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  A=0.5d0*A*A*DEXP(SPVM(1,KV,KS)/VAR(8,NN))
                                  B=5.d0-2.d0*SPM(3,LS,MS)
                                  C=B/(B+A) !C is the Zcorrection
                                ELSE
                                  A=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  B=2.d0*(SPVM(1,KV,KS)/VAR(10,NN))/(DEXP(SPVM(1,KV,KS)/VAR(10,NN))-1.d0)
                                  IF(VAR(10,NN) < 1.d-3) B=0.d0
                                  C=5.d0-2.d0*SPM(3,LS,MS)
                                  D=B*VAR(10,NN)+C*VAR(8,NN)
                                  CALL ROOTF_SECANT(2,KS,LS,MS,VAR(11,NN),D,E)
                                  C=(C*(VAR(8,NN)-E))/(A*VAR(8,NN)-B*VAR(10,NN)) !C is the Zcorrection
                                END IF
!write(*,*) VAR(10,NN)/SPVM(1,KV,KS),VAR(8,NN)/SPVM(1,KV,KS),E/SPVM(1,KV,KS),C
!end do
!stop

                                ! =========== Israel's correction of VT ====================
                                ! O2-O2 omega = 0.71, dref=3.985
                                ! O2-O  omega = 0.75, dref=3.442
                                ! IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) THEN
                                !   !calibrated with O2-O MeQct model
                                !   IF (VAR(8,NN)<= 1.d3) ZV=C*(60.01d0+0.04268d0*VAR(8,NN)-2.29d-5*VAR(8,NN)**2.d0)
                                !   IF (VAR(8,NN) > 1.d3) ZV=C*DEXP(-0.1372d0*DLOG(VAR(8,NN))+5.328d0)
                                ! END IF
                                ! IF (LS==3.AND.MS==3) THEN
                                !   !calibrated with O2-O2 ibraguimova2013 data
                                !   A=VAR(8,NN)**(-1.d0/3.d0)                          !T^(-1/3)
                                !   B=7.206d0-289.4d0*A+4919.d0*A**2.d0-2.23d4*A**3.d0 !log10(Zv)
                                !   ZV=C*10.d0**B
                                ! END IF

                                ! ============ Han's correction ===========================
                                ! MF-DSMC-Correct by Han Luo
                                IF (LS==3 .and. MS==3 .and. IMF .ne. 0) THEN
                                  IF ( IMF == 1) THEN
                                    ! MF-DSMC-SHO model, the following tries to match Ibraguimova 2013's VT relaxation time
                                    ! should be used with corrected total cross sections, check around L 11187
                                    ! and omega = 0.99531, dref = 17.474 A
                                    A=VAR(8,NN)**(-1.d0/3.d0)                          !T^(-1/3)
                                    IF (VAR(8,NN) .ge. 6.0d3) THEN
                                        B=5.2258d0-232.766d0*A+4.6605d3*A**2 - 2.22176d4*A**3
                                        ZV = 10.0d0**B
                                    ELSE
                                        ZV=4.5134d5*VAR(8,NN)**(-1.795d0)/(1.0D0-DEXP(-2238.D0/VAR(8,NN)))*DEXP(144.123d0*A)
                                    END IF
                                  ELSE IF (IMF == 3) THEN
                                    ! MF-DSMC-AHO Morse
                                    ! dref = 4.1515, omega = 0.7318
                                    A=VAR(8,NN)**(-1.d0/3.d0)                          !T^(-1/3)
                                    B=1.767075d5*A**4 - 6.347191d4*A**3 + 8.422653d3*A**2 - 4.179710d2*A + 8.891025d0
                                    ZV = 10.0d0**B
                                  END IF
                                  ZV=C*ZV
                                END IF
                              END IF
                              CALL ZGF(RANF,IDT)
                              IF ((1.D00/ZV > RANF).OR.(IKA /= 0)) THEN
!-- a vibrational redistribution is always made after a chemical reaction
                                II=0
                                nstep=0
                                IF ((IKA /= 0).AND.(IPRS == 2)) THEN  !modified LB for reacting collisions
                                  B=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  REST_DOF=REST_DOF-B
                                  A=REST_DOF    !remaining degrees of freedom
                                ELSE
                                  A=5.d0-2.d0*SPM(3,KS,JS)              !nonreacting collisions
                                END IF

                                DO WHILE ((II == 0).and.(nstep < 100000))
                                  nstep=nstep+1
                                  if (nstep > 99000) then
                                    write (*,*) nstep,ecc,maxlev
                                    stop
                                  end if
                                  CALL ZGF(RANF,IDT)
                                  IV=RANF*(MAXLEV+0.99999999D00)
                                  IPVIB(KV,K)=IV
                                  CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
                                  !EVIB=DFLOAT(IV)*BOLTZ*SPVM(1,KV,KS)
                                  IF (EVIB < ECC) THEN
                                    CALL INELASTIC_VT(KS, JS, EVIB, ECC, A, PROB)
                                    ! PROB=(1.D00-EVIB/ECC)**(A*0.5d0-1.d0)
 !--PROB is the probability ratio of eqn (5.61) and (5.45)
                                    CALL ZGF(RANF,IDT)
                                    IF (PROB > RANF) II=1
                                  END IF
                                END DO
                                IF (IRELAX /= 0) ECT=ECC-EVIB !remaining energy available for redistribution; comment this line for isothermal relaxations
                              END IF
                            END IF
                          END DO
                        END IF
                      END IF
                END DO
              END IF  !--end of VT exchange
!
!--sample post-collisional vibrational levels
!-------------------------------------------------
              IF (IKA /= 0) THEN
                    DO KV=1,MMVM
                      IF (IVDC0 == 1) THEN
                        SELECT CASE(JX)
                        CASE(2) !recombination and exchange
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,M)
                          K=0
                        CASE(3) !dissociation
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,LM) !new born
                          K=IPVIB(KV,M)  !third body
                        END SELECT
                      ELSE
                        SELECT CASE(JX)
                        CASE(2) !recombination and exchange
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,M)
                          K=0
                        CASE(3) !dissociation
                          I=IPVIB(KV,M)
                          J=IPVIB(KV,LM) !new born
                          K=IPVIB(KV,L)  !third body
                        END SELECT
                      END IF
                      IF (I > 100) I=100
                      IF (J > 100) J=100
                      IF (K > 100) K=100
                      NPVIB(2,IKA,1,KV,I)=NPVIB(2,IKA,1,KV,I)+1  !bin counter
                      NPVIB(2,IKA,2,KV,J)=NPVIB(2,IKA,2,KV,J)+1
                      NPVIB(2,IKA,3,KV,K)=NPVIB(2,IKA,3,KV,K)+1
                    END DO
              END IF
              ! up to now
              ! nonreacting case: ECT collisional energy
              ! reacting case: ECT collisional energy + rotational energy to distribute
!
!--Check RT energy redistribution (LB model)
!-------------------------------------------------
              DO NSP=IAX,KK,IVM
                    IF  (IVDC0 < 2) THEN
                      IF (NSP == 1) THEN
                        K=L ; KS=LS ; JS=MS
                      ELSE IF (NSP == 2) THEN
                        K=M ; KS=MS ; JS=LS
                      ELSE
                        K=LM ; KS=LMS ; JS=MS  !new born molecule from dissociation reaction
                      END IF
                    ELSE
                      IF (NSP == 1) THEN
                        K=M ; KS=MS ; JS=LS
                      ELSE IF (NSP == 2) THEN
                        K=L ; KS=LS ; JS=MS
                      ELSE
                        K=LM ; KS=LMS ; JS=LS  !new born molecule from dissociation reaction
                      END IF
                    END IF
                    IF (ISPR(1,KS) > 0) THEN
                      IF (ISPR(2,KS) == 0) THEN
                        B=SPM(7,KS,JS) !Zr
                      ELSE
                        COLT=100.d0/VAR(8,NN) !T*/Ttrans; using SMILE data
                        B=20.d0/(1.d0+(0.5d0*PI**(3.d0/2.d0))*DSQRT(COLT)+PI*(1.d0+PI/4.d0)*COLT) !Zr
                      END IF
                      A=5.d0-2.d0*SPM(3,KS,LS)
                      C=A/(A+ISPR(1,KS)+ISPR(1,JS)) !Zr correction
                      A=B*C                         !Corrected Zr
                      CALL ZGF(RANF,IDT)
                      IF ((1.d0/A > RANF).OR.(IKA /= 0)) THEN
                        ! we don't distinguish for the value of IKA because PROT is REALLY zero
                      ! for reacting case
                        ECC=ECT+PROT(K)
                        CALL INELASTIC_RT(KS, JS, ECC, ERM, IDT)
                        PROT(K)=ERM*ECC
                        ECT=ECC-PROT(K)
                      END IF
                    END IF
              END DO
              ! At this point, energy left in ECT is just collisional energy
!
!--adjust VR after energy redistribution
!-------------------------------------------------
              IF (JX == 2) THEN
                VR=DSQRT(2.D00*ECT/SPM(1,LS,MS))
              ELSE
                ! Reaction off
                VR=DSQRT(2.D00*ECT/SPM(1,IREA(1,IKA),IREA(2,IKA))) !same as in SPARTA
              END IF
!
!--calculate new velocities
!-------------------------------------------------
              ! At this point for exchange or dissocation reaction
              ! VRC and VRI are still the original relative velocity of the two molecule
              ! For recombination, it's the relative velocity of the
              ! new molecule and the third body
              ! The values are changed in CHECK_REACTION
              VRC = VRC/VRI*VR             ! rescale to match current magnitude
              IF ( .not. IVHS ) THEN
                CALL ZGF(RANF,IDT)
                IF (RANF .le. (CVR2*1d20/VRI/TXSECTION)) THEN
                   IVHS = .true.
                END IF
                ! either by probability or if CVR2>CVR
              END IF

              ! The current implementation may not handle recombination well for EXP scattering
              IF (IVHS) THEN
                IF (JX == 2) THEN
                  CALL SCATTER_MOL(LS, MS, BMAX, ET0, VRI, VR, VRC, VRCP, IDT)
                ELSE
                  CALL SCATTER_MOL(IREA(1,IKA), IREA(2,IKA), BMAX,ET0, VRI, VR, VRC, VRCP, IDT)
                END IF
              ELSE
                VRCP = VRC
              END IF

              PV(1:3,L)=VCM(1:3)+RMM*VRCP(1:3)
              PV(1:3,M)=VCM(1:3)-RML*VRCP(1:3)
              IF (JX == 3) THEN
                IF (IVDC0 == 1) THEN
                  PV(1:3,LM)=PV(1:3,L)
                  IPCP(LM)=L
                ELSE
                  PV(1:3,LM)=PV(1:3,M)
                  IPCP(LM)=M
                END IF
              END IF
              IPCP(L)=M
              IPCP(M)=L
!
              !--to sustain a bimodal translational-rotatioanal distribution
              !IF (IGS==4) THEN
              !  DO NSP=1,2
              !    IF (NSP == 1) THEN
              !      K=L ; KS=LS
              !    ELSE
              !      K=M ; KS=MS
              !    END IF
              !    CALL ZGF(RANF,IDT)
              !    IF (RANF > 0.5d0) THEN
              !      A=FTMP(1)   !based on translational reference temperature
              !    ELSE
              !      A=FRTMP(1)  !based on rotational reference temperature
              !    END IF
              !    B=SQRT(2.D00*BOLTZ*A/SP(5,KS))
              !    DO J=1,3
              !      CALL RVELC(PV(J,K),C,B,IDT) !for bimodal distribution
              !    END DO
              !    PV(1,K)=PV(1,K)+VFX(1)
              !    PV(2,K)=PV(2,K)+VFY(1)
              !    IF (ISPR(1,KS) > 0) CALL SROT(KS,A,PROT(K),IDT) !for bimodal distribution
              !  END DO
              !END IF
!
            END IF   !--end of collision
          END IF   !--separate simple gas / mixture cases
        END DO   !--do for all candidate collisions
        ! MVELC = 0.0d0
        ! MSVELC = 0.0d0
        ! DO K = 1, ICCELL(2,N)
          ! M = ICREF(ICCELL(1,N) + K)
          ! LS = ABS(IPSP(M))
          ! DO J = 1,3
            ! MVELC(J) = MVELC(J) + PV(J,M)
            ! MSVELC(J) = MSVELC(J) + PV(J,M)**2
          ! END DO
        ! END DO
        ! DO J=1,3
          ! MVELC(J) = MVELC(J) / ICCELL(2,N)
          ! MSVELC(J) = MSVELC(J) / ICCELL(2,N)
          ! VARVELC(J) = MSVELC(J) - MVELC(J)**2
        ! END DO
        ! WRITE(19,'(I3,1X,20G13.5)') N, (MVELC(J),J=1,3), (MSVELC(J),J=1,3), (VARVELC(J),J=1,3)
      END IF
    END IF
  END IF
END DO
! WRITE(9,*)
! WRITE(9,*)
! CLOSE(9)
!$omp end do
!$omp end parallel
!
!--remove any recombined atoms
N=0
DO WHILE (N < NM)
  N=N+1
  IF (IPCELL(N) < 0) THEN
    CALL REMOVE_MOL(N)
    N=N-1
  END IF
END DO
!
!--dissociate marked molecules
!IF (MNRE > 0) CALL DISSOCIATION
!
RETURN
END SUBROUTINE COLLISIONS
