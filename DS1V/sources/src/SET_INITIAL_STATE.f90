!
!***************************************************************************
!
SUBROUTINE SET_INITIAL_STATE
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE EXPCOL, only:EXPCOL_INIT
!
IMPLICIT NONE
!
INTEGER :: I,J,L,K,KK,KN,NMI,II,III,INC,NC,NDES,IDES(10000),IDT=0 !--isebasti: included NC,IDT
INTEGER(KIND=8) :: N,M
REAL(KIND=8) :: AA,A,B,BB,SN,XMIN,XMAX,WFMIN,DENG,TDF,DNF,VMP2,RANF,EVIB,TVAR(10000),BVMP,BTEMP !--isebasti: included TDF,DNF,VMP2,RANF,EVIB
REAL(KIND=8), DIMENSION(3) :: DMOM
REAL(KIND=8), DIMENSION(2) :: ROTE
REAL(8),EXTERNAL :: ERF,GAM
!
!--NSET the alternative set numbers in the setting of exact initial state
!--DMOM(N) N=1,2,3 for x,y and z momentum sums of initial molecules
!--DENG the energy sum of the initial molecules
!--ROTE alternative sets of rotational energy
!--NMI the initial number of molecules
!--INC counting increment
!
DMOM=0.D00
DENG=0.D00
!--set the number of molecules, divisions etc. based on stream 1
!
NMI=10000*AMEG+0.999999999D00
NDIV=NMI/MOLSC !--MOLSC molecules per division
WRITE (9,*) 'The number of divisions is',NDIV
!
MDIV=NDIV
ILEVEL=0
!
ALLOCATE (JDIV(0:ILEVEL,MDIV),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR JDIV ARRAY',ERROR
ENDIF
!
DDIV=(XB(2)-XB(1))/DFLOAT(NDIV)
NCELLS=NDIV
WRITE (9,*) 'The number of sampling cells is',NCELLS
NCCELLS=NCIS*NDIV
WRITE (9,*) 'The number of collision cells is',NCCELLS
!
!IF (IFX == 0) XS=0.D00  !--isebasti: do not reset XS read from READ_DATA
!
IF (ISECS == 0) THEN
  IF (IFX == 0) FNUM=(XB(2)-XB(1))*FND(1)/DFLOAT(NMI)
  IF (IFX == 1) FNUM=PI*(XB(2)**2-XB(1)**2)*FND(1)/DFLOAT(NMI)
  IF (IFX == 2) FNUM=1.3333333333333333333333D00*PI*(XB(2)**3-XB(1)**3)*FND(1)/DFLOAT(NMI)
ELSE
  IF (IFX == 0) FNUM=((XS-XB(1))*FND(1)+(XB(2)-XS)*FND(2))/DFLOAT(NMI)
  IF (IFX == 1) FNUM=PI*((XS**2-XB(1)**2)*FND(1)+(XB(2)**2-XS**2)*FND(2))/DFLOAT(NMI)
  IF (IFX == 2) FNUM=1.3333333333333333333333D00*PI*((XS**3-XB(1)**3)*FND(1)-(XB(2)**3-XS**3)*FND(2))/DFLOAT(NMI)
END IF
!
FTIME=0.D00
TLIM=1.D20
!
TOTMOV=0.D00
TOTCOL=0.D00
PCOLLS=0.D00
TOTDUP=0.D00
TCOL=0.D00
CSCR=0.D00
TDISS=0.D00
TRECOMB=0.D00
TFOREX=0.D00
TREVEX=0.D00
TREACG=0
TREACL=0
TNEX=0.D00
!
DO N=1,NDIV
  JDIV(0,N)=-N
END DO
!
ALLOCATE (CELL(4,NCELLS),ICELL(NCELLS),CCELL(5,NCCELLS),ICCELL(3,NCCELLS),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR CELL ARRAYS',ERROR
ENDIF
!
ALLOCATE (COLLS(NCELLS),WCOLLS(NCELLS),CLSEP(NCELLS),REAC(MNRE),SREAC(MNSR),VAR(21,NCELLS), &
          VARSP(0:11,NCELLS,MSP),VARS(0:32+MSP,2),CS(0:8+MMVM,NCELLS,MSP),CSS(0:8,2,MSP,2), &  !-isebasti: correcting CS allocation
          CSSS(6,2),CST(0:4,NCELLS),BINS(0:NBINS,5,MSP),BIN(0:NBINS,5),EVREM(MNRE),&
          PDFS(0:NBINS,5,MSP),PDF(0:NBINS,5), &
          NDROT(MSP,100),NDVIB(NSCELLS,0:MMVM,MSP,0:100),STAT=ERROR) !--isebasti: CST,PDFS,PDF,BINS,BIN included
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR SAMPLING ARRAYS',ERROR
ENDIF
!
!
CALL INITIALISE_SAMPLES
!
!--Set the initial cells

DO N=1,NCELLS
  CELL(2,N)=XB(1)+DFLOAT(N-1)*DDIV
  CELL(3,N)=CELL(2,N)+DDIV
  CELL(1,N)=CELL(2,N)+0.5D00*DDIV
  IF (IFX == 0) CELL(4,N)=CELL(3,N)-CELL(2,N)    !--calculation assumes unit cross-section
  IF (IFX == 1) CELL(4,N)=PI*(CELL(3,N)**2-CELL(2,N)**2)  !--assumes unit length of full cylinder
  IF (IFX == 2) CELL(4,N)=1.33333333333333333333D00*PI*(CELL(3,N)**3-CELL(2,N)**3)    !--flow is in the full sphere
  ICELL(N)=NCIS*(N-1)
  DO M=1,NCIS
    L=ICELL(N)+M
    XMIN=CELL(2,N)+DFLOAT(M-1)*DDIV/DFLOAT(NCIS)
    XMAX=XMIN+DDIV/DFLOAT(NCIS)
    IF (IFX == 0) CCELL(1,L)=XMAX-XMIN
    IF (IFX == 1) CCELL(1,L)=PI*(XMAX**2-XMIN**2)  !--assumes unit length of full cylinder
    IF (IFX == 2) CCELL(1,L)=1.33333333333333333333D00*PI*(XMAX**3-XMIN**3)    !--flow is in the full sphere
    CCELL(2,L)=0.D00
    ICCELL(3,L)=N
  END DO
  VAR(8,N)=FTMP(1)
  VAR(9,N)=FRTMP(1)
  VAR(10,N)=FVTMP(1)
  VAR(11,N)=FTMP(1)
END DO
!
IF (IWF == 0) AWF=1.D00
IF (IWF == 1) THEN
!--FNUM must be reduced to allow for the weighting factors
  A=0.D00
  B=0.D00
  DO N=1,NCELLS
    A=A+CELL(4,N)
    B=B+CELL(4,N)/(1.+WFM*CELL(1,N)**IFX)
  END DO
  AWF=A/B
  FNUM=FNUM*B/A
END IF
!
WRITE (9,*) 'FNUM is',FNUM
!
!--set the information on the molecular species
!
A=0.D00
B=0.D00
DO L=1,MSP
  A=A+SP(5,L)*FSP(L,1)
  B=B+(3.+ISPR(1,L))*FSP(L,1)
  VMP(L,1)=SQRT(2.D00*BOLTZ*FTMP(1)/SP(5,L))
  IF ((ITYPE(2)== 0).OR.(ISECS == 1)) VMP(L,2)=SQRT(2.D00*BOLTZ*FTMP(2)/SP(5,L))
  VNMAX(L)=3.*VMP(L,1)
  IF (L == 1) THEN
    VMPM=VMP(L,1)
  ELSE
    IF (VMP(L,1) > VMPM) VMPM=VMP(L,1)
  END IF
END DO
WRITE (9,*) 'VMPM =',VMPM
FDEN=A*FND(1)
FPR=FND(1)*BOLTZ*FTMP(1)
FMA=VFX(1)/DSQRT((B/(B+2.D00))*BOLTZ*FTMP(1)/A)
DO L=1,MSP
  DO M=1,MSP
    SPM(4,L,M)=0.5D00*(SP(1,L)+SP(1,M))  !d_ref
    SPM(3,L,M)=0.5D00*(SP(3,L)+SP(3,M))  !omega_ref
    SPM(5,L,M)=0.5D00*(SP(2,L)+SP(2,M))  !t_ref
    SPM(1,L,M)=SP(5,L)*(SP(5,M)/(SP(5,L)+SP(5,M)))
    SPM(2,L,M)=0.25D00*PI*(SP(1,L)+SP(1,M))**2
    AA=2.5D00-SPM(3,L,M)
    A=GAM(AA)
    SPM(6,L,M)=1.D00/A
    SPM(8,L,M)=0.5D00*(SP(4,L)+SP(4,M))
    IF ((ISPR(1,L) > 0).AND.(ISPR(1,M) > 0)) THEN
      SPM(7,L,M)=(SPR(1,L)+SPR(1,M))*0.5D00
    END IF
    IF ((ISPR(1,L) > 0).AND.(ISPR(1,M) == 0)) THEN
      SPM(7,L,M)=SPR(1,L)
    END IF
    IF ((ISPR(1,M) > 0).AND.(ISPR(1,L) == 0)) THEN
      SPM(7,L,M)=SPR(1,M)
    END IF
  END DO
END DO
!
IF (GASCODE==8) THEN !special code with wysong2014 data for O2-O collisions
    ! special settings for GASCODE==8
    ! O-O2
    SPM(3,3,4)=0.75d0      !omega_ref
    SPM(4,3,4)=3.442d-10  !d_ref
    SPM(5,3,4)=273.d0     !t_ref
    SPM(2,3,4)=PI*SPM(4,3,4)**2
    SPM(6,3,4)=1.D00/GAM(2.5D0 - SPM(3,3,4))
    DO L = 2,6
      SPM(L,4,3) = SPM(L,3,4)
    ENDDO

    ! N2-O2
    L = 3; M = 1;
    SPM(3,L,M)=0.7343d0      !omega_ref
    SPM(4,L,M)=4.0512d-10  !d_ref
    SPM(5,L,M)=273.d0     !t_ref
    SPM(2,L,M)=PI*SPM(4,L,M)**2
    SPM(6,L,M)=1.D00/GAM(2.5D0 - SPM(3,L,M))
    DO N = 2,6
      SPM(N,M,L) = SPM(N,L,M)
    ENDDO

    ! NO-N2
    L = 1; M = 5;
    SPM(3,L,M)=0.7292554d0      !omega_ref
    SPM(4,L,M)=4.34104d-10  !d_ref
    SPM(5,L,M)=273.d0     !t_ref
    SPM(2,L,M)=PI*SPM(4,L,M)**2
    SPM(6,L,M)=1.D00/GAM(2.5D0 - SPM(3,L,M))
    DO N = 2,6
      SPM(N,M,L) = SPM(N,L,M)
    ENDDO

    IF (nonVHS == 2) THEN
      CALL EXPCOL_INIT()
    ELSE IF (nonVHS == 1) THEN
      CCELL(4,:) = CCELL(4,:)*1.1D0
      INONVHS(1,4) = 1
      INONVHS(4,1) = 1
    ELSE IF (nonVHS == 3) THEN
      INONVHS = 0
    END IF
END IF

!IF (GASCODE == 8 )THEN
  !SPM(3,3,4) = 0.75d0
  !SPM(4,3,4) = 9.61877d-10
  !SPM(5,3,4) = 273.d0
  !SPM(2,3,4) = PI*SPM(4,3,4)**2
  !SPM(3,4,3) = SPM(3,3,4)
  !SPM(4,4,3) = SPM(4,3,4)
  !SPM(5,4,3) = SPM(5,3,4)
  !SPM(2,4,3) = SPM(2,4,3)

  !SPM(3,3,3) = 1.678d0
  !SPM(4,3,3) = 24.926d-10
  !SPM(5,3,3) = 273.d0
  !SPM(2,3,3) = PI*SPM(4,3,3)**2
!END IF



!
IF (MNRE > 0) CALL INIT_REAC
!
IF (MMVM > 0) THEN  !--isebasti: included
  IDL=0
  DO L=1,MSP
    IF (ISPV(L) > 0) THEN
      DO K=1,ISPV(L)
        IDL(K,L)=SPVM(4,K,L)/SPVM(1,K,L)  !assuming SHO
      END DO
    END IF
  END DO
END IF
!
IF (MNRE > 0) THEN
  DO L=1,MSP
    IF (ISPV(L) == 0) DPER(L)=0   !number of different vibrational energy levels
    IF (ISPV(L) == 1) DPER(L)=IDL(K,L)+1.d0
  END DO
!
  CPER(:,:)=1.d0  !CPER
  DO KK=1,MNRE
    IF (REA(5,KK) < 0.d0) THEN !consider only endothermic (e.g., dissocation) reactions
      DO KN=1,2
        L=IREA(KN,KK)
!
        IF (ISPV(L) == 1) THEN
          A=0.d0
          B=0.d0
          DO I=0,IDL(1,L)
            B=B+1.d0
            EVIB=DFLOAT(I)*BOLTZ*SPVM(1,1,L)
            IF (EVIB < DABS(REA(5,KK))) A=A+1.d0
          END DO
          CPER(KK,KN)= DPER(L)/A
        END IF
!
        IF (ISPV(L) == 3) THEN
          II=0
          DO I=0,IDL(1,L)
            DO J=0,IDL(2,L)
             DO K=0,IDL(3,L)
                II=II+1.d0    !number of possible states
                !stores the energy level for vibrational state (I,J,K)
                TVAR(II)=DFLOAT(I)*BOLTZ*SPVM(1,1,L)+DFLOAT(J)*BOLTZ*SPVM(1,2,L)+DFLOAT(K)*BOLTZ*SPVM(1,3,L)
              END DO
            END DO
          END DO
!
          NDES=0
          DO I=1,II
            DO J=1,II
              A=0.999999d0*TVAR(J)
              B=1.000001d0*TVAR(J)
              IF ((TVAR(I) > A ).AND.(TVAR(I) < B )) THEN
                NDES=NDES+1    !number of different energy states
                IDES(NDES)=I   !a different energy state
              END IF
            END DO
          END DO
          DPER(L)=NDES         !number of different vibrational energy levels
!
          A=0.d0
          DO I=1,NDES
            EVIB=TVAR(IDES(I))
            IF (EVIB < DABS(REA(5,KK))) A=A+1.d0  !number of different energy states bellow REA(5,K)
          END DO
          !this is the C constant in Bondar's paper required for the normalization condition (Sum[CPER/DPER]=A*CPER/DPER=1)
          !note that only states with Evib < REA(5,K) contribute to the summation; others states are prohibited
          CPER(KK,KN)=DPER(L)/A
        END IF
        !write(*,*) KK,L,IDL(:,L),II,DPER(L),A,CPER(KK,KN)
      END DO
    END IF
  END DO
!
END IF
!pause
!
IF (MSP == 1) THEN
  RMAS=SPM(1,1,1)
  CXSS=SPM(2,1,1)
  RGFS=SPM(6,1,1)
END IF
!
DO L=1,MSP
  CR(L)=0.
  DO M=1,MSP
    CR(L)=CR(L)+2.D00*SPI*SPM(4,L,M)**2*FND(1)*FSP(M,1)*(FTMP(1)/SPM(5,L,M))** &
        (1.-SPM(3,L,M))*DSQRT(2.*BOLTZ*SPM(5,L,M)/SPM(1,L,M))
  END DO
END DO
A=0.D00
FP=0.D00
DO L=1,MSP
  A=A+FSP(L,1)*CR(L)
  FP=FP+(2./SPI)*FSP(L,1)*VMP(L,1)/CR(L)
END DO
CTM=1.D00/A
WRITE (9,*) 'Approximate collision time in the stream is',CTM
!
WRITE (9,*) 'Approximate mean free path in the stream is',FP
!
!--set the initial time step
DTM=CTM*CPDTM
IF (DABS(VFX(1)) > 1.D-6) THEN
  A=(0.5D00*DDIV/VFX(1))*TPDTM
ELSE
  A=0.5D00*DDIV/VMPM
END IF
IF (IVB == 1) THEN
  B=0.25D00*DDIV/(ABS(VELOB)+VMPM)
  IF (B < A) A=B
END IF
IF (DTM > A) DTM=A
DTSAMP=SAMPRAT*DTM
DTOUT=OUTRAT*DTSAMP
TSAMP=0.D00
TOUT=0.D00
TPOUT=OUTRAT
NOUT=-1
ENTMASS=0.D00
!
WRITE (9,*) 'The initial value of the overall time step is',DTM
!
!--initialise cell quantities associated with collisions
!
DO N=1,NCCELLS
  CCELL(3,N)=DTM/2.D00
  CCELL(4,N)=2.D00*VMPM*SPM(2,1,1)
  CALL ZGF(RANF,IDT)
  CCELL(2,N)=RANF
  CCELL(5,N)=0.D00
END DO

!
!--set the entry quantities
!
DO K=1,2
  IF (ITYPE(K) == 0) THEN
    DO L=1,MSP
      IF (K == 1) SN=VFX(1)/VMP(L,1)
      IF (K == 2) SN=-VFX(2)/VMP(L,2)
      AA=SN
      A=1.D00+ERF(AA)
      BB=DEXP(-SN**2)
      ENTR(3,L,K)=SN
      ENTR(4,L,K)=SN+SQRT(SN**2+2.D00)
      ENTR(5,L,K)=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,K)))
      ENTR(6,L,K)=3.D00*VMP(L,K)
      B=BB+SPI*SN*A
      ENTR(1,L,K)=(FND(K)*FSP(L,K)*VMP(L,K))*B/(FNUM*2.D00*SPI)
      ENTR(2,L,K)=0.D00
    END DO
  END IF
END DO
!
!--initialize variables for updating boundary conditions
!
UVFX(:)=VFX(:)  !--isebasti
UFND(:)=FND(:)
UFTMP(:)=FTMP(:)
UVMP(:,:)=VMP(:,:)
!
!--Set the uniform stream
!
MNM=1.1D00*NMI
!
IF (MMVM > 0) THEN
  ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),     &
       IPVIB(0:MMVM,MNM),STAT=ERROR)
  IPVIB=0
ELSE
  IF (MMRM > 0) THEN
    ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),STAT=ERROR)
  ELSE
    ALLOCATE (PX(MNM),PTIM(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),STAT=ERROR)
  END IF
END IF
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR MOLECULE ARRAYS',ERROR
ENDIF
!
NM=0
!
IF ((IGS==1).OR.(IGS==3)) THEN
  WRITE (*,*) 'Setting the initial gas'
  DO L=1,MSP
    ROTE=0.
    DO K=1,ISECS+1
      IF (ISECS == 0) THEN         !--no secondary stream
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)
        XMIN=XB(1)
        XMAX=XB(2)
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (K == 1) THEN
          M=IDINT(DFLOAT(NMI)*((XS**JFX-XB(1)**JFX)*FND(1)/A)*FSP(L,1))
          XMIN=XB(1)
          XMAX=XS
        ELSE
          M=IDINT(DFLOAT(NMI)*((XB(2)**JFX-XS**JFX)*FND(2)/A)*FSP(L,2))
          XMIN=XS
          XMAX=XB(2)
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,VMP(L,K),IDT)
            END DO
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,FRTMP(K),PROT(NM),IDT)
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,FVTMP(K),IPVIB(J,NM),J,IDT)
                  IF (IGS==3) THEN
                    !forcing a specific vibrational distribution
                    CALL ZGF(RANF,IDT)
                    !IF (IPVIB(J,NM)==5  .AND. RANF <= 0.6d-0) IPVIB(J,NM)=10 !Tv at 2930K
                    !IF (IPVIB(J,NM)>=15 .AND. RANF <= 0.2d-0) IPVIB(J,NM)=10 !Tv at 2930K
                    !IF (IPVIB(J,NM)>=5) IPVIB(J,NM)=10 !Tv at 2820K
                    !IF (RANF <= 1d-6) IPVIB(J,NM)=59 !trick
                  END IF
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
    END DO
  END DO
END IF
!
!--isebasti: included IGS=2 option
IF (IGS==2) THEN
  WRITE (*,*) 'Setting the initial gas with given temp TDF(NC) and number density DNF(NC) distributions'
  DO L=1,MSP
    ROTE=0.
    DO NC=1,NCELLS
   !DO K=1,ISECS+1
      XMIN=CELL(2,NC)
      XMAX=CELL(3,NC)
      TDF=FTMP(1)+2.d3*DEXP(-((CELL(1,NC)-XS)/0.2d-3)**2.d0)
      DNF=FND(1)*(FTMP(1)/TDF)     !nden is Gaussian
      TDF=FTMP(1)                  !T is uniform
      IF (ISECS == 0) THEN         !--no secondary stream
        K=1
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)*(CELL(4,NC)/(XB(2)-XB(1)))*(DNF/FND(1))
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (CELL(1,NC) < XS) THEN
          K=1
          M=IDINT(DFLOAT(NMI)*((XMAX**JFX-XMIN**JFX)*DNF/A)*FSP(L,1))
        ELSE
          K=2
          M=IDINT(DFLOAT(NMI)*((XMAX**JFX-XMIN**JFX)*DNF/A)*FSP(L,2))
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            VMP2=DSQRT(2.D00*BOLTZ*TDF/SP(5,L))  !based on TDF
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,VMP2,IDT)
            END DO
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,TDF,PROT(NM),IDT)   !based on TDF
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,TDF,IPVIB(J,NM),J,IDT)    !based on TDF
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
   !END DO
    END DO
  END DO
!
  WRITE (9,*) 'DMOM',DMOM
  WRITE (9,*) 'DENG',DENG
END IF
!
!--isebasti: included IGS=4 option
IF (IGS==4) THEN
  WRITE (*,*) 'Setting the initial gas with bimodal energy distribution'
  DO L=1,MSP
    ROTE=0.
    DO K=1,ISECS+1
      IF (ISECS == 0) THEN         !--no secondary stream
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)
        XMIN=XB(1)
        XMAX=XB(2)
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (K == 1) THEN
          M=IDINT(DFLOAT(NMI)*((XS**JFX-XB(1)**JFX)*FND(1)/A)*FSP(L,1))
          XMIN=XB(1)
          XMAX=XS
        ELSE
          M=IDINT(DFLOAT(NMI)*((XB(2)**JFX-XS**JFX)*FND(2)/A)*FSP(L,2))
          XMIN=XS
          XMAX=XB(2)
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            !--for bimodal distribution
            CALL ZGF(RANF,IDT)
            IF (RANF > 0.5d0) THEN
              BTEMP=FTMP(1)   !based on translational reference temperature
            ELSE
              BTEMP=FRTMP(1)  !based on rotational reference temperature
            END IF
            BVMP=SQRT(2.D00*BOLTZ*BTEMP/SP(5,L))
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,BVMP,IDT) !for bimodal distribution
            END DO
            !--end
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,BTEMP,PROT(NM),IDT) !for bimodal distribution
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,FVTMP(K),IPVIB(J,NM),J,IDT) !for Boltzmann distribution
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
    END DO
  END DO
END IF
!
!--SPECIAL CODING FOR INITIATION OF EXPLOSION IN H2-02 MIXTURE (FORCED IGNITION CASES)
IF ((GASCODE == 8).AND.(IFI > 0)) THEN !--isebasti: IFC>0 allow forced ignition
  IF(IGS == 1) A=0.5D00 !0.5 for p=0.1 atm;
  IF(IGS == 2) A=1.0D00
  !
  !--set the vibrational levels of A% random molecules to 5
  IF (IFI == 1) THEN
    M=0.01D00*A*NM
    DO N=1,M
      CALL ZGF(RANF,IDT)
      K=INT(RANF*DFLOAT(NM))+1
      L=IPSP(K)       !--isebasti: modified to increase all vibrational degrees of freedom
      DO KK=1,ISPV(L) !--isebasti
        IPVIB(KK,K)=5 !--isebasti
      END DO          !--isebasti
    END DO
  END IF
  !
  !--isebasti: set vibrational levels of molecules within XMIN and XMAX to dissociation level
  IF (IFI == 2) THEN
    XMIN=XS-0.005D00*A*(XB(2)-XB(1))
    XMAX=XS+0.005D00*A*(XB(2)-XB(1))
    DO N=1,NM
      IF ((PX(N) >= XMIN).AND.(PX(N) <= XMAX)) THEN
        L=IPSP(N)
        IF (ISPV(L) > 0) IPVIB(1,N)=IDL(1,L) !setting only one mode
      END IF
    END DO
  END IF
  !
  !--isebasti: set molecules within XMIN and XMAX to energy states consistent with temperature TDF
  IF (IFI == 3) THEN
    XMIN=XS-0.005D00*A*(XB(2)-XB(1))
    XMAX=XS+0.005D00*A*(XB(2)-XB(1))
    DO N=1,NM
      IF ((PX(N) >= XMIN).AND.(PX(N) <= XMAX)) THEN
        L=IPSP(N)
        TDF=5.D03
        B=SQRT(2.D00*BOLTZ*TDF/SP(5,L))  !B=VMP at temperature TDF
        DO K=1,3
          CALL RVELC(PV(K,N),BB,B,IDT)
        END DO
        IF (PX(N) < XS) PV(1,N)=PV(1,N)+VFX(1)
        IF (PX(N) > XS) PV(1,N)=PV(1,N)+VFX(2)
        IF (ISPR(1,L) > 0) CALL SROT(L,TDF,PROT(N),IDT)
        IF (MMVM > 0) THEN
          DO K=1,ISPV(L)
            CALL SVIB(L,TDF,IPVIB(K,N),K,IDT)
          END DO
        END IF
      END IF
    END DO
END IF
!
END IF

CLOSE(9)
!
!
RETURN
END SUBROUTINE SET_INITIAL_STATE
