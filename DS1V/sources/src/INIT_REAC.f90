!
!*****************************************************************************
!
SUBROUTINE INIT_REAC
!
!--initialise the reaction variables (in the old reaction model)
!
USE MOLECS
!
USE GAS
USE CALC
USE GEOM
USE OUTPUT
USE MFDSMC, only : IMF,IMFS,MFRMASS, NMFANG,NMFpair, IMFpair,MF_SET_AHO, IMFdia
!
IMPLICIT NONE
!
INTEGER :: I,K,L,M,N,LS,MS,NNRE,NTBP,II
REAL(KIND=8) :: B,C,D,EPS,AL,X,BB
REAL(8),EXTERNAL :: GAM
!
!--I,K,L,M,N working integers
!--A,B working variables
!--LS,MS species
!--EPS symmetry factor
!
!
!--convert the reaction data into the required form
!
IF (MNRE == 0) THEN
  NNRE=1
ELSE
  NNRE=MNRE
END IF
IF (MTBP == 0) THEN
  NTBP=1
ELSE
  NTBP=MTBP
END IF
!
IF (MNRE > 0) THEN
  REA=0.; IREA=0; NREA=0; JREA=0; NRSP=0; IRCD=0 ; REAC=0.; EVREM = 0.d0
  IF ( IMF .ne. 0) THEN
     MFRMASS=-100.d0; NMFANG = 0;
  END IF
!
  OPEN(10, FILE="ChemicalReaction.TXT")
  WRITE(10, "(A, I3)") "Number of reactions",MNRE
  DO  N=1,MNRE
    WRITE(10, "(A,I3)") "Reaction ",N
    L=LE(N)
    M=ME(N)
!--the pre-reaction species codes are L and M
    IF ((L > 0).AND.(M > 0)) THEN
      IREA(1,N)=L
      IREA(2,N)=M
      IF (M < L) THEN
        K = L
        L = M
        M = K
      END IF
      NRSP(L,M) = NRSP(L,M) + 1
      IF (L /= M) NRSP(M,L) = NRSP(M,L) + 1
      K = NRSP(L,M)
!--this is the Kth reaction involving species L and M
      IRCD(K,M,L)=N
      IRCD(K,L,M)=N
    END IF

!============= set post reaction info
    K=KP(N)
    L=LP(N)
    M=MP(N)
!--three post-collision particles (dissociation)
    IF (K > 0) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Dissociation'
      WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3,' + ',I3)") LE(N),ME(N),K,L,M
      WRITE(10,*)
      NREA(1,N)=2
      NREA(2,N)=1
      JREA(1,N,1)=K
      JREA(1,N,2)=L
      JREA(2,N,1)=M
    END IF
!--two post-collision particles (exchange)
    IF (K == 0) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Exchange'
      WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3)") LE(N),ME(N),L,M
      WRITE(10,*)
      NREA(1,N)=1
      NREA(2,N)=1
      JREA(1,N,1)=L
      JREA(2,N,1)=M
    END IF
!--one post collision particle (recombination)
    IF (K == -1) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Recombination'
      WRITE(10, "(I3,' + ',I3,' + ',I3,' -> ',I3,' + ',I3)") LE(N),ME(N),M,L,M
      WRITE(10,*)
      NREA(1,N)=1
      NREA(2,N)=0
      JREA(1,N,1)=L
      JREA(2,N,1)=M
    END IF
!
    REA(1,N)=CI(N)
    REA(2,N)=AE(N)
    REA(3,N)=AC(N)
    REA(4,N)=BC(N)
    REA(5,N)=ER(N)
  END DO
  CLOSE(10)
!
END IF
!
DO L=1,MNRE
  DO M=1,MSP
    DO N=1,MSP
      IF (IRCD(L,N,M) == 0) IRCD(L,N,M)=IRCD(L,M,N)
    END DO
  END DO
END DO
!
!--calculate the coefficients of the energy terms in steric factors (see SMILE documentation)
!
NMFpair = 0
IF (IMF .ne. 0 .and. IMFS == 1 .and. MNRE >0) THEN
  DO I=1,MSP
    DO K = I,MSP
      NMFpair = NMFpair+1
      IMFpair(I,K) = NMFpair
      IMFpair(K,I) = NMFpair
    END DO
  END DO
END IF


IF (IMF >= 2 .and. MNRE>0) THEN
  CALL MF_SET_AHO()
END IF

IF (IMF .ne. 0) THEN
  OPEN(10, FILE="ChemicalReaction.TXT", ACCESS="APPEND")
  WRITE(10,*)
  WRITE(10,*)
  IF (IMF == 1) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-SHO  Dia:',IMFdia
  ELSE IF (IMF == 2) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-AHO QCT vphase  Dia:',IMFdia
  ELSE IF (IMF == 3 .or. IMF == 4) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-AHO Morse vphase  Dia:',IMFdia
  END IF

  IF (IMF == 4) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: Impact parameter coulpling'
  END IF

  WRITE(10,'(A)') "Macheret-Fridman model parameter:"
END IF

DO K=1,MNRE
  LS=IREA(1,K)
  MS=IREA(2,K)
  EPS=1.D00
  IF (LS == MS) EPS=2.D00
  I=NREA(1,K)+NREA(2,K)
  AL=SPM(3,LS,MS)-0.5d0  !alpha
  X=REA(4,K)-0.5d0+AL
!
  IF (I == 1) THEN
!--recombination reaction
    C=EPS*SPI*REA(3,K)*SPM(1,LS,MS)**0.5d0
    D=(2.d0**1.5d0)*(BOLTZ**REA(4,K))*SPM(2,LS,MS)*((2.d0-AL)*BOLTZ*SPM(5,LS,MS))**AL
    REA(6,K)=(C/D)
  END IF
  IF (I == 2) THEN
!--exchange reaction
    BB=2.d0-AL
    B=GAM(BB)
    C=EPS*SPI*REA(3,K)*(SPM(1,LS,MS)/(8.d0*BOLTZ))**0.5d0
!   D=SPM(2,LS,MS)*(BOLTZ**X)*(2.d0-AL*SPM(5,LS,MS))**AL  !Alexeenko-SMILE report with bug
!   REA(6,K)=(C/D)*(1.d0/B)                               !Alexeenko-SMILE report with bug
    D=SPM(2,LS,MS)*(BOLTZ**X)*(SPM(5,LS,MS))**AL          !Ivanov-SMILE report
    REA(6,K)=(C/D)                                        !Ivanov-report
  END IF
!--dissociation reaction
  IF (I == 3) THEN
    IF (REA(1,K) < 0.) THEN !TCE
      BB=2.d0-AL
      B=GAM(BB)
      C=EPS*SPI*REA(3,K)*(SPM(1,LS,MS)/(8.d0*BOLTZ))**0.5d0
      D=SPM(2,LS,MS)*(BOLTZ**X)*(SPM(5,LS,MS))**AL        !Ivanov-SMILE report
      REA(6,K)=(C/D)                                      !Ivanov-SMILE report
    ELSE  !VDC
      BB=2.d0-AL
      B=GAM(BB)
      C=EPS*REA(3,K)*SPI*0.5d0*(SPM(1,LS,MS)/(2.d0*BOLTZ))**(0.5d0-AL)
      D=SPM(2,LS,MS)*B
      REA(6,K)=(C/D)
    END IF

    IF (IMF .ne. 0) THEN  !pre set for Macheret Fridman model
      IF (GASCODE .ne. 8)THEN
        STOP "Only gascode = 8 works with MF model"
      END IF
      ! LS is always the one to be dissociated
      ! MFRAMSS(1): reduced mass of whole system
      IF ((LS .gt. 5) .or. (MS .gt. 5)) THEN
        WRITE(10, "(A, I3, A, I3,A)") "Reaction: ",LS, ' + ',MS, ' is not supported by MF model'
        WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3,' + ',I3)") LE(K),ME(K),KP(K),LP(K),MP(K)
        WRITE(10,*)
        MFRMASS(:,K) = 0.0d0
        NMFANG(K) = 0
      ELSE
      !----- calculated reduced mass
      ! reduced mass of A+B
        MFRMASS(1,K) = SP(5,LS)*SP(5,MS)/(SP(5,LS)+SP(5,MS))
        DO II = 1,2
          IF (ISPV(IREA(II,K)) == 0) THEN
          ! single atom
            MFRMASS(II+1,K) = SP(5,IREA(II,K))
          ELSE IF (IREA(II,K) == 5) THEN
          ! NO molecule
            MFRMASS(II+1,K) = SP(5,2)*SP(5,4)/(SP(5,2)+SP(5,4))
          ELSE IF (IREA(II,K) == 3 .or. IREA(II,K) == 1) THEN
          ! nitogen, oxygen and other homonuclear molecule
          ! reduced mass should be 1/4 of the mass of molecule
            MFRMASS(II+1,K) = SP(5,IREA(II,K))/4.0d0
          ELSE
            WRITE(*,*) "Not supported MFRMASS", IREA(II,K)
            STOP
          END IF
        END DO
        WRITE(10, "(A, I3, A, I3)") "Reaction: ",LS, ' + ',MS
        WRITE(10, "(A, 3(G10.3, 1X))") "Reduced mass (kg): ", MFRMASS(1,K), MFRMASS(2,K),MFRMASS(3,K)
      !
      !----- determine reaction type
        NMFANG(K) = 4 ! default: homo-atom
        IF (ISPV(MS) == 1) NMFANG(K) = 8  !default: homo-homo
        IF (LS == 5) NMFANG(K) = NMFANG(K)+1 !hetero-homo or hetero-atom: 5,9
        IF (MS == 5) NMFANG(K) = NMFANG(K)+2
      ! ----- homo   hetero   !dissociated molecule
      !atom    4       5
      !homo    8       9
      !hetero  10     11
      !
      END IF
    END IF
  END IF
!
END DO
IF (IMF .ne. 0) THEN
  CLOSE(10)
END IF
!
RETURN
!
END SUBROUTINE INIT_REAC
