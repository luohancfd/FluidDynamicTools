!
!*****************************************************************************
SUBROUTINE OUTPUT_RESULTS
!
!--calculate the surface and flowfield properties
!--generate TECPLOT files for displaying these properties
!--calculate collisiion rates and flow transit times and reset time intervals
!--add molecules to any flow plane molecule output files
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE MFDSMC
!
IMPLICIT NONE
!
INTEGER :: I,IJ,J,JJ,K,KK,JR,KR,L,M,N,NN,NMCR,CTIME,NMS(MSP),ISNAP(1,NSNAP),ZCHECK,NSS,MCT,MAXLEVEL,IDT=0 !--isebasti: included ZCHECK,NSS,IDT
INTEGER(KIND=8) :: NNN
REAL(KIND=8) :: AA,BB,BB2,CC,AS,AT,AU,AQ  !--isebasti: included RANF,GAM
REAL(KIND=8) :: A,B,C,C2,D,SDTM,SMCR,DOF,AVW,UU,VDOFM,TVIBM,DTMI,TT,EVIBT,RANF,F(MMVM,0:100,MSP),TVPDF(MMVM,MSP)  !--isebasti: included RANF,F,TVPDF
REAL(KIND=8) :: VPF(MMVM,MSP),EQROT(MSP,100),EQVIBVT(MMVM,MSP,0:100),EQVIBOT(MMVM,MSP,0:100)  !--isebsati: included
REAL(KIND=8) :: DSUM(0:14),SUMS(0:8,2),TV(MMVM,MSP),TVIB(MSP),DF(NCELLS,MMVM,MSP),VDOF(MSP),PPA(MSP),&
                THCOL(MSP,MSP),PSNAP(1,NSNAP),VSNAP(3,NSNAP),SDOF(MSP),SRR(MNRE),EVIB,QVIBVT,QVIBOT,COEF(0:9),ROOTF
INTEGER,PARAMETER :: PQKIND = SELECTED_REAL_KIND(8)
REAL(KIND=PQKIND) :: AQUAD
CHARACTER (LEN=32) :: FILENAME,TNAME
CHARACTER (LEN=4) :: E
REAL(8) :: COLRR(MNRE)
REAL(8),EXTERNAL :: ERF,gamain,GAM
!
!--CTIME  computer time (microseconds)
!--SUMS(N,L) sum over species of CSS(N,J,L,M) for surface properties
!
!--For flowfield properties,where <> indicates sampled sum
!--DSUM(0) the molecular number sum over all species
!--DSUM(1) the weighted number sum over all species
!--DSUM(2) the weighted sum of molecular masses
!--DSUM(3),(4),(5) the weighted sum over species of m*<u>,<v>,<w>
!--DSUM(6) the weighted sum over species of m*(<u**2>+<v**2>+<w**2>)
!--DSUM(7) the weighted sum over species of <u**2>+<v**2>+<w**2>
!--DSUM(8) the weighted sum of rotational energy
!--DSUM(9) the weighted sum of rotational degrees of freedom
!--DSUM(10) the weighted sum over species of m*<u**2>
!--DSUM(11) the weighted sum over species of m*<v**2>
!--DSUM(12) sum over species of m*<w**2>
!--DSUM(13) the weighted sum over species of m*u*v
!--DSUM(14) the weighted sum over species of m*(u^2+v^2+w^2)*u
!--UU velocity squared
!--DOF degrees of freedom
!--AVW the average value of the viscosity-temperature exponent
!--DVEL velocity difference
!--TVEL thermal speed
!--SMCR sum of mcs/mfp over cells
!--NMCR number in the sum
!--VDOFM effective vibrational degrees of freedom of mixture
!--TVIB(L)
!--VDOF(L)
!--TV(K,L) the temperature of vibrational mode K of species L
!--SDOF(L) total degrees of freedom for species L
!--PPA particles per atom
!--NMS number per species
!--NSNAP number of particles considered in SNAPSHOT files
!--SRR sampled reaction rate
!
!----------------------------------------------------------------------------
!--set variables
!
NOUT=NOUT+1
IF (NOUT > 9999) NOUT=NOUT-9999
CALL NUMCHAR4 (NOUT,E)
WRITE (*,*) 'Generating files for output interval',NOUT
WRITE (*,*) 'ISF,FTIME,Number of samples',ISF,FTIME,NSAMP
!
!----------------------------------------------------------------------------
!--compute surface properties
!
VARSP=0.D00
IF (IFX == 0) A=FNUM/(FTIME-TISAMP)   !--flow X-section area = unity for 1-D flow
DO JJ=1,2
  IF (IFX == 1) A=FNUM/(2.D00*PI*XB(JJ)*(FTIME-TISAMP))
  IF (IFX == 2) A=FNUM/(4.D00*PI*XB(JJ)*XB(JJ)*(FTIME-TISAMP))
  !--JJ=1 for surface at XB(1), JJ=2 for surface at XB(2)
  IF (ITYPE(JJ) == 2) THEN
    SUMS=0.D00
    DO L=1,MSP
      DO J=0,8
        DO IJ=1,2
          SUMS(J,IJ)=SUMS(J,IJ)+CSS(J,JJ,L,IJ)
        END DO
      END DO
    END DO
!
    VARS(0,JJ)=SUMS(0,1)
    VARS(1,JJ)=SUMS(1,1)
    VARS(2,JJ)=SUMS(1,2)
    VARS(3,JJ)=SUMS(1,1)*A
    VARS(4,JJ)=SUMS(1,2)*A
    VARS(5,JJ)=SUMS(2,1)*A
    VARS(6,JJ)=SUMS(2,2)*A
    VARS(7,JJ)=SUMS(3,1)*A
    VARS(8,JJ)=SUMS(3,2)*A
    VARS(9,JJ)=SUMS(4,1)*A
    VARS(10,JJ)=SUMS(4,2)*A
    VARS(11,JJ)=SUMS(5,1)*A
    VARS(12,JJ)=SUMS(5,2)*A
    VARS(13,JJ)=SUMS(6,1)*A
    VARS(14,JJ)=SUMS(6,2)*A
    VARS(15,JJ)=SUMS(7,1)*A
    VARS(16,JJ)=SUMS(7,2)*A
    VARS(17,JJ)=SUMS(8,1)*A
    VARS(18,JJ)=SUMS(8,2)*A
    IF (CSSS(1,JJ) > 1.D-6) THEN
      VARS(19,JJ)=CSSS(3,JJ)/CSSS(2,JJ)
      VARS(20,JJ)=(CSSS(4,JJ)-CSSS(2,JJ)*VARS(19,JJ)*VARS(19,JJ))/(CSSS(1,JJ)*3.D00*BOLTZ)-TSURF(JJ)
      VARS(19,JJ)=VARS(19,JJ)-VSURF(JJ)
      IF (CSSS(6,JJ) > 0.D00) THEN
        VARS(21,JJ)=(2.D000/BOLTZ)*(CSSS(5,JJ)/CSSS(6,JJ))-TSURF(JJ)
      ELSE
        VARS(21,JJ)=0.D00
      END IF
    ELSE
      VARS(19,JJ)=0.D00
      VARS(20,JJ)=0.D00
      VARS(21,JJ)=0.D00
    END IF
    VARS(22,JJ)=(SUMS(2,1)+SUMS(2,2))*A
    VARS(23,JJ)=(SUMS(3,1)+SUMS(3,2))*A
    VARS(24,JJ)=(SUMS(4,1)+SUMS(4,2))*A
    VARS(25,JJ)=(SUMS(5,1)+SUMS(5,2))*A
    VARS(26,JJ)=(SUMS(6,1)+SUMS(6,2))*A
    VARS(27,JJ)=(SUMS(7,1)+SUMS(7,2))*A
    VARS(28,JJ)=(SUMS(8,1)+SUMS(8,2))*A
    VARS(29,JJ)=VARS(11,JJ)+VARS(13,JJ)+VARS(15,JJ)
    VARS(30,JJ)=VARS(12,JJ)+VARS(14,JJ)+VARS(16,JJ)
    VARS(31,JJ)=VARS(29,JJ)+VARS(30,JJ)
    DO L=1,MSP
      IF (SUMS(1,1) > 0) THEN
        VARS(32+L,JJ)=100.*CSS(1,JJ,L,1)/SUMS(1,1)
      ELSE
        VARS(32+L,JJ)=0.
      END IF
    END DO
  END IF
END DO
!
!----------------------------------------------------------------------------
!--compute flowfield properties
VAR=0.D00
VARSP=0.
SMCR=0
NMCR=0
VDOFM=0.
DO N=1,NCELLS
  A=FNUM/(CELL(4,N)*NSAMP)
  IF (IVB == 1) A=A*((XB(2)-XB(1))/(XB(2)+VELOB*0.5D00*(FTIME+TISAMP)-XB(1)))**(IFX+1)
!--check the above for non-zero XB(1)
  DSUM=0.
  NMCR=NMCR+1
  DO L=1,MSP
    IF (CS(0,N,L) > 0.99d0) THEN
      DSUM(0)=DSUM(0)+CS(0,N,L)
      DSUM(1)=DSUM(1)+CS(1,N,L)
      DSUM(2)=DSUM(2)+SP(5,L)*CS(1,N,L)
      DO K=1,3
        DSUM(K+2)=DSUM(K+2)+SP(5,L)*CS(K+1,N,L)
        IF (CS(1,N,L) >= 0.1D00) THEN
          VARSP(K+1,N,L)=CS(K+4,N,L)/CS(1,N,L)
  !--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
          VARSP(K+8,N,L)=CS(K+1,N,L)/CS(1,N,L)
  !--VARSP(9,10,11 are temporarily the mean of the velocities
        END IF
      END DO
      DSUM(6)=DSUM(6)+SP(5,L)*(CS(5,N,L)+CS(6,N,L)+CS(7,N,L))
      DSUM(7)=DSUM(7)+CS(5,N,L)+CS(6,N,L)+CS(7,N,L)
      IF (ISPR(1,L) > 0) THEN
        DSUM(8)=DSUM(8)+CS(8,N,L)
        DSUM(9)=DSUM(9)+CS(1,N,L)*ISPR(1,L)
      END IF
      DSUM(10)=DSUM(10)+SP(5,L)*CS(5,N,L)
      DSUM(11)=DSUM(11)+SP(5,L)*CS(6,N,L)
      DSUM(12)=DSUM(12)+SP(5,L)*CS(7,N,L)

      DSUM(13)=DSUM(13)+SP(5,L)*CSH(1,N,L)
      DSUM(14)=DSUM(14)+SP(5,L)*CSH(2,N,L)
    END IF
  END DO
  AVW=0.
  DO L=1,MSP
    IF (CS(0,N,L) > 0.99D0) THEN
      VARSP(0,N,L)=CS(1,N,L)
      VARSP(1,N,L)=0.D00
      VARSP(6,N,L)=0.
      VARSP(7,N,L)=0.
      VARSP(8,N,L)=0.
      IF (DSUM(1) > 0.1) THEN
        VARSP(1,N,L)=CS(1,N,L)/DSUM(1)  !isebasti: deleted 100* factor
        AVW=AVW+SP(3,L)*CS(1,N,L)/DSUM(1)
        IF ((ISPR(1,L) > 0).AND.(CS(1,N,L) > 0.5)) VARSP(6,N,L)=(2.D00/BOLTZ)*CS(8,N,L)/(DFLOAT(ISPR(1,L))*CS(1,N,L))
      END IF
      VARSP(5,N,L)=0.
      DO K=1,3
        VARSP(K+1,N,L)=(SP(5,L)/BOLTZ)*(VARSP(K+1,N,L)-VARSP(K+8,N,L)**2)
        VARSP(5,N,L)=VARSP(5,N,L)+VARSP(K+1,N,L)
      END DO
      VARSP(5,N,L)=VARSP(5,N,L)/3.D00
      VARSP(8,N,L)=(3.D00*VARSP(5,N,L)+DFLOAT(ISPR(1,L))*VARSP(6,N,L))/(3.D00+DFLOAT(ISPR(1,L))) !isebasti: included according to DSMC.f90
    END IF
  END DO
!
  IF (IVB == 0) VAR(1,N)=CELL(1,N)
  IF (IVB == 1) THEN
    C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV)      !--new DDIV
    VAR(1,N)=XB(1)+(DFLOAT(N-1)+0.5)*C
  END IF
  VAR(2,N)=DSUM(0)
  IF (DSUM(1) > 0.5) THEN
    VAR(3,N)=DSUM(1)*A    !--number density Eqn. (4.28)
    VAR(4,N)=VAR(3,N)*DSUM(2)/DSUM(1)   !--density Eqn. (4.29)
    VAR(5,N)=DSUM(3)/DSUM(2)    !--u velocity component Eqn. (4.30)
    VAR(6,N)=DSUM(4)/DSUM(2)    !--v velocity component
    VAR(7,N)=DSUM(5)/DSUM(2)    !--w velocity component
    UU= VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2

    VAR(22,N)=-A*DSUM(13)  ! shear tau_xy Eqn(1.52)
    ! for 1D case u*v should be zero, i.e. no motion in one direction
    ! Full expression:
    !VAR(22,N)=-A*(DSUM(13) - DSUM(2)*VAR(5,N)*VAR(4,N))
    ! tau = - n(\bar{muv} - \bar{m} u0*v0)

    VAR(23,N)= A*(0.5d0*DSUM(14) - 0.5d0*VAR(5,N)*DSUM(6) &
                  & - VAR(5,N)*DSUM(10) - VAR(6,N)*DSUM(13) &
                  & + VAR(5,N)**3*DSUM(2))
    ! VAR(23,N)= A*(0.5d0*DSUM(14) - 0.5d0*VAR(5,N)*DSUM(6) &
    !              & - VAR(5,N)*DSUM(10) - VAR(6,N)*DSUM(13) &
    !              & + VAR(5,N)**3*DSUM(2) + VAR(5,N)*VAR(6,N)**2*DSUM(2))
    ! Tensor notation
    ! Denote operator Z as the summation over all molecule,
    ! V_j = Z(m v_j) / Z(m)
    ! q_i = A * Z( 1/2 * m (v_j-V_j)*(v_j-V_j)*(v_i - V_i))
    !     = A * { 1/2*Z(m*v_i*v_j*v_j) - 1/2*V_i*Z(m*v_j*v_j) - V_j*Z(m*v_i*v_j) + V_i*V_j*V_j*Z(m)
    ! For 1D case, W_0 = 0
    ! For Couette flow, VAR(5,N) = 0, then you can get Eqn 12.14 in Bird's book
    IF (DSUM(1) > 1) THEN
      VAR(8,N)=(ABS(DSUM(6)-DSUM(2)*UU))/(3.D00*BOLTZ*DSUM(1))  !Eqn. (1.51)
!--translational temperature
      VAR(19,N)=(ABS(DSUM(10)-DSUM(2)*VAR(5,N)**2))/(BOLTZ*DSUM(1))
      VAR(20,N)=(ABS(DSUM(11)-DSUM(2)*VAR(6,N)**2))/(BOLTZ*DSUM(1))
      VAR(21,N)=(ABS(DSUM(12)-DSUM(2)*VAR(7,N)**2))/(BOLTZ*DSUM(1))
    ELSE
      VAR(8,N)=1.
      VAR(19,N)=1.
      VAR(20,N)=1.
      VAR(21,N)=1.
    END IF
!--rotational temperature
    IF (DSUM(9) > 0.01D00) THEN
      VAR(9,N)=(2.D00/BOLTZ)*DSUM(8)/DSUM(9)    !Eqn. (4.36)
    ELSE
      VAR(9,N)=0.
    END IF
    DOF=(3.D00+DSUM(9)/DSUM(1))
!--vibration temperature default
    VAR(10,N)=FTMP(1)
!--overall temperature based on translation and rotation
    VAR(11,N)=(3.*VAR(8,N)+(DSUM(9)/DSUM(1))*VAR(9,N))/DOF
!--scalar pressure (now (from V3) based on the translational temperature)
    VAR(18,N)=VAR(3,N)*BOLTZ*VAR(8,N)
!
!--Tvib calculations according to DSMC.f90
    IF (MMVM > 0) THEN
      TV = 0.0d0
      TVIB=0.0d0
      DO L=1,MSP
        VDOF(L)=0.
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (CS(K+8,N,L) > 0.d0 .and. CS(0,N,L) >= 1.0d0) THEN
              AQUAD=CS(K+8,N,L)/CS(1,N,L)   !average vibrational energy (J)
              EVIB=AQUAD/EVOLT              !average vibrational energy (eV)
              IF (IVMODEL(L,1) == 0) THEN
                TV(K,L)=SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/AQUAD)  !--Eqn.(4.45) - assuming SHO
              ELSE
!do i=0,100
!aquad=VIBEN(0,L)*.9999+(0.09745000001d0*EVOLT-VIBEN(0,L))*DFLOAT(i)/100.d0
                  IF (AQUAD >= VIBEN(0,L)) THEN
                    !use secant method to calculate Tvib(evib) based on Boltzmann distribution
                    CALL ROOTF_SECANT(1,L,L,L,1.d3,AQUAD,TV(K,L))
                    !DF(N,K,L)=2.d0*AQUAD/(BOLTZ*TV(K,L)) !--Eqn. (11.28) Bird94 - general definition
                    DF(N,K,L)=2.d0*(SPVM(1,K,L)/TV(K,L))/(DEXP(SPVM(1,K,L)/TV(K,L))-1.d0) !SHO model
                  ELSE
                    TV(K,L)=1.d0 !to avoid division by zero
                    DF(N,K,L)=0.d0
                  END IF
!write(*,*) AQUAD/EVOLT,TV(K,L),SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/DFLOAT(AQUAD))
!pause
!end do
!stop
              END IF
            ELSE
              TV(K,L)=0.
              DF(N,K,L)=0.
            END IF
            VDOF(L)=VDOF(L)+DF(N,K,L)  !--Eqn.(4.49)
          END DO
          TVIB(L)=0.
          DO K=1,ISPV(L)
            IF (VDOF(L) > 1.D-6) THEN
              TVIB(L)=TVIB(L)+TV(K,L)*DF(N,K,L)/VDOF(L)  !--Eqn.(4.50)
            ELSE
              TVIB(L)=SUM(TV(:,L))/DFLOAT(ISPV(L))
            END IF
          END DO
        ELSE
          TVIB(L)=0. !TREF  !--isebasti: TREF is not defined
          VDOF(L)=0.
        END IF
        VARSP(7,N,L)=TVIB(L)
      END DO
      VDOFM=0.
      TVIBM=0.
      A=1.D00 !--isebasti: instead of 0
      DO L=1,MSP
        IF (ISPV(L) > 0) A=A+CS(1,N,L)
      END DO
      DO L=1,MSP
        IF (ISPV(L) > 0) THEN
          VDOFM=VDOFM+VDOF(L)*CS(1,N,L)/A  !--Eqn.(4.51)
          TVIBM=TVIBM+TVIB(L)*CS(1,N,L)/A  !--Eqn.(4.52)
        END IF
      END DO
      VAR(10,N)=TVIBM
    END IF
!
!--convert the species velocity components to diffusion velocities
    DO L=1,MSP
      IF (VARSP(0,N,L) > 0.5) THEN
        DO K=1,3
          VARSP(K+8,N,L)=VARSP(K+8,N,L)-VAR(K+4,N)
        END DO
      ELSE
        DO K=1,3
          VARSP(K+8,N,L)=0.D00
        END DO
      END IF
    END DO
!
!--reset the overall temperature and degrees of freedom (now including vibrational modes)
    IF (MMVM > 0) THEN
      DO L=1,MSP
        SDOF(L)=3.D00+ISPR(1,L)+VDOF(L)
        VARSP(8,N,L)=(3.*VARSP(5,N,L)+ISPR(1,L)*VARSP(6,N,L)+VDOF(L)*VARSP(7,N,L))/SDOF(L)  !species overall T
        ! Not verified yet
        IF (VARSP(0,N,L) > 0.5) THEN
          VARSP(12,N,L)= (CSH(3,N,L) - CS(8,N,L)*VAR(5,N))*FNUM/(CELL(4,N)*NSAMP)
          VARSP(13,N,L)= (CSH(5,N,L) - CS(4,N,L)*VAR(5,N))*FNUM/(CELL(4,N)*NSAMP)
        END IF
      END DO
      A=0.D00
      B=0.D00
      DO L=1,MSP
        A=A+SDOF(L)*VARSP(8,N,L)*CS(1,N,L)
        B=B+SDOF(L)*CS(1,N,L)
      END DO
      VAR(11,N)=A/B !mixture overall T
      DOF=DOF+VDOFM !--isebasti: included
    END IF
!
!--Mach number
    VAR(17,N)=DSQRT(VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2)
    VAR(12,N)=VAR(17,N)/SQRT((DOF+2.D00)*VAR(11,N)*(DSUM(1)*BOLTZ/DSUM(2))/DOF)
!--average number of molecules in (collision) cell
    VAR(13,N)=DSUM(0)/DBLE(NSAMP)/DFLOAT(NCIS)
    IF (COLLS(N) > 2.) THEN
!--mean collision time (see page 17 of my Bird94 book)
      VAR(14,N)=0.5D00*(FTIME-TISAMP)*(DSUM(1)/NSAMP)/WCOLLS(N)
!--mean free path (based on r.m.s speed with correction factor based on equilibrium)
      VAR(15,N)=0.92132D00*DSQRT(DABS(DSUM(7)/DSUM(1)-UU))*VAR(14,N)
      VAR(16,N)=CLSEP(N)/(COLLS(N)*VAR(15,N))
    ELSE
!--m.f.p set by nominal values
      VAR(14,N)=1.D10
      VAR(15,N)=1.D10/VAR(3,N)
    END IF
  ELSE
    DO L=3,19
      VAR(L,N)=0.
    END DO
  END IF
!
  IF(N == NSPDF) THEN
    TVPDF(:,:)=TV(:,:)                        !store TV for the cell with pdf sampling
    IF (DSUM(1) > 0.5d0 )THEN
      A = 2.*PCOLLS*DFLOAT(NSAMP)/DSUM(1)         !mean collision times
    ELSE
      A = 0.0d0
    END IF
    MCT=A                                     !integer value
!
    IF(VARSP(1,N,1) > VARSP(1,N,3)) THEN
      L=1                                     !for N2 species
    ELSE
      L=3                                     !for O2 species
    END IF
!
    B=0.d0
    QVIBVT=0.d0
    DO M=0,IVMODEL(L,2)
      CALL VIB_ENERGY(EVIBT,M,1,L)
      B=B+EVIBT*DEXP(-EVIBT/(BOLTZ*VAR(8,N)))
      QVIBVT=QVIBVT+DEXP(-EVIBT/(BOLTZ*VAR(8,N)))
    END DO
    EVIBT=(B/QVIBVT)/EVOLT                     !equilibrium vibrational temperature based Tt (ev)
!
    OPEN (7,FILE='RELAX.DAT',ACCESS='APPEND') !write temperature relaxation
    IF (VARSP(1,N,4)>0.0d0) THEN
      B = 1.d0/(VAR(3,N)*VARSP(1,N,4)*CSCR(L,4)/TCOL(L,4)) !mct between collisions of one particular L molecule with any O
    ELSE
      B = -1.0d0
    END IF
    IF (VARSP(1,N,L)>0.0d0) THEN
      C = 1.d0/(VAR(3,N)*VARSP(1,N,L)*CSCR(L,L)/TCOL(L,L)) !mct between collisions of one particular L molecule with any L
    ELSE
      C = -1.0d0
    END IF
    IF (DSUM(1) > 0.5d0 ) THEN
      WRITE (7,993) FTIME,VAR(8:11,N),VAR(18,N),DSUM(9)/DSUM(1),VDOFM,NM,A,&
        DABS(0.5*(VAR(8,N)+VAR(9,N))-VAR(10,N)),DABS(VAR(8,N)-VAR(9,N)),DABS(VAR(8,N)-VAR(10,N)),&
        (TV(:,L),L=1,MSP),EVIB,EVIBT,VAR(14,N),B,C
    ELSE
      WRITE (7,993) FTIME,VAR(8:11,N),VAR(18,N),0.0d0,VDOFM,NM,A,&
        DABS(0.5*(VAR(8,N)+VAR(9,N))-VAR(10,N)),DABS(VAR(8,N)-VAR(9,N)),DABS(VAR(8,N)-VAR(10,N)),&
        (TV(:,L),L=1,MSP),EVIB,EVIBT,VAR(14,N),B,C
    END IF
    CLOSE (7)
  END IF
!
END DO
!
!----------------------------------------------------------------------------
!--write general/surface properties
!
OPEN (3,FILE='DS1GEN.DAT')
!
IF (IFX == 0) WRITE (3,*) 'DSMC program DS1 for a one-dimensional plane flow'
IF (IFX == 1) WRITE (3,*) 'DSMC program DS1 for a cylindrical flow'
IF (IFX == 2) WRITE (3,*) 'DSMC program DS1 for a spherical flow'
WRITE (3,*)
!
WRITE (3,*) 'Interval',NOUT,'Time ',FTIME, ' with',NSAMP,' samples from',TISAMP
WRITE (*,*) 'TOTAL MOLECULES = ',NM
!
NMS=0
DO N=1,NM
  M=IPSP(N)
  NMS(M)=NMS(M)+1
END DO
WRITE (3,*) 'Total simulated molecules =',NM
DO N=1,MSP
  WRITE (*,*) 'SPECIES ',N,' TOTAL = ',NMS(N)
  WRITE (3,*) 'Species ',N,' total = ',NMS(N)
END DO
!
NNN=DINT(TOTMOV)
WRITE (3,*) 'Total molecule moves   =',NNN
NNN=DINT(TOTCOL)
WRITE (3,*) 'Total collision events =',NNN
NNN=DINT(TOTDUP)
WRITE (3,*) 'Total duplicate collisions =',NNN
!
WRITE (3,*) 'Species dependent collision numbers in current sample'
DO N=1,MSP
  WRITE (3,901) TCOL(N,1:MSP)
END DO
901 FORMAT(20G13.5)
CTIME=MCLOCK()
WRITE (3,*) 'Computation time',FLOAT(CTIME)/1000.,'seconds'
WRITE (3,*) 'Collision events per (cpu) second',(TOTCOL-TOTCOLI)*1000.D00/DFLOAT(CTIME)
WRITE (3,*) 'Molecule moves per (cpu) second',(TOTMOV-TOTMOVI)*1000.D00/DFLOAT(CTIME)
WRITE (3,*)
!
IF (IPDF > 0) THEN
  WRITE (3,*) 'Distribution function sampling started at',TISAMP
  WRITE (3,*) 'Total dissociations',NDISSOC
  WRITE (3,*) 'Total recombinations',NRECOMB
  WRITE (3,*) 'Gas temperature',VAR(11,1),' K'
  WRITE (3,*) 'Gas translational temperature',VAR(8,1),' K'
  WRITE (3,*) 'Gas rotational temperature',VAR(9,1),' K'
  WRITE (3,*) 'Gas vibrational temperature',VAR(10,1),' K'
  DO L=1,MSP
    WRITE (3,*) 'Species',L,' overall temperature',VARSP(8,1,L),' K'
    WRITE (3,*) 'Species',L,' translational temperature',VARSP(5,1,L),' K'
    WRITE (3,*) 'Species',L,' rotational temperature',VARSP(6,1,L),' K'
    WRITE (3,*) 'Species',L,' vibrational temperature',VARSP(7,1,L),' K'
  END DO
WRITE (3,*)
END IF
!
IF ((ITYPE(1) == 2).OR.(ITYPE(2) == 2)) WRITE (3,*) 'Surface quantities'
DO JJ=1,2
  IF (ITYPE(JJ) == 2) THEN
    WRITE (3,*)
    WRITE (3,*) 'Surface at',XB(JJ)
    WRITE (3,*) 'Incident sample',VARS(0,JJ)
    WRITE (3,*) 'Number flux',VARS(3,JJ),' /sq m/s'
    WRITE (3,*) 'Inc pressure',VARS(5,JJ),' Refl pressure',VARS(6,JJ)
    WRITE (3,*) 'Pressure', VARS(5,JJ)+VARS(6,JJ),' N/sq m'
    WRITE (3,*) 'Inc y shear',VARS(7,JJ),' Refl y shear',VARS(8,JJ)
    WRITE (3,*) 'Net y shear',VARS(7,JJ)-VARS(8,JJ),' N/sq m'
    WRITE (3,*) 'Net z shear',VARS(9,JJ)-VARS(10,JJ),' N/sq m'
    WRITE (3,*) 'Incident translational heat flux',VARS(11,JJ),' W/sq m'
    IF (MMRM > 0) WRITE (3,*) 'Incident rotational heat flux',VARS(13,JJ),' W/sq m'
    IF (MMVM > 0) WRITE (3,*) 'Incident vibrational heat flux',VARS(15,JJ),' W/sq m'
    WRITE (3,*) 'Total incident heat flux',VARS(29,JJ),' W/sq m'
    WRITE (3,*) 'Reflected translational heat flux',VARS(12,JJ),' W/sq m'
    IF (MMRM > 0) WRITE (3,*) 'Reflected rotational heat flux',VARS(14,JJ),' W/sq m'
    IF (MMVM > 0) WRITE (3,*) 'Reflected vibrational heat flux',VARS(16,JJ),' W/sq m'
    WRITE (3,*) 'Total reflected heat flux',VARS(30,JJ),' W/sq m'
    WRITE (3,*) 'Net heat flux',VARS(31,JJ),' W/sq m'
    WRITE (3,*) 'Slip velocity (y direction)',VARS(19,JJ),' m/s'
    WRITE (3,*) 'Translational temperature slip',VARS(20,JJ),' K'
    IF (MMRM > 0) WRITE (3,*) 'Rotational temperature slip',VARS(21,JJ),' K'
    IF (MSP > 1) THEN
      DO L=1,MSP
        WRITE (3,*) 'Species',L,' percentage',VARS(L+32,JJ)
      END DO
    END IF
    WRITE (3,*)
  END IF
END DO
!
WRITE (3,994) ' Macro or Coll Temps (ITCV,IEAA,IZV):',ITCV,IEAA,IZV
!
PPA=0
DO N=1,NCELLS
  DO M=1,MSP
    PPA(M)=PPA(M)+VARSP(0,N,M)
  END DO
END DO
CLOSE(3)
!
!----------------------------------------------------------------------------
!--write number of reactions
!
IF (MNRE > 0) THEN
  OPEN (3,FILE='DS1REAC.DAT',ACCESS='APPEND')
!
  IF (JCD == 1) THEN
    WRITE (3,*) 'New integrated DSMC chemistry model with no experimentally based rates'
    WRITE (3,*)
    WRITE(3,*) 'GAINS FROM REACTIONS'
    WRITE(3,*) '                          Dissoc.     Recomb. Endo. Exch.  Exo. Exch.'
    DO M=1,MSP
      WRITE (3,*) ' SPECIES',M,TREACG(1,M),TREACG(2,M),TREACG(3,M),TREACG(4,M)
    END DO
    WRITE (3,*)
    WRITE(3,*) 'LOSSES FROM REACTIONS'
    WRITE(3,*) '                          Dissoc.     Recomb. Endo. Exch.  Exo. Exch.'
    DO M=1,MSP
      WRITE (3,*) ' SPECIES',M,TREACL(1,M),TREACL(2,M),TREACL(3,M),TREACL(4,M)
    END DO
    WRITE (3,*)
    WRITE (3,*) 'TOTALS'
    DO M=1,MSP
      WRITE (3,*) ' SPECIES',M,' GAINS',TREACG(1,M)+TREACG(2,M)+TREACG(3,M)+TREACG(4,M),&
                  ' LOSSES',TREACL(1,M)+TREACL(2,M)+TREACL(3,M)+TREACL(4,M)
    END DO
  END IF
!
  IF (JCD == 0) THEN
    AS=0.
    SRR=0.
    DO M=1,MNRE
      AS=AS+REAC(M)
    END DO
    IF (AS == 0.) AS=1.d0
    IF (IREAC == 0) THEN
      WRITE (3,993) FTIME,AS,REAC(1:MNRE)/AS  !no reaction rate sampling
    ELSE
      DO K=1,MNRE
        A=REAC(K)*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP) !# of reaction/V/t
        C=TCOL(LE(K),ME(K))*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP)
        IF(KP(K) >= 0) THEN
          ! nA * nB
          B=(VARSP(1,NSPDF,LE(K))*VARSP(1,NSPDF,ME(K))*VAR(3,NSPDF)**2.d0)/(AVOG*1.d3)                             !exchange/dissociation
        ELSE
          B=(VARSP(1,NSPDF,LE(K))*VARSP(1,NSPDF,ME(K))*VARSP(1,NSPDF,MP(K))*VAR(3,NSPDF)**3.d0)/(AVOG*1.d3)**2.d0  !recombination
        END IF
        SRR(K)=A/B !sampled reaction rate (cm3/mol/s)
        COLRR(K) = C/B
        ! COLRR is the collision rate, but not the collision frequency
        ! For different specie, collision frequency = COLRR * n1 * n2
        ! For the same specie, collision frequency = COLRR * n1 * n2 *0.5


       ! IF (VAR(10,NSPDF) < 6.d3)  SRR(K)=SRR(K)/1.d3        !trick to sample low T rates
      END DO
      WRITE (3,993) FTIME,VAR(8:11,NSPDF),VAR(18,NSPDF),AS,REAC(1:MNRE)/AS,SRR(1:MNRE),COLRR(1:MNRE)
    END IF
  END IF
!
  CLOSE(3)

7787 FORMAT('#Reaction: ',I3,' Sp: ',I2,'+',I2,'->',I2,'+',I2,'+',I2)
7786 FORMAT("#",A13,3X,4(A14,3X))
  DO K=1,MNRE
    IF (KP(K) > 0) THEN !only for dissociation reaction
      WRITE(FILENAME,'("REAC_COUNT_",I3.3,".DAT")') K
      IF (K == 1) FILENAME='REAC_COUNT.DAT'
      OPEN(3,FILE=trim(FILENAME),ACCESS = 'APPEND')
      IF (NOUT < 1) THEN
        WRITE(3,7787) K, LE(K),ME(K),KP(K),LP(K),MP(K)
        WRITE(3,7786) 'Nreac','Evrem (K)','Evrem/D','Evrem (eV)','Time'
      END IF
      IF (REAC(K) > 1.0d0) THEN
        WRITE(3,'(5(E14.6,3x))') REAC(K), EVREM(K)/REAC(K)/BOLTZ, EVREM(K)/REAC(K)/REA(2,K), &
          EVREM(K)/EVOLT/REAC(K),FTIME
      ELSE
        WRITE(3,'(5(E14.6,3x))') 0.d0, 0.d0, 0.d0, 0.d0, FTIME
      END IF
      CLOSE(3)
    END IF
  END DO

END IF
!
!
!---------------------------------------------------------------------------
!--write debug information for imf method
7788 FORMAT('ZONE I=',I6,', T="REAC',I3.3,'", SOLUTIONTIME=',G14.6)
7789 FORMAT(F6.2,3X,21(G14.6,3X))
7790 FORMAT('ZONE I = ',i6,', J = ',i6,', T= "REAC',I3.3,'_',I1,'", DATAPACKING = POINT')
7785 FORMAT(I3,2X,F10.6,2X,7(G14.6,2X))
7784 FORMAT(I3,2X,F10.6,2X,F10.6,2X,14(G14.6,2X))
IF (IREAC == 2 .AND. IMF .ne. 0 .AND. MNRE <= 2.AND. IMFS == 1 .AND. MNRE>0) THEN
  !
  !----- Distribution of Et and Er------------------------
  !
  OPEN (3,FILE='IMF_ETR.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Nomralize by T: ",FTMP0,"K, Current T: ",VAR(8,NSPDF),"K "
  WRITE(3,"(5A)") 'VARIABLES = "E/kT",',&
    '"ET0_N","ET0_P","ET_N","ET_P","ETR_N","ETR_P",',&
    '"ER0_N","ER0_P","ER_N","ER_P","ERR_N","ERR_P",',&
    '"ER01_N","ER01_P","ER1_N","ER1_P","ERR1_N","ERR1_P",',&
    '"PEQ","PCOL","PBOLTZ"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      WRITE(3,7788) 1000,L,FTIME
      K = IMFpair(IREA(1,L),IREA(2,L))

      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE IF (IREA(1,L) > IREA(2,L)) THEN
        IJ = 2; JJ = 1
      END IF

      AA = SUM(NMFET0(:,K))
      BB = SUM(NMFER0(:,IJ,K)); BB2 = SUM(NMFER0(:,JJ,K))

      A = SUM(NMFET(:,K))
      C = SUM(NMFER(:,IJ,K));   C2 = SUM(NMFER(:,JJ,K))

      AS = SUM(NMFETR(:,L))
      AU = SUM(NMFERR(:,1,L)); AT = SUM(NMFERR(:,2,L))


      CC = 2.5d0-SPM(3,IREA(1,L),IREA(2,L))
      DO N = 1,1000
        B = DFLOAT(N)*0.01d0
        D = 0.5d0*SPI*(ERF(dsqrt(B))-ERF(dsqrt(B-0.01d0)))
        !WRITE(3,7789) B,NMFET0(N),NMFET0(N)/AA,&
        !NMFER0(N), NMFER0(N)/BB,&
        !NMFET(L,N),DFLOAT(NMFET(L,N))/A,&
        !NMFER(L,N),DFLOAT(NMFER(L,N))/C, &
        !NMFETR(L,N),DFLOAT(NMFETR(L,N))/AS, &
        !NMFERR(L,N),DFLOAT(NMFERR(L,N))/SRR,&
        !(0.99d0+B)*dexp(-B+0.01d0)-(1.0d0+B)*dexp(-B),&
        !dcosh(B-0.01d0)-dcosh(B)-dsinh(B-0.01d0)+dsinh(B)
        WRITE(3,7789) B,&
          NMFET0(N,K),   NMFET0(N,K)/AA,  NMFET(N,K),    NMFET(N,K)/A,&
          NMFETR(N,L),   NMFETR(N,L)/AS,  NMFER0(N,IJ,K),NMFER0(N,IJ,K)/BB, &
          NMFER(N,IJ,K), NMFER(N,IJ,K)/C, NMFERR(N,1,L), NMFERR(N,1,L)/AU, &
          NMFER0(N,JJ,K),NMFER0(N,JJ,K)/BB2, &
          NMFER(N,JJ,K), NMFER(N,JJ,K)/C2,NMFERR(N,2,L), NMFERR(N,2,L)/AT, &
          gamain(B,1.5d0,NN) - gamain(B-0.01d0,1.5d0,NN),&
          gamain(B,CC,NN) - gamain(B-0.01d0,CC,NN),&
          dcosh(B-0.01d0)-dcosh(B)-dsinh(B-0.01d0)+dsinh(B)
      END DO
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
  !
  ! -- Distribution of Ev
  !
  OPEN (3,FILE='IMF_EV.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Initial T: ",FVTMP0,"K Current T: ",VAR(10,NSPDF),"K"
  WRITE(3,"(3A)")'VARIABLES = "v","Ev (eV)","Ev2 (eV)",', &
    '"Ev0_N","Ev0_P","Ev_N","Ev_P","EvR_N","EvR_P","P0",',&
    '"Ev01_N","Ev01_P","Ev1_N","Ev1_P","EvR1_N","EvR1_P","P1"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      K = IMFpair(IREA(1,L),IREA(2,L))
      ! write chemical reaction
      WRITE(3,*)
      WRITE(3,'("#",I3,"+",I3,"->",I3,"+",I3,"+",I3)') LE(L),ME(L),KP(L),LP(L),MP(L)

      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE
        IJ = 2; JJ = 1
      END IF

      BB = SUM(NMFEV0(:,IJ,K)); AA = SUM(NMFEV(:,IJ,K)); C = SUM(NMFEVR(:,1,L))
      IF (BB .ge. 1.0d0) THEN
        BB = 1.0d0/BB
      ELSE
        BB = 0.0d0
      END IF
      IF (AA .ge. 1.0d0) THEN
        AA = 1.0d0/AA
      ELSE
        AA = 0.0d0
      END IF
      IF (C  .ge. 1.0d0) THEN
        C = 1.0d0/C
      ELSE
        C = 0.0d0
      END IF

      KK = 100
      IF (ISPV(IREA(2,L)) == 0) THEN
        ! this is an atom-diatom dissociation
        ! KK is the maximum vibrational level
        IF (IVMODEL(IREA(1,L),1) == 1) KK = IVMODEL(IREA(1,L),2)
        !--- write zone header
        WRITE(3,7788) KK+1,L,FTIME  !write zone header
        WRITE(3,"(A)") 'PASSIVEVARLIST = [3,11-17]'
        !--- write data
        DO N = 0,KK
          CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
          WRITE(3,7785) N,EVIB/EVOLT, NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB,&
            NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,NMFEVR(N,1,L),NMFEVR(N,1,L)*C,&
            DEXP(-EVIB/BOLTZ/VAR(10,NSPDF))
        END DO
      ELSE
        BB2 = SUM(NMFEV0(:,JJ,K));
        IF (BB2 .ge. 1.0d0) THEN
          BB2 = 1.0d0/BB2
        ELSE
          BB2 = 0.0d0
        END IF
        A = SUM(NMFEV(:,JJ,K))
        IF (A .ge. 1.0d0) THEN
          A = 1.0d0/A
        ELSE
          A = 0.0d0
        END IF
        CC = SUM(NMFEVR(:,2,L))
        IF (CC .ge. 1.0d0) THEN
          CC = 1.0d0/CC
        ELSE
          CC = 0.0d0
        END IF

        KR = 100; JR = 100;
        IF (IVMODEL(IREA(1,L),1) == 1)    KR = IVMODEL(IREA(1,L),2)
        IF (IVMODEL(IREA(2,L),1) == 1)    JR = IVMODEL(IREA(2,L),2)
        KK = MAX(KR,JR)

        !--- write zone header
        WRITE(3,7788) KK+1,L,FTIME
        !--- write data
        DO N=0,MIN(KR,JR)
          CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
          CALL VIB_ENERGY(AS,N,1,IREA(2,L))
          WRITE(3,7784) N,EVIB/EVOLT,AS/EVOLT,&
            NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB, NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,&
            NMFEVR(N,1,L),NMFEVR(N,1,L)*C, DEXP(-EVIB/BOLTZ/VAR(10,NSPDF)),&
            NMFEV0(N,JJ,K),NMFEV0(N,JJ,K)*BB2, NMFEV(N,JJ,K),NMFEV(N,JJ,K)*A,&
            NMFEVR(N,2,L),NMFEVR(N,2,L)*CC, DEXP(-AS/BOLTZ/VAR(10,NSPDF))
        END DO

        IF (KR < JR) THEN
          DO N=KR+1,JR
            CALL VIB_ENERGY(AS,N,1,IREA(2,L))
            WRITE(3,7784) N,0.0d0,AS/EVOLT,&
              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0 ,&
              NMFEV0(N,JJ,K),NMFEV0(N,JJ,K)*BB2, NMFEV(N,JJ,K),NMFEV(N,JJ,K)*A,&
              NMFEVR(N,2,L),NMFEVR(N,2,L)*CC, DEXP(-AS/BOLTZ/VAR(10,NSPDF))
          END DO
        ELSE IF (KR > JR) THEN
          DO N=JR+1,KR
            CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
            WRITE(3,7784) N,EVIB/EVOLT,0.0d0,&
              NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB, NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,&
              NMFEVR(N,1,L),NMFEVR(N,1,L)*C, DEXP(-EVIB/BOLTZ/VAR(10,NSPDF)),&
              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0
          END DO
        END IF
      END IF
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
  !
  ! -- Reaction probability
  OPEN (3,FILE='IMF_Prob.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Nomralize by T: ",FTMP0,"K, Current T:",VAR(8,NSPDF),"K"
  WRITE(3,"(A)")'VARIABLES = "Et/kT", "Ev (eV)","N","NR1","NR2","P1","P2"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      K = IMFpair(IREA(1,L),IREA(2,L))
      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE
        IJ = 2; JJ = 1
      END IF

      ! write data for the dissociating particle
      KK = 100
      IF (IVMODEL(IREA(1,L),1) == 1)     KK = IVMODEL(IREA(1,L),2)
      WRITE(3,7790) 1000, KK+1,L,1

      DO N = 0,KK
        DO M = 1,1000
          CALL VIB_ENERGY(C,N,1,IREA(1,L))
          A = 0.0D0;       B = 0.0D0
          IF (NMFVT0(N,M,IJ,K) .NE. 0.0d0)  A = NMFVT(N,M,IJ,K)/NMFVT0(N,M,IJ,K)
          IF (NMFVT(N,M,IJ,K) .NE. 0.0d0)   B = NMFVTR(N,M,1,L) / NMFVT(N,M,IJ,K)

          WRITE(3,7789) DFLOAT(M)*0.01d0, C/EVOLT,&
            NMFVT0(N,M,IJ,K),NMFVT(N,M,IJ,K),NMFVTR(N,M,1,L),A,B
        END DO
      END DO
      WRITE(3,*)

      ! write data for the collider
      KK = 100
      IF (IVMODEL(IREA(2,L),1) == 1)  KK = IVMODEL(IREA(2,L),2)
      WRITE(3,7790) 1000, KK+1,L,2

      DO N = 0,KK
        DO M = 1,1000
          CALL VIB_ENERGY(C,N,1,IREA(2,L))
          A = 0.0D0;       B = 0.0D0
          IF (NMFVT0(N,M,JJ,K) .NE. 0.0d0)  A = NMFVT(N,M,JJ,K)/NMFVT0(N,M,JJ,K)
          IF (NMFVT(N,M,JJ,K) .NE. 0.0d0)   B = NMFVTR(N,M,2,L) / NMFVT(N,M,JJ,K)

          WRITE(3,7789) DFLOAT(M)*0.01d0, C/EVOLT,&
            NMFVT0(N,M,JJ,K),NMFVT(N,M,JJ,K),NMFVTR(N,M,2,L),A,B
        END DO
      END DO
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
END IF


!
!----------------------------------------------------------------------------
!--write flowfield overall properties
!
IF (ISF >= 1) OPEN (3,FILE='DS1FP00_'//E//'.DAT')
IF (ISF == 0) OPEN (3,FILE='DS1FP00.DAT')

!
WRITE(TNAME,7791) NOUT
7791 FORMAT('ZONE T = "'i4.4,'"')
!
WRITE (3,994) '# Flowfield Overall Properties at FTIME:',FTIME
WRITE (3,994) '# NCELLS:', NCELLS
WRITE (3,994) '# NSAMP: ', NSAMP
WRITE (3,994) '# X-coord.     Cell   Sample     Number Dens.   Density   u velocity &
            &  v velocity   w velocity   Trans. Temp.   Rot. Temp.   Vib. Temp. &
            &  Temperature  Mach no.     Mols/cell    m.c.t        m.f.p     &
            &  mcs/mfp        speed      Pressure        TTX          TTY          TTZ     &
            &  Tau_xy      Qx &
            &  dtm/mct     <dx/mfp>      Fx          Fy            Fz         Qtransfer &
            &  Species Fractions'
WRITE (3,994,advance='no') 'VARIABLES = "x(m)","ID","NSample","n","rho","u","v","w",&
            & "Tt","Tr","Tv","T","Ma","#Mol/Cell","mct","mfp","mcs/mfp","Speed","p=nkT",&
            & "Ttx","Tty","Ttz","Tau_xy","Qx","dtm/mct","dx/mfp","Fx","Fy","Fz","Qtransfer",'
DO N=1,MSP-1
  WRITE(3,"(A3,I2,A2)",advance='no') '"Sp',N,'",'
END DO
WRITE(3,"(A3,I2,A1)") '"Sp',MSP,'"'

WRITE (3,994) TRIM(TNAME)
WRITE (3,994) 'SOLUTIONTIME =',FTIME
WRITE (3,992) 'STRANDID = ',QCTMODEL
!
DO N=1,NCELLS
  WRITE (3,996) VAR(1,N),N,VAR(2:23,N),DTM/VAR(14,N),(CELL(3,N)-CELL(2,N))/DFLOAT(NCIS)/VAR(15,N),&
                CST(1:4,N)/CST(0,N),VARSP(1,N,1:MSP)
END DO
CLOSE(3)
!
!----------------------------------------------------------------------------
!--write flowfield properties per species
!
DO L=1,MSP
    WRITE(FILENAME,775) L
775 FORMAT('DS1FP',i2.2)
    OPEN (3,FILE=TRIM(FILENAME)//'.DAT')
!
  WRITE (3,994) '# Flowfield Properties Species:',L
  WRITE (3,994) '# NCELLS:', NCELLS
  WRITE (3,994) '# NSAMP: ', NSAMP
  WRITE (3,994) '# X-coord.     Cell   Sample    Fraction     Species TTx   Species TTy  Species TTz &
             & Trans. Temp.  Rot. Temp.   Vib. Temp.   Spec. Temp  u Diff Vel   v Diff Vel   w Diff Vel. Q_rot_x  Q_vib_x'
  DO N=1,NCELLS
    WRITE (3,996) VAR(1,N),N,VARSP(0,N,L),VARSP(1,N,L),VARSP(2:13,N,L)
  END DO
CLOSE(3)
END DO
!
!----------------------------------------------------------------------------
!--write composition of a reacting gas as a function of time
!
OPEN (10,FILE='COMPOSITION.DAT',ACCESS='APPEND')
A=0.d0
AS=NM
IF (IENERS > 0) THEN
  IF (NOUT <= 1) ENERS0=ENERS
  A=ENERS/ENERS0-1.d0
  ENERS0=ENERS
END IF
WRITE (10,993) FTIME,VAR(8:11,NSPDF),VAR(18,NSPDF),NMS(1:MSP)/AS,A,NM,VAR(5,1),VAR(5,NCELLS)
CLOSE (10)
!
!----------------------------------------------------------------------------
!--write pdf files
!
IF(IPDF > 0) THEN
!
!--compute normalized velocity and speed pdfs
!
  PDF=0.d0
  PDFS=0.d0
  DO N=1,NBINS
    DO L=1,MSP
      DO K=1,3
        PDF(N,K)=PDF(N,K)+BINS(N,K,L)/DBINV(3)/BIN(0,K)
        PDFS(N,K,L)=BINS(N,K,L)/DBINV(3)/BINS(0,K,L)
      END DO
      PDF(N,4)=PDF(N,4)+BINS(N,4,L)/DBINC(3)/BIN(0,4)
      PDF(N,5)=PDF(N,5)+BINS(N,5,L)/DBINE(3)/BIN(0,5)
      PDFS(N,4,L)=BINS(N,4,L)/DBINC(3)/BINS(0,4,L)
      PDFS(N,5,L)=BINS(N,5,L)/DBINE(3)/BINS(0,5,L)
    END DO
  END DO
!
!--velocity components
  OPEN (7,FILE='DS1DVEL.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Interval, Equilibrium pdf, U V W pdfs total, U V W pdfs spec1, U V W pdfs spec2 ... '
  DO N=1,NBINS
    A=DBINV(1)+DBINV(3)*(DFLOAT(N)-0.5d0)                !normalized U/VMP
    B=(1.d0/SPI)/DEXP(A*A)                               !Boltzmann distribution
    WRITE (7,993) A,B,PDF(N,1:3),(PDFS(N,1:3,L),L=1,MSP) !pdfs
  END DO
  CLOSE(7)
!
!--speed
  OPEN (7,FILE='DS1DSPD.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Interval, Equilibrium pdf, C pdf total, C pdf spec1, C pdf spec2 ... '
  DO N=1,NBINS
    A=DBINC(1)+DBINC(3)*(DFLOAT(N)-0.5d0)                !normalized C/VMP
    B=((4.d0/SPI)*A*A)/DEXP(A*A)                         !Boltzmann distribution; ftr- Eqn 15.25 or N.5 from Laurendeau
    WRITE (7,993) A,B,PDF(N,4),(PDFS(N,4,L),L=1,MSP)     !pdfs
  END DO
  CLOSE(7)
!
!--thermal (translational) energy
  OPEN (7,FILE='DS1DTEN.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Energy (eV), U V W pdfs total, U V W pdfs spec1, U V W pdfs spec2 ... '
  DO N=1,NBINS
    A=DBINE(1)+DBINE(3)*(DFLOAT(N)-0.5d0)                !energy in eV
    WRITE (7,993) A,PDF(N,5),(PDFS(N,5,L),L=1,MSP) !pdfs
  END DO
  CLOSE(7)
!
!--rotational energy
  IF (MMRM > 0) THEN
    DO L=1,MSP
      IF (ISPR(1,L) > 0) THEN
        EQROT=1.d0
        DO M=1,100
          A=(DFLOAT(M)-0.5D00)*0.1D00   !--rotational values
          BB=.5*ISPR(1,L)
          !Boltzmann distribution; frot - Eqn 11.21 from Bird94
          EQROT(L,M)=(1.d0/GAM(BB))*(A**(BB-1.d0))*DEXP(-A)*0.1D00
        END DO
        WRITE(FILENAME,777) L
        777 FORMAT('DS1DROT',i2.2)
        OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
        WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
        WRITE (7,994) '# Equilibrium based on the rotational temperature'
        WRITE (7,994) '# Species and number of modes:',L,ISPR(1,L)
        WRITE (7,994) '# erot/kT       sample      f_rot    f_rot_equilib    ratio'
        DO M=1,100
          A=(DFLOAT(M)-0.5D00)*0.1D00
          B=DFLOAT(NDROT(L,M))/BINS(0,1,L) !NDROT and EQROT are fractions in interval 0.1, so multiply by 10 to obtain pdfs
          WRITE (7,993) A,DFLOAT(NDROT(L,M)),B*10.D00,EQROT(L,M)*10.D00,B/EQROT(L,M)
        END DO
        CLOSE (7)
      END IF
    END DO
  END IF
!
!--vibrational levels
  IF (MMVM > 0) THEN
    DO I=1,NSCELLS
      J=NSVEC(I) !sampling cell
      F =0.0d0
      DO L=1,MSP
        IF (ISPV(L) > 0 .AND. (L==1 .or. L==3 .or. L==5) ) THEN !plotting only O2 or N2 populations
          K=1 !single vibrational mode
!
          !calculate the vibrational partition function
          EQVIBVT=1.d0
          EQVIBOT=1.d0
          QVIBVT=0.d0
          QVIBOT=0.d0
          MAXLEVEL=IVMODEL(L,2)
          DO M=0,MAXLEVEL
            CALL VIB_ENERGY(EVIB,M,K,L)
            QVIBVT=QVIBVT+DEXP(-EVIB/(BOLTZ*VARSP(7,J,L)))       !based on Tv_mode
            QVIBOT=QVIBOT+DEXP(-EVIB/(BOLTZ*VARSP(8,J,L))) !based on Tov
          END DO
!
          !calculate Boltzmann distribution function
          DO M=0,MAXLEVEL
            CALL VIB_ENERGY(EVIB,M,K,L)
            EQVIBVT(K,L,M)=DEXP(-EVIB/(BOLTZ*VARSP(7,J,L)))/QVIBVT       !based on Tv_mode
            EQVIBOT(K,L,M)=DEXP(-EVIB/(BOLTZ*VARSP(8,J,L)))/QVIBOT !based on Tov
          END DO
!
          ! TAG: SHOCK SAMPLE
          ! Use the following naming of file for vibrational level sampling
          ! WRITE(FILENAME,778) L,J
          ! 778 FORMAT('DS1DVIB',i2.2,'_',i4.4)

          ! Regular naming
          WRITE(FILENAME,879) L
          879 FORMAT('DS1DVIB',i2.2)

          WRITE(TNAME,779) L,J !NOUT
          779 FORMAT('ZONE T = "',i2.2,'_',i4.4,'"')
          OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
          WRITE (7,994) '# NS, X-coord, Tt, Tr, Tv, To, t, MCT: ', J,VAR(1,J),VAR(8:11,J),FTIME,MCT
          WRITE (7,994) '# Equilibrium based on both vibrational temperature and overall gas temperature'
          WRITE (7,994) '# Species and number of modes:',L,ISPV(L)
          WRITE (7,994) '# level      Ev       sample      f_vib       f_equil(Tvib) ratio Eq.   f_equil(Tov)  ratio'
          WRITE (7,994) TRIM(TNAME)
          WRITE (7,994) 'SOLUTIONTIME =',FTIME
          WRITE (7,992) 'STRANDID = ',QCTMODEL
          !DO M=0,MAXLEVEL
          DO M=0,IVMODEL(L,2)
            F(1:ISPV(L),M,L)=DFLOAT(NDVIB(J,1:ISPV(L),L,M))/DFLOAT(NDVIB(J,0,L,0))
            CALL VIB_ENERGY(EVIB,M,K,L)
            A=EVIB/EVOLT !energy in eV
            WRITE (7,993) M,A,DFLOAT(NDVIB(J,1:ISPV(L),L,M)),F(1:ISPV(L),M,L),&
                          EQVIBVT(1:ISPV(L),L,M),F(1:ISPV(L),M,L)/EQVIBVT(1:ISPV(L),L,M),&
                          EQVIBOT(1:ISPV(L),L,M),F(1:ISPV(L),M,L)/EQVIBOT(1:ISPV(L),L,M)
          END DO
          CLOSE (7)
!
        END IF
      END DO
    END DO
  END IF
!
END IF

!-- vibrational state-specific rates
IF (IREAC == 2 .AND. IMF .ne. 0 .AND. MNRE <= 2.AND. IMFS == 1 .AND. MNRE>0 .AND. NSCELLS==1) THEN
  OPEN(3, FILE='IMF_Vrate.DAT')
  WRITE(3,"(A,G14.6,A,F10.3)") '# time:',FTIME, ' VT: ',VAR(10,NSPDF)
  WRITE(3,"(A)") 'VARIABLES = "Evib(eV)","v","Nreac","Rate (cm3/mol/s)"'
  DO L = 1,MNRE
    NN = LE(L)
    IF (KP(L) > 0 .and. IVMODEL(NN,1) == 1) THEN ! dissociation alone
      WRITE(3,*)
      WRITE(3,'("#",I3,"+",I3,"->",I3,"+",I3,"+",I3)') LE(L),ME(L),KP(L),LP(L),MP(L)
      WRITE(3,7788) IVMODEL(NN,2)+1,L,FTIME  !write zone header
      DO N=0,IVMODEL(NN,2)
        CALL VIB_ENERGY(EVIB,N,1,NN)
        A=NMFEVR(N,1,L)*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP) !# of reaction/V/t
        B=(VARSP(1,NSPDF,LE(L))*F(1,N,NN)*VARSP(1,NSPDF,ME(L))*VAR(3,NSPDF)**2.d0)/(AVOG*1.d3)    !number density

        IF (NMFEVR(N,1,L) == 0) THEN
          WRITE(3,'(F12.6,2X,I3,2X,E14.6,2X,E14.6)') EVIB/EVOLT,N,NMFEVR(N,1,L),0.0d0
        ELSE
          WRITE(3,'(F12.6,2X,I3,2X,E14.6,2X,E14.6)') EVIB/EVOLT,N,NMFEVR(N,1,L),A/B
        END IF
        ! IF (VAR(10,NSPDF) < 6.d3)  SRR(K)=SRR(K)/1.d3        !trick to sample low T rates
      END DO
      write(3,*)
    END IF
  END DO
  CLOSE(3)
END IF
!
!--sampled pre- and post-reaction vibrational levels
IF (MNRE > 0) THEN
!
  J=1 !pre-reaction
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
    B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
    C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
    IF(A == 0) A=1.d0
    IF(B == 0) B=1.d0
    IF(C == 0) C=1.d0
    WRITE(FILENAME,780) KK
    780 FORMAT('DS1VIBL1_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled pre-reaction vibrational levels'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,B,C,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# I, Molec1 Sample and Fraction, Molec2 Sample and Fraction, , Molec3 Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NPVIB(J,KK,1,1:3,I),NPVIB(J,KK,1,1:3,I)/A,&
                      NPVIB(J,KK,2,1:3,I),NPVIB(J,KK,2,1:3,I)/B,&
                      NPVIB(J,KK,3,1:3,I),NPVIB(J,KK,3,1:3,I)/C
    END DO
    CLOSE (7)
  END DO
!
  J=2 !post-reaction
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
    B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
    C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
    IF(A == 0) A=1.d0
    IF(B == 0) B=1.d0
    IF(C == 0) C=1.d0
    WRITE(FILENAME,782) KK
    782 FORMAT('DS1VIBL2_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled post-reaction vibrational levels'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,B,C,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# I, Molec1 Sample and Fraction, Molec2 Sample and Fraction, , Molec3 Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NPVIB(J,KK,1,1:3,I),NPVIB(J,KK,1,1:3,I)/A,&
                      NPVIB(J,KK,2,1:3,I),NPVIB(J,KK,2,1:3,I)/B,&
                      NPVIB(J,KK,3,1:3,I),NPVIB(J,KK,3,1:3,I)/C
    END DO
    CLOSE (7)
  END DO
!
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NEVIB(KK,1,:))
    IF(A == 0) A=1.d0
    WRITE(FILENAME,783) KK
    783 FORMAT('DS1VIBE1_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled pre-reaction total vibrational energiy'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# Temperatures:',VAR(8:11,NSPDF)
    WRITE (7,994) '# [Ev/(BOLTZ*1e3)], Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NEVIB(KK,1,I),NEVIB(KK,1,I)/A  !pre-reaction total vibrational energy
    END DO
    CLOSE (7)
  END DO
!
  IF (IREAC == 2) THEN
    FEVIB=0.d0
    DO KK=1,MNRE
      A=SUM(NEVIB(KK,1,:))
      IF(A == 0) A=1.d0
      DO I=0,100
        FEVIB(KK,1,1,I)=NEVIB(KK,1,I)/A  !pre-reaction total vibrational energies
      END DO
    END DO
!
    J=1
    FPVIB=0.d0
    DO KK=1,MNRE
      A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
      B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
      C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
      IF(A == 0) A=1.d0
      IF(B == 0) B=1.d0
      IF(C == 0) C=1.d0
      DO I=0,100
        FPVIB(KK,1,1,1:3,I)=NPVIB(J,KK,1,1:3,I)/A  !pre-reaction vibrational levels pdfs
        FPVIB(KK,1,2,1:3,I)=NPVIB(J,KK,2,1:3,I)/B
        FPVIB(KK,1,3,1:3,I)=NPVIB(J,KK,3,1:3,I)/C
      END DO
      ZCHECK=1234567
      WRITE(FILENAME,786) KK,J
      786 FORMAT('DS1VIBF1_',i2.2,'_',i2.2)
      787 CONTINUE
      OPEN (7,FILE=TRIM(FILENAME)//'.BIN',FORM='UNFORMATTED',ERR=787)
      WRITE (7) VAR(11,NSPDF),FEVIB(KK,J,:,:),FPVIB(KK,J,:,:,:),ZCHECK
      CLOSE (7)
    END DO
  END IF
END IF
!
!----------------------------------------------------------------------------
!--sample molecule properties for SNAPSHOT files
!
NSS=0
DO K=1,NSNAP
  CALL ZGF(RANF,IDT)
  N=INT(RANF*DFLOAT(NM))+1
  NSS=NSS+1
  ISNAP(1,NSS)=IPSP(N)
  PSNAP(1,NSS)=PX(N)
  VSNAP(:,NSS)=PV(:,N)
END DO
!
!----------------------------------------------------------------------------
!--write binary files for post processing of multiple run cases and snapshot

ZCHECK=1234567
!
102 CONTINUE
  !OPEN (7,FILE='DS1OUT_'//E//'.BIN',FORM='UNFORMATTED',ERR=102)
  OPEN (7,FILE='DS1OUT.BIN',FORM='UNFORMATTED',ERR=102)
  WRITE (7)IFX,NOUT,FTIME,NSAMP,TISAMP,NM,TOTMOV,TOTCOL,PCOLLS,TOTDUP,TCOL,CSCR,&
          CTIME,TOTCOLI,TOTMOVI,NDISSOC,NRECOMB,&
          ITYPE,XB,VARS,NCELLS,VARSP,&  !end of general properties
          JCD,TREACG,TREACL,ITCV,IEAA,IZV,REAC,VAR,DTM,CELL,NCIS,CST,UVFX,NMS,& !end of reactions,flow properties, and composition
          IPDF,IPRS,BIN,BINS,DBINV,SPI,PDF,PDFS,NSPDF,NDROT,NDVIB,DBINC,DBINE,ISNAP,PSNAP,VSNAP,SP,ZCHECK    !end of sample pdf and snapshot
  CLOSE(7)

WRITE (*,*) 'Output and binary files are written'
!
!----------------------------------------------------------------------------
!--I/O format
!
999 FORMAT (I5,50G13.5)
998 FORMAT (G270.0)
997 FORMAT (G175.0)
996 FORMAT (G13.5,I5,50G13.5)
995 FORMAT (22G13.5)
994 FORMAT (A,22G13.5)
993 FORMAT (99G13.5)
992 FORMAT (A,I2.2)
!
!----------------------------------------------------------------------------
!--reset collision and transit times etc.
!
DTMI=DTM
DTM=DTM*2.
!--this makes it possible for DTM to increase, it will be reduced as necessary
DO N=1,NCCELLS
!
  NN=ICCELL(3,N)
  B=(CELL(3,NN)-CELL(2,NN))/DFLOAT(NCIS)     !--collision cell width
!
  IF ((NN > 0).AND.(NN <= NCELLS)) THEN
!
    IF (VAR(13,NN) > 20.D00) THEN
!--consider the local collision rate
      CCELL(3,N)=0.5D00*VAR(14,NN)*CPDTM
!--look also at collision cell transit time based on the local flow speed
      A=0.5D00*(B/(ABS(VAR(5,NN))))*TPDTM
      IF (A < CCELL(3,N)) CCELL(3,N)=A
      IF (2.D00*CCELL(3,N) < DTM) DTM=2.D00*CCELL(3,N)
    ELSE
!-- base the time step on a collision cell transit time at the refence vmp
      A=TPDTM*B/VMPM
      IF (A < CCELL(3,N)) CCELL(3,N)=A
    END IF
    IF (1.9D00*CCELL(3,N) < DTM) DTM=1.9D00*CCELL(3,N)
!
  END IF
END DO
!
WRITE (9,*) 'DTM changes  from',DTMI,' to',DTM
WRITE (*,*) 'DTM changes  from',DTMI,' to',DTM
!
!----------------------------------------------------------------------------
!--update output interval
!
TPOUT=OUTRAT
IF (ISF > 0) THEN
  IF ((NOUT >= 0  ).AND.(NOUT < 100)) TPOUT=OUTRAT*.1d0
  IF ((NOUT >= 0  ).AND.(NOUT < 100)) TPOUT=OUTRAT*.1d0
  IF ((NOUT >= 100).AND.(NOUT < 150)) TPOUT=OUTRAT*.2d0
  IF ((NOUT >= 150).AND.(NOUT < 170)) TPOUT=OUTRAT*.5d0
  ! IF ((NOUT >= 200)) TPOUT=OUTRAT*INT(.9999999+(NOUT-190)/10) !comment for RELAX.DAT
!
  IF (IRELAX .eq. 0) THEN
  ! for special case, accelerate convergence
    IF (NOUT >= 50 ) CPDTM=0.005
    IF (NOUT >= 100) CPDTM=0.010
    IF (NOUT >= 150) CPDTM=0.050
    IF (NOUT >= 200) CPDTM=0.100
  END IF
!  IF ((CPDTM < 0.02  ).AND.(NOUT >= 200)) CPDTM=0.02
!  IF ((CPDTM < 0.2   ).AND.(NOUT >= 250)) CPDTM=0.2
END IF
!
! IF (ISF==0) THEN !special coding for shockwave sampling
 ! TPOUT=OUTRAT*2.d0
! END IF
!
IF (MNRE > 0) THEN
  IF ((IRM == 201).AND.(FTIME > 4.d-5)) ISF=0 !steady state is achieved
  IF ((IRM == 202).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 206).AND.(FTIME > 4.d-5)) ISF=0
  IF ((IRM == 213).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 214).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 218).AND.(FTIME > 4.d-5)) ISF=0
  IF ((IRM == 227).AND.(FTIME > 2.d-6)) ISF=0
  IF ((IRM == 229).AND.(FTIME > 2.d-6)) ISF=0
  IF ((IRM == 231).AND.(FTIME > 1.d-1)) ISF=0
END IF
IF ((IGS == 3).AND.(FTIME > 3.d-3)) ISF=0
!
DTSAMP=DTSAMP*DTM/DTMI
DTOUT=DTSAMP*TPOUT
TOUT=FTIME
!
WRITE (9,*) 'NOUT:',NOUT,' OUTRAT:',TPOUT
WRITE (*,*) 'NOUT:',NOUT,' OUTRAT:',TPOUT
!
!----------------------------------------------------------------------------
!
RETURN
!
END SUBROUTINE OUTPUT_RESULTS
