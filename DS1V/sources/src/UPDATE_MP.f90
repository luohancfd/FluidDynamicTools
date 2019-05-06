!
!*****************************************************************************
!
SUBROUTINE UPDATE_MP
!
!--author: isebasti
!--update macroscopic properties
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: NS,K,KK,L,N,NMCR,NBC,KCELLS(2,10)     !need to set NBC and KCELLS in some commom block
REAL(KIND=8) :: A,B,C,SMCR,DOF,AVW,UU,SVDF,VDOFM,TVIBM,EVIBM,DSUM(0:12),AA,BB,SN,UBMEAN(2)
REAL(KIND=8), DIMENSION(MSP) :: TVIB,VDOF,SDOF
REAL(KIND=8), DIMENSION(MMVM,MSP) :: TV,THCOL
REAL(KIND=8), DIMENSION(NCELLS,MMVM,MSP) :: DF
!
!--set cells close to boundaries and some inner cells
!
!IF (NBC < NCELLS) NBC=NCELLS
!DO N=1,NBC
!  KCELLS(1,N)=N
!  KCELLS(2,N)=NCELLS-(N-1)
!END DO
!
!----------------------------------------------------------------------------
!
VAR=0.D00
VARSP=0.
SMCR=0
NMCR=0
VDOFM=0.
!
!--can use the lines below to consider only cells close to boundaries
!NBC=1  !consider a maximum of NBC sampling cells from each boundary
!DO NS=1,2
! DO KK=1,NBC
!  N=KCELLS(NS,KK)
!
!--consider all NCELLS sampling cells
DO N=1,NCELLS  !loop is the same as in OUTPUT subroutine
!
  A=FNUM/(CELL(4,N)*NSAMP)
  IF (IVB == 1) A=A*((XB(2)-XB(1))/(XB(2)+VELOB*0.5D00*(FTIME+TISAMP)-XB(1)))**(IFX+1)
!--check the above for non-zero XB(1)
  DSUM=0.
  NMCR=NMCR+1
  DO L=1,MSP
    DSUM(0)=DSUM(0)+CS(0,N,L)
    DSUM(1)=DSUM(1)+CS(1,N,L)
    DSUM(2)=DSUM(2)+SP(5,L)*CS(1,N,L)
    DO K=1,3
      DSUM(K+2)=DSUM(K+2)+SP(5,L)*CS(K+1,N,L)
      IF (CS(1,N,L) > 0.1D00) THEN
        VARSP(K+1,N,L)=CS(K+4,N,L)/CS(1,N,L)  !--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
        VARSP(K+8,N,L)=CS(K+1,N,L)/CS(1,N,L)  !--VARSP(9,10,11 are temporarily the mean of the velocities
      END IF
    END DO
    DSUM(6)=DSUM(6)+SP(5,L)*(CS(5,N,L)+CS(6,N,L)+CS(7,N,L))
    DSUM(10)=DSUM(10)+SP(5,L)*CS(5,N,L)
    DSUM(11)=DSUM(11)+SP(5,L)*CS(6,N,L)
    DSUM(12)=DSUM(12)+SP(5,L)*CS(7,N,L)
    IF (CS(1,N,L) > 0.5D00) THEN
      DSUM(7)=DSUM(7)+CS(5,N,L)+CS(6,N,L)+CS(7,N,L)
    END IF
    IF (ISPR(1,L) > 0) THEN
      DSUM(8)=DSUM(8)+CS(8,N,L)
      DSUM(9)=DSUM(9)+CS(1,N,L)*ISPR(1,L)
    END IF
  END DO
  AVW=0.
  DO L=1,MSP
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
  END DO
!
  IF (IVB == 0) VAR(1,N)=CELL(1,N)
  IF (IVB == 1) THEN
    C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV)   !--new DDIV
    VAR(1,N)=XB(1)+(DFLOAT(N-1)+0.5)*C
  END IF
  VAR(2,N)=DSUM(0)
  IF (DSUM(1) > 0.5) THEN
    VAR(3,N)=DSUM(1)*A               !--number density Eqn. (4.28)
    VAR(4,N)=VAR(3,N)*DSUM(2)/DSUM(1) !--density Eqn. (4.29)
    VAR(5,N)=DSUM(3)/DSUM(2)          !--u velocity component Eqn. (4.30)
    VAR(6,N)=DSUM(4)/DSUM(2)          !--v velocity component
    VAR(7,N)=DSUM(5)/DSUM(2)          !--w velocity component
    UU= VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2
    IF (DSUM(1) > 1) THEN
      VAR(8,N)=(ABS(DSUM(6)-DSUM(2)*UU))/(3.D00*BOLTZ*DSUM(1))  !--translational temperature Eqn. (4.39)
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
      DO L=1,MSP
        VDOF(L)=0.
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (CS(K+8,N,L) > 0.) THEN
              A=CS(K+8,N,L)/CS(1,N,L)
              TV(K,L)=SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/A)  !--Eqn.(4.45) - assuming SHO
              DF(N,K,L)=2.d0*A/(BOLTZ*TV(K,L)) !--Eqn. (11.28) Bird94 - general definition
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
              TVIB(L)=FVTMP(1)
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
    VAR(13,N)=DSUM(0)/NSAMP/DFLOAT(NCIS)
    IF (COLLS(N) > 2.) THEN
!--mean collision time
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
! END DO  !loop is the same as in OUTPUT
END DO
!
RETURN
END SUBROUTINE UPDATE_MP
