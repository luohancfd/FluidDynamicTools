!
!*****************************************************************************
!
SUBROUTINE UPDATE_BC
!
!--author: isebasti
!--update inflow boundary conditions
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: NS,K,KK,L,N,NBC,KCELLS(2,10)      !need to set NBC and KCELLS in some commom block
REAL(KIND=8) :: A,B,AA,BB,SN,UBMEAN(2)
!
!--set cells close to boundaries
!
NBC=1  !consider a maximum of NBC sampling cells from each boundary
IF (NCELLS < NBC) NBC=NCELLS
DO N=1,NBC
  KCELLS(1,N)=N
  KCELLS(2,N)=NCELLS-(N-1)
END DO
!
!--take mean velocities weighted by density
!
UBMEAN=0.
!
DO NS=1,2
  A=0.d0
  DO KK=1,NBC
    N=KCELLS(NS,KK)
    UBMEAN(NS)=UBMEAN(NS)+VAR(4,N)*VAR(5,N)
    A=A+VAR(4,N)
  END DO
  UBMEAN(NS)=UBMEAN(NS)/A !=DSUM(rho*u)/DSUM(rho)
END DO
!
!--update inlet macroscopic properties velocities
!
IF (IUBC == 1) THEN !trying to ensure boundary velocities near to VFX (assumed UFND and UFTMP are equal to freestream values)
  A=UBMEAN(1)-VFX(1)
  B=UBMEAN(2)-VFX(2)
!
  IF (DABS(A) >= 50.d0 .AND. DABS(A) < 80.d0) UVFX(1)=UVFX(1)-0.2d0*A
  IF (DABS(B) >= 50.d0 .AND. DABS(B) < 80.d0) UVFX(2)=UVFX(2)-0.2d0*B
!
  IF (DABS(A) >= 10.d0 .AND. DABS(A) < 50.d0) UVFX(1)=UVFX(1)-0.1d0*A
  IF (DABS(B) >= 10.d0 .AND. DABS(B) < 50.d0) UVFX(2)=UVFX(2)-0.1d0*B
!
  IF (DABS(A) >= 00.d0 .AND. DABS(A) < 10.d0) UVFX(1)=UVFX(1)-0.05d0*A
  IF (DABS(B) >= 00.d0 .AND. DABS(B) < 10.d0) UVFX(2)=UVFX(2)-0.05d0*B
END IF
!
IF (IUBC == 2) THEN !pressure based boundary conditions (p and T are known/fixed at boundaries); see pg 66 in my Master's thesis
  DO K=1,2
    IF(K==1) N=1
    IF(K==2) N=NCELLS
!
    A=VAR(17,N)/VAR(12,N)     !interior sound speed
    UFND(K)=FND(K)            !inlet number density
    UFTMP(K)=FTMP(K)          !inlet temperature
    B=UFND(K)*BOLTZ*UFTMP(K)  !inlet pressure
!
    IF (K == 1) UVFX(K)=VAR(5,N)+(B-VAR(18,N))/(VAR(4,N)*A)   !upstream velocity
    IF (K == 2) UVFX(K)=VAR(5,N)-(B-VAR(18,N))/(VAR(4,N)*A)   !downstream velocity
  END DO
END IF
!
!--update inlet most probable speeds
!
DO K=1,2
  DO L=1,MSP
    UVMP(L,K)=SQRT(2.D0*BOLTZ*UFTMP(K)/SP(5,L)) !most probable speed
  END DO
END DO
!
!--update the entry quantities
!
DO K=1,2
  IF (ITYPE(K) == 0) THEN
    DO L=1,MSP
      IF (K == 1) SN=UVFX(1)/UVMP(L,1)
      IF (K == 2) SN=-UVFX(2)/UVMP(L,2)
      AA=SN
      A=1.D00+ERF(AA)
      BB=DEXP(-SN**2)
      ENTR(3,L,K)=SN
      ENTR(4,L,K)=SN+SQRT(SN**2+2.D00)
      ENTR(5,L,K)=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,K)))
      ENTR(6,L,K)=3.D00*UVMP(L,K)
      B=BB+SPI*SN*A
      ENTR(1,L,K)=(UFND(K)*FSP(L,K)*UVMP(L,K))*B/(FNUM*2.D00*SPI)
      ENTR(2,L,K)=0.D00
    END DO
  END IF
END DO
!
RETURN
END SUBROUTINE UPDATE_BC
