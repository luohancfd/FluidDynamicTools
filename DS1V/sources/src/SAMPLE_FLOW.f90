!
!*****************************************************************************
SUBROUTINE SAMPLE_FLOW
!
!--sample the flow properties
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: NC,NCC,LS,N,M,K,L,I,KV
REAL(KIND=8) :: A,TE,TT,WF,EVIB
!
!--NC the sampling cell number
!--NCC the collision cell number
!--LS the species code
!--N,M,K working integers
!--TE total translational energy
!
NSAMP=NSAMP+1
TSAMP=FTIME
!
!$omp parallel &
!$omp private(n,ncc,nc,wf,ls,m,k,evib) &
!$omp reduction(+:cs)
!$omp do schedule(static)
DO N=1,NM
  NCC=IPCELL(N)
  NC=ICCELL(3,NCC)
  WF=1.D00
  IF (IWF == 1) WF=1.D00+WFM*PX(N)**IFX
  IF ((NC > 0).AND.(NC <= NCELLS)) THEN
    IF (MSP > 1) THEN
      LS=ABS(IPSP(N))
    ELSE
      LS=1
    END IF
    CS(0,NC,LS)=CS(0,NC,LS)+1.D00
    CS(1,NC,LS)=CS(1,NC,LS)+WF
    DO M=1,3
      CS(M+1,NC,LS)=CS(M+1,NC,LS)+WF*PV(M,N)
      CS(M+4,NC,LS)=CS(M+4,NC,LS)+WF*PV(M,N)**2
    END DO
    IF (MMRM > 0) CS(8,NC,LS)=CS(8,NC,LS)+WF*PROT(N)
    IF (MMVM > 0) THEN
      IF (ISPV(LS) > 0) THEN
        DO K=1,ISPV(LS)
          CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,LS)
          CS(K+8,NC,LS)=CS(K+8,NC,LS)+WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,LS) !vibrational energy accumulation
        END DO
      END IF
    END IF
  ELSE
    WRITE (*,*) 'Illegal sampling cell',NC,NCC,' for MOL',N,' at',PX(N)  !;STOP
  END IF
END DO
!$omp end do
!$omp end parallel
!
A=0.d0
IF ((IENERS > 0).AND.(GASCODE == 8)) CALL ENERGY(GASCODE,0,A)
ENERS=ENERS+A
!
RETURN
END SUBROUTINE SAMPLE_FLOW
