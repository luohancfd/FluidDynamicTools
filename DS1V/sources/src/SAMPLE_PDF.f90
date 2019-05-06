!
!*****************************************************************************
!
SUBROUTINE SAMPLE_PDF
!
!--author: isebasti
!--sample pdfs
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
REAL(KIND=8) :: A,VMPL,SC(5)
INTEGER(KIND=8) :: NS,NCI,NCF,NC,NMC,J,K,N,L,M,I,JB,KV
!
!--NS cell to be sampled
!--NBINS number of bins
!--BINS(0,I,L) sampled sum for property I and species L
!--BIN(0,I) total sampled sum for property I
!
!--set auxiliar variables
!
DBINV(1)=-3.d0 ; DBINV(2)=3.d0  !normalized velocity bin lower/upper boundaries
DBINC(1)= 0.d0 ; DBINC(2)=3.d0  !normalized speed
DBINE(1)= 0.d0 ; DBINE(2)=10.d0 !thermal energy in eV
!
DBINV(3)=(DBINV(2)-DBINV(1))/DFLOAT(NBINS) !bin intervals
DBINC(3)=(DBINC(2)-DBINC(1))/DFLOAT(NBINS)
DBINE(3)=(DBINE(2)-DBINE(1))/DFLOAT(NBINS)
!
!--sum to corresponding bin
!
NSPDF=NCELLS/2.+0.9999999d0       !user defined; must be consistent with grid
NS=NSPDF                          !sampling cell
NCI=ICELL(NS)+1                   !initial collision cell
IF(NS < NCELLS)  NCF=ICELL(NS+1)  !final collision cell
IF(NS == NCELLS) NCF=NCCELLS
IF(NCELLS == 1) NCI=1             !debuging; it can be improved
!
DO NC=NCI,NCF          !do over collision cells
  NMC=ICCELL(2,NC)     !#of molecules in collision cell NC
!
  DO J=1,NMC           !do over molecules in collision cell NC
    K=J+ICCELL(1,NC)   !adress of molecule j in ICREF
    N=ICREF(K)         !molecule number
    L=IPSP(N)          !species of molecule n
!
    VMPL=SQRT(2.D00*BOLTZ*VARSP(5,NS,L)/SP(5,L)) !most probable velocity species L
!
    DO I=1,3
      SC(I)=(PV(I,N)-VAR(I+4,NS))/VMPL                      !normalized thermal velocity
      JB=1+(SC(I)-DBINV(1))/DBINV(3)                        !bin index
      IF (JB>0.AND.JB<NBINS) BINS(JB,I,L)=BINS(JB,I,L)+1.d0 !bin accumulation per species
    END DO
!
    SC(4)=DSQRT(SC(1)*SC(1)+SC(2)*SC(2)+SC(3)*SC(3))        !normalized thermal speed
    JB=1+(SC(4)-DBINC(1))/DBINC(3)                          !bin index
    IF (JB>0.AND.JB<NBINS) BINS(JB,4,L)=BINS(JB,4,L)+1.d0   !bin accumulation per species
!
    A=0.d0
    DO I=1,3
      A=A+(PV(I,N)-VAR(I+4,NS))**2.d0
    END DO
    SC(5)=0.5d0*SP(5,L)*A/EVOLT                    !translational thermal energy (eV)
    JB=1+(SC(5)-DBINE(1))/DBINE(3)                 !bin index
    IF (JB>NBINS) JB=NBINS
    BINS(JB,5,L)=BINS(JB,5,L)+1.d0                 !bin accumulation per species
!
    IF (ISPR(1,L) > 0) THEN                        !rotational energy
      M=1+(PROT(N)/(BOLTZ*VARSP(6,NS,L)))/0.1D00   !bin index
      IF (M < 101) NDROT(L,M)=NDROT(L,M)+1         !bin accumulation per species
    END IF
!
    BIN(0,:)=BIN(0,:)+1.d0          !total counter
    BINS(0,:,L)=BINS(0,:,L)+1.d0    !species counter
  END DO
END DO
!
!--sum to corresponding bin (vibrational levels)
!
DO I=1,NSCELLS
  NS=NSVEC(I)                       !sampling cell
  NCI=ICELL(NS)+1                   !initial collision cell
  IF(NS < NCELLS)  NCF=ICELL(NS+1)  !final collision cell
  IF(NS == NCELLS) NCF=NCCELLS
  IF(NCELLS == 1) NCI=1             !debuging; it can be improved
!
  DO NC=NCI,NCF          !do over collision cells
    NMC=ICCELL(2,NC)     !#of molecules in collision cell NC
!
    DO J=1,NMC           !do over molecules in collision cell NC
      K=J+ICCELL(1,NC)   !adress of molecule j in ICREF
      N=ICREF(K)         !molecule number
      L=IPSP(N)          !species of molecule n
!
      IF (ISPV(L) > 0) THEN
        DO KV=1,ISPV(L)                                           !vibrational mode
          NDVIB(NS,KV,L,IPVIB(KV,N))=1+NDVIB(NS,KV,L,IPVIB(KV,N)) !bin accumulation per species
        END DO
      END IF
!
      NDVIB(NS,0,L,0)=1+NDVIB(NS,0,L,0)   !total counter
    END DO
!
  END DO
END DO
!
RETURN
END SUBROUTINE SAMPLE_PDF
