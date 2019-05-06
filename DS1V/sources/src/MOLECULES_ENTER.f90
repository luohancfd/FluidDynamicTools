!
!*****************************************************************************
!
SUBROUTINE MOLECULES_ENTER
!
!--molecules enter boundary at XB(1) and XB(2) and may be removed behind a wave
!
USE MOLECS
USE GAS
USE CALC
USE GEOM
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: K,L,M,N,NENT,II,J,JJ,KK,NTRY,IDT=0 !included IDT
REAL(KIND=8) :: A,B,AA,BB,U,VN,XI,X,DX,DY,DZ,RANF !--isebasti: included RANF
!
!--NENT number to enter in the time step
!
DO J=1,2       !--J is the end
  IF (ITYPE(J) == 0) THEN
    KK=1 !--the entry surface will normally use the reference gas (main stream) properties
    IF ((J == 2).OR.((ISECS == 1).AND.(XB(2) > 0.D00))) KK=2  !--KK is 1 for reference gas 2 for the secondary stream  !--isebasti: AND replaced by OR
    DO L=1,MSP
      A=ENTR(1,L,J)*DTM+ENTR(2,L,J)
      NENT=A
      ENTR(2,L,J)=A-NENT
      IF (NENT > 0) THEN
        DO M=1,NENT
          IF (NM >= MNM) THEN
            WRITE (*,*) 'EXTEND_MNM from MOLECULES_ENTER'
            CALL EXTEND_MNM(1.1d0)
          END IF
          NM=NM+1
          AA=DMAX1(0.D00,ENTR(3,L,J)-3.D00)
          BB=DMAX1(3.D00,ENTR(3,L,J)+3.D00)
          II=0
          DO WHILE (II == 0)
            CALL ZGF(RANF,IDT)
            B=AA+(BB-AA)*RANF
            U=B-ENTR(3,L,J)
            A=(2.D00*B/ENTR(4,L,J))*DEXP(ENTR(5,L,J)-U*U)
            CALL ZGF(RANF,IDT)
            IF (A > RANF) II=1
          END DO
          PV(1,NM)=B*UVMP(L,KK)  !--isebasti: VMP replaced by UVMP
          IF (J == 2) PV(1,NM)=-PV(1,NM)
!
          CALL RVELC(PV(2,NM),PV(3,NM),UVMP(L,KK),IDT)  !--isebasti: VMP replaced by UVMP
          PV(2,NM)=PV(2,NM)+VFY(J)
!
          IF (ISPR(1,L) > 0) CALL SROT(L,UFTMP(KK),PROT(NM),IDT)  !--isebasti: UFTMP replaced by UFTMP
!
          IF (MMVM > 0) THEN
            IPVIB(0,NM)=0
            DO K=1,ISPV(L)
              CALL SVIB(L,UFTMP(KK),IPVIB(K,NM),K,IDT)  !--isebasti: FTMP replaced by UFTMP
            END DO
          END IF
          IPSP(NM)=L
!--advance the molecule into the flow
          CALL ZGF(RANF,IDT)
          XI=XB(J)
          DX=DTM*RANF*PV(1,NM)
          IF (IFX == 0) X=XI+DX
          IF (IFX > 0) DY=DTM*RANF*PV(2,NM)
          DZ=0.D00
          IF (IFX == 2) DZ=DTM*RANF*PV(3,NM)
          IF (IFX > 0) CALL AIFX(XI,DX,DY,DZ,X,PV(1,NM),PV(2,NM),PV(3,NM),IDT)
          IF (IFX == 0) PX(NM)=X
          PTIM(NM)=FTIME
          IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),JJ)
          IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),JJ,PTIM(NM))
          IPCP(NM)=0
          IF (XREM > XB(1)) ENTMASS=ENTMASS+SP(5,L)
        END DO
      END IF
    END DO
  END IF
END DO
!
!--stagnation streamline molecule removal
IF (XREM > XB(1)) THEN
  ENTMASS=FREM*ENTMASS
  NTRY=0
  DO WHILE ((ENTMASS > 0.D00).AND.(NTRY < 10000))
    NTRY=NTRY+1
    IF (NTRY == 10000) THEN
      WRITE (*,*) 'Unable to find molecule for removal'
      ENTMASS=0.D00
      VNMAX=0.D00
    END IF
    CALL ZGF(RANF,IDT)
    N=NM*RANF+0.9999999D00
    IF (PX(N) > XREM) THEN
      CALL ZGF(RANF,IDT)
      !IF (RANF < ((PX(N)-XREM)/(XB(2)-XREM))**2) THEN
      IF (DABS(VFY(1)) < 1.D-3) THEN
        VN=DSQRT(PV(2,N)*PV(2,N)+PV(3,N)*PV(3,N))   !--AXIALLY SYMMETRIC STREAMLINE
      ELSE
        VN=DABS(PV(3,N))   !--TWO-DIMENSIONAL STREAMLINE
      END IF
      L=IPSP(N)
      IF (VN > VNMAX(L)) VNMAX(L)=VN
      CALL ZGF(RANF,IDT)
      IF (RANF < VN/VNMAX(L)) THEN
        CALL REMOVE_MOL(N)
        ENTMASS=ENTMASS-SP(5,L)
        NTRY=0
      END IF
      !END IF
    END IF
  END DO
END IF

!
END SUBROUTINE MOLECULES_ENTER
