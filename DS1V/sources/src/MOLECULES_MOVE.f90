!
!*****************************************************************************
!
SUBROUTINE MOLECULES_MOVE
!
!--molecule moves appropriate to the time step
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: N,L,M,K,NCI,J,II,JJ,IDT=0 !included IDT,iflag
REAL(KIND=8) :: A,B,X,XI,XC,DX,DY,DZ,DTIM,S1,XM,R,TI,DTC,POB,UR,WFI,WFR,WFRI,RANF !--isebasti: included RANF
REAL(KIND=8) :: LOCAL_CSS(0:8,2), LOCAL_CSSS(6)
INTEGER :: LOCAL_ISP
!
!--N working integer
!--NCI initial collision cell
!--DTIM time interval for the move
!--POB position of the outer boundary
!--TI initial time
!--DTC time interval to collision with surface
!--UR radial velocity component
!--WFI initial weighting factor
!--WFR weighting factor radius
!--WFRI initial weighting factor radius
!
!$OMP PARALLEL &
!$OMP& PRIVATE(IDT,N,NCI,DTIM,WFI,II,TI,XI,DX,X,DZ,DY,R,J,L,DTC,XC,S1,WFR,WFRI,K,M,JJ) &
!$OMP& PRIVATE(LOCAL_CSS,LOCAL_CSSS,LOCAL_ISP) &
!$OMP& DEFAULT(SHARED) &
!$OMP& REDUCTION(+:TOTMOV,ENTMASS,CSS,CSSS)
!$    IDT=OMP_GET_THREAD_NUM()  !THREAD ID
!$OMP DO SCHEDULE (STATIC)
!
DO N=1,NM
  LOCAL_CSS=0.0d0
  LOCAL_CSSS=0.0d0
!
  NCI=IPCELL(N)
  LOCAL_ISP=IPSP(N)
  IF (IMTS == 0) DTIM=DTM
  IF (IMTS == 1) DTIM=2.D00*CCELL(3,NCI)
  IF (FTIME-PTIM(N) > 0.5*DTIM) THEN
    WFI=1.D00
    IF (IWF == 1) WFI=1.D00+WFM*PX(N)**IFX
    II=0 !--becomes 1 if a molecule is removed
    TI=PTIM(N)
    PTIM(N)=TI+DTIM
    TOTMOV=TOTMOV+1
!
    XI=PX(N)
    DX=DTIM*PV(1,N)
    X=XI+DX
!
    IF (IFX > 0) THEN
      DZ=0.D00
      DY=DTIM*PV(2,N)
      IF (IFX == 2) DZ=DTIM*PV(3,N)
      R=DSQRT(X*X+DY*DY+DZ*DZ)
    END IF
!
    IF (IFX == 0) THEN
      DO J=1,2    ! 1 for minimum x boundary, 2 for maximum x boundary
        IF (II == 0) THEN
          IF (((J == 1).AND.(X < XB(1))).OR.((J == 2).AND.(X > (XB(2)+VELOB*PTIM(N))))) THEN  !--molecule crosses a boundary
            IF ((ITYPE(J) == 0).OR.(ITYPE(J) == 3)) THEN
              IF (XREM > XB(1)) THEN
                ENTMASS=ENTMASS-SP(5,LOCAL_ISP)
              END IF
              IPCELL(N)=-IPCELL(N) !molecule is marked for removel !--isebasti: original use CALL REMOVE_MOL(N); !N=N-1
              II=1
            END IF
!
            IF (ITYPE(J) == 1) THEN
              IF ((IVB == 0).OR.(J == 1)) THEN
                X=2.D00*XB(J)-X
                PV(1,N)=-PV(1,N)
              ELSE IF ((J == 2).AND.(IVB == 1)) THEN
                DTC=(XB(2)+TI*VELOB-XI)/(PV(1,N)-VELOB)
                XC=XI+PV(1,N)*DTC
                PV(1,N)=-PV(1,N)+2.*VELOB
                X=XC+PV(1,N)*(DTIM-DTC)
              END IF
            END IF
!
            IF (ITYPE(J) == 2) THEN
              CALL REFLECT(N,J,X,IDT,LOCAL_CSS,LOCAL_CSSS)
              DO L=1,2
                DO K=0,8
                  CSS(K,J,LOCAL_ISP,L) = CSS(K,J,LOCAL_ISP,L) + LOCAL_CSS(K,L)
                END DO
              END DO
              CSSS(:,J) = CSSS(:,J) + LOCAL_CSSS

                IF((J == 2).AND.(X < XB(1))) THEN !--isebasti: included to fix bug
                  DO WHILE ((X<XB(1)).OR.(X>XB(2)))
                   !WRITE(9,*) 'REFLECTED MOLEC CROSSING OPPOSITE BOUNDARY:',FTIME,N,X,PV(1,N)
                    IF (X < XB(1)) THEN
                      CALL REFLECT(N,1,X,IDT,LOCAL_CSS, LOCAL_CSSS) !molecules reflected at xb2 crossed xb1
                      DO L=1,2
                          CSS(:,1,LOCAL_ISP,L) = CSS(:,1,LOCAL_ISP,L) + LOCAL_CSS(:,L)
                      END DO
                      CSSS(:,1) = CSSS(:,1) + LOCAL_CSSS
                    END IF
                    IF (X > XB(2)) THEN
                      CALL REFLECT(N,2,X,IDT,LOCAL_CSS,LOCAL_CSSS) !molecules reflected at xb1 crossed xb2
                      DO L=1,2
                          CSS(:,2,LOCAL_ISP,L) = CSS(:,2,LOCAL_ISP,L) + LOCAL_CSS(:,L)
                      END DO
                      CSSS(:,2) = CSSS(:,2) + LOCAL_CSSS
                    END IF
                   !WRITE(9,*) 'REFLECTED MOLEC CROSSING OPPOSITE BOUNDARY:',FTIME,N,X,PV(1,N)
                  END DO
                END IF
            END IF
          END IF
        END IF
      END DO
    ELSE         !--cylindrical or spherical flow
!--check boundaries
      IF ((X < XB(1)).AND.(XB(1) > 0.D00)) THEN
        CALL RBC(XI,DX,DY,DZ,XB(1),S1)
        IF (S1 < 1.D00) THEN     !--intersection with inner boundary
          IF (ITYPE(1) == 2) THEN !--solid surface
            DX=S1*DX
            DY=S1*DY
            DZ=S1*DZ
            CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
            CALL REFLECT(N,1,X,IDT,LOCAL_CSS, LOCAL_CSSS)
            DO L=1,2
                CSS(:,1,LOCAL_ISP,L) = CSS(:,1,LOCAL_ISP,L) + LOCAL_CSS(:,L)
            END DO
            CSSS(:,1) = CSSS(:,1) + LOCAL_CSSS
          ELSE
            IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
            II=1
          END IF
        END IF
      ELSE IF ((IVB == 0).AND.(R > XB(2))) THEN
        CALL RBC(XI,DX,DY,DZ,XB(2),S1)
        IF (S1 < 1.D00) THEN     !--intersection with outer boundary
          IF (ITYPE(2) == 2) THEN !--solid surface
            DX=S1*DX
            DY=S1*DY
            DZ=S1*DZ
            CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
            X=1.001D00*XB(2)
            DO WHILE (X > XB(2))
              CALL REFLECT(N,2,X,IDT,LOCAL_CSS,LOCAL_CSSS)
              DO L=1,2
                  CSS(:,2,LOCAL_ISP,L) = CSS(:,2,LOCAL_ISP,L) + LOCAL_CSS(:,L)
              END DO
              CSSS(:,2) = CSSS(:,2) + LOCAL_CSSS
            END DO
          ELSE
            IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
            II=1
          END IF
        END IF
      ELSE IF ((IVB == 1).AND.(R > (XB(2)+PTIM(N)*VELOB))) THEN
        IF (IFX == 1) UR=DSQRT(PV(1,N)**2+PV(2,N)**2)
        IF (IFX == 2) UR=DSQRT(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        DTC=(XB(2)+TI*VELOB-XI)/(UR-VELOB)
        S1=DTC/DTIM
        DX=S1*DX
        DY=S1*DY
        DZ=S1*DZ
        CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
        PV(1,N)=-PV(1,N)+2.*VELOB
        X=X+PV(1,N)*(DTIM-DTC)
      ELSE
        CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
      END IF
!
!--DIAGNOSTIC
      IF (X > XB(2)+PTIM(N)*VELOB) THEN
        WRITE (*,*) N,FTIME,X,XB(2)+PTIM(N)*VELOB
      END IF
!
!--Take action on weighting factors
      IF ((IWF == 1).AND.(II >= 0)) THEN
        WFR=WFI/(1.D00+WFM*X**IFX)
        L=0
        WFRI=WFR
        IF (WFR >= 1.D00) THEN
          DO WHILE (WFR >= 1.D00)
            L=L+1
            WFR=WFR-1.D00
          END DO
        END IF
        CALL ZGF(RANF,IDT)
        IF (RANF <= WFR) L=L+1
        IF (L == 0) THEN
          IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
          II=1
        END IF
        L=L-1
        IF (L > 0) THEN
          DO K=1,L
!$omp critical
            IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
            NM=NM+1  !--isebasti: for spherical and cylindrical geometries, it will cause a problem because NM value changes.
            PX(NM)=X
            DO M=1,3
              PV(M,NM)=PV(M,N)
            END DO
            IF (MMRM > 0) PROT(NM)=PROT(N)
            IPCELL(NM)=ABS(IPCELL(N))  !--isebasti: I think we should remove ABS
            IPSP(NM)=IPSP(N)
            IPCP(NM)=IPCP(N)
            IF (MMVM > 0) THEN
              DO M=1,MMVM
                IPVIB(M,NM)=IPVIB(M,N)
              END DO
            END IF
            PTIM(NM)=PTIM(N)    !+5.D00*DFLOAT(K)*DTM
!--note the possibility of a variable time advance that may take the place of the duplication buffer in earlier programs
!
            IF (PX(NM) > XB(2)+PTIM(NM)*VELOB) THEN
              WRITE (*,*) 'DUP',NM,FTIME,PX(NM),XB(2)+PTIM(NM)*VELOB
            END IF
!$omp end critical
          END DO
        END IF
      END IF
    END IF
!
    IF (II == 0) THEN
      PX(N)=X
        if ((px(n) > xb(1)).and.(px(n) < xb(2))) then
          continue
        else
!$omp critical(move)
          write (*,*) n,'Outside flowfield at',px(n)
!$omp end critical(move)
        end if
      IF (IVB == 0) CALL FIND_CELL(PX(N),IPCELL(N),JJ)
      IF (IVB == 1) CALL FIND_CELL_MB(PX(N),IPCELL(N),JJ,PTIM(N))
    END IF
!
  END IF
!
END DO
!$omp end do
!$omp end parallel
!
!--isebasti: remove marked molecules
N=0
DO WHILE (N < NM)
  N=N+1
  IF (IPCELL(N) < 0) THEN
    CALL REMOVE_MOL(N)
    N=N-1
  END IF
END DO
!
RETURN
!
END SUBROUTINE MOLECULES_MOVE
