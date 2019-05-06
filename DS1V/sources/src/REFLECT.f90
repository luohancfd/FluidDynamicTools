!
!*****************************************************************************
!
SUBROUTINE REFLECT(N,J,X,IDT,LOCAL_CSS, LOCAL_CSSS)
!
!--reflects molecule N and samples the surface J properties
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
INTEGER :: N,J,L,K,IDT !included IDT
REAL(KIND=8) :: A,B,VMPS,DTR,X,XI,DX,DY,DZ,WF,EVIB,RANF !--isebasti: included RANF
REAL(KIND=8),INTENT(OUT) :: LOCAL_CSS(0:8,2), LOCAL_CSSS(6)
!
!--VMPS most probable velocity at the surface temperature
!--DTR time remaining after molecule hits a surface
!
!iflag=0; if ((n == 66430).and.(ftime > 1.5255121136)) iflag=1 !--isebasti: lines for debugging
!

L=IPSP(N)
WF=1.D00
IF (IWF == 1) WF=1.D00+WFM*X**IFX

! A, B are pre-reflection values
A=(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
B=DABS(PV(1,N))
!remove!$omp critical(reflect1)
LOCAL_CSS(0,1)=1.D00
LOCAL_CSS(1,1)=WF
LOCAL_CSS(2,1)=WF*PV(1,N)*SP(5,L)
LOCAL_CSS(3,1)=WF*(PV(2,N)-VSURF(J))*SP(5,L)
LOCAL_CSS(4,1)=WF*PV(3,N)*SP(5,L)
LOCAL_CSS(5,1)=WF*0.5D00*SP(5,L)*A
! CSS(0,J,L,1)=CSS(0,J,L,1)+1.D00
! CSS(1,J,L,1)=CSS(1,J,L,1)+WF
! CSS(2,J,L,1)=CSS(2,J,L,1)+WF*PV(1,N)*SP(5,L)
! CSS(3,J,L,1)=CSS(3,J,L,1)+WF*(PV(2,N)-VSURF(J))*SP(5,L)
! CSS(4,J,L,1)=CSS(4,J,L,1)+WF*PV(3,N)*SP(5,L)
! CSS(5,J,L,1)=CSS(5,J,L,1)+WF*0.5D00*SP(5,L)*A

IF (ISPR(1,L) > 0) THEN
  LOCAL_CSS(6,1)=WF*PROT(N)
ELSE
  LOCAL_CSS(6,1)=0.0d0
END IF
! IF (ISPR(1,L) > 0) CSS(6,J,L,1)=CSS(6,J,L,1)+WF*PROT(N)

LOCAL_CSS(7,1) = 0.0d0
IF (MMVM > 0) THEN
  IF (ISPV(L) > 0) THEN
    DO K=1,ISPV(L)
      CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,L)
      LOCAL_CSS(7,1)=LOCAL_CSS(7,1)+WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
      ! CSS(7,J,L,1)=CSS(7,J,L,1)+WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
    END DO
  END IF
END IF
LOCAL_CSS(8,1) = 0.0d0


LOCAL_CSSS(1)=WF/B
LOCAL_CSSS(2)=WF*SP(5,L)/B
LOCAL_CSSS(3)=WF*SP(5,L)*PV(2,N)/B
!--this assumes that any flow normal to the x direction is in the y direction
LOCAL_CSSS(4)=WF*SP(5,L)*A/B
! CSSS(1,J)=CSSS(1,J)+WF/B
! CSSS(2,J)=CSSS(2,J)+WF*SP(5,L)/B
! CSSS(3,J)=CSSS(3,J)+WF*SP(5,L)*PV(2,N)/B
! !--this assumes that any flow normal to the x direction is in the y direction
! CSSS(4,J)=CSSS(4,J)+WF*SP(5,L)*A/B


IF (ISPR(1,L) > 0) THEN
  LOCAL_CSSS(5)=WF*PROT(N)/B
  LOCAL_CSSS(6)=WF*ISPR(1,L)/B
ELSE
  LOCAL_CSSS(5)=0.0d0
  LOCAL_CSSS(6)=0.0d0
END IF
! IF (ISPR(1,L) > 0) THEN
!   CSSS(5,J)=CSSS(5,J)+WF*PROT(N)/B
!   CSSS(6,J)=CSSS(6,J)+WF*ISPR(1,L)/B
! END IF
!remove!$omp end critical(reflect1)
!
CALL ZGF(RANF,IDT)
IF (FSPEC(J) > RANF) THEN !--specular reflection
  X=2.D00*XB(J)-X
  PV(1,N)=-PV(1,N)
  DTR=(X-XB(J))/PV(1,N)
ELSE                      !--diffuse reflection
  VMPS=SQRT(2.D00*BOLTZ*TSURF(J)/SP(5,L))
  DTR=(X-XB(J))/PV(1,N)   !--isebasti: original line DTR=(XB(J)-PX(N))/PV(1,N) was corrected
  !write(*,*) n,px(n),x-px(n),x,dtm,dtr
  CALL ZGF(RANF,IDT)
  PV(1,N)=SQRT(-LOG(RANF))*VMPS
  IF (J == 2) PV(1,N)=-PV(1,N)
  CALL RVELC(PV(2,N),PV(3,N),VMPS,IDT)
  PV(2,N)=PV(2,N)+VSURF(J)
  IF (ISPR(1,L) > 0) CALL SROT(L,TSURF(J),PROT(N),IDT)
  IF (MMVM > 0) THEN
    DO K=1,ISPV(L)
      CALL SVIB(L,TSURF(J),IPVIB(K,N),K,IDT)
    END DO
  END IF
END IF
!
! A, B are now post reflection value
A=(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
B=DABS(PV(1,N))
!remove!$omp critical(reflect2)
LOCAL_CSS(2,2)=-WF*PV(1,N)*SP(5,L)
LOCAL_CSS(3,2)=-WF*(PV(2,N)-VSURF(J))*SP(5,L)
LOCAL_CSS(4,2)=-WF*PV(3,N)*SP(5,L)
LOCAL_CSS(5,2) = -WF*0.5D00*SP(5,L)*A
! CSS(2,J,L,2)=CSS(2,J,L,2)-WF*PV(1,N)*SP(5,L)
! CSS(3,J,L,2)=CSS(3,J,L,2)-WF*(PV(2,N)-VSURF(J))*SP(5,L)
! CSS(4,J,L,2)=CSS(4,J,L,2)-WF*PV(3,N)*SP(5,L)
! CSS(5,J,L,2)=CSS(5,J,L,2)-WF*0.5D00*SP(5,L)*A

IF (ISPR(1,L) > 0) THEN
  LOCAL_CSS(6,2) = -WF*PROT(N)
ELSE
  LOCAL_CSS(6,2) = 0.0d0
END IF
! IF (ISPR(1,L) > 0) CSS(6,J,L,2)=CSS(6,J,L,2)-WF*PROT(N)

LOCAL_CSS(7,2)=0.0D0
IF (MMVM > 0) THEN
  IF (ISPV(L).GT.0) THEN
    DO K=1,ISPV(L)
      CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,L)
      LOCAL_CSS(7,2)=LOCAL_CSS(7,2)-WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
      ! CSS(7,J,L,2)=CSS(7,J,L,2)-WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
    END DO
  END IF
END IF

LOCAL_CSSS(1)=LOCAL_CSSS(1)+WF/B
LOCAL_CSSS(2)=LOCAL_CSSS(2)+WF*SP(5,L)/B
LOCAL_CSSS(3)=LOCAL_CSSS(3)+WF*SP(5,L)*PV(2,N)/B
!--this assumes that any flow normal to the x direction is in the y direction
LOCAL_CSSS(4)=LOCAL_CSSS(4)+WF*SP(5,L)*A/B
! CSSS(1,J)=CSSS(1,J)+WF/B
! CSSS(2,J)=CSSS(2,J)+WF*SP(5,L)/B
! CSSS(3,J)=CSSS(3,J)+WF*SP(5,L)*PV(2,N)/B
! !--this assumes that any flow normal to the x direction is in the y direction
! CSSS(4,J)=CSSS(4,J)+WF*SP(5,L)*A/B

IF (ISPR(1,L) > 0) THEN
  LOCAL_CSSS(5)=LOCAL_CSSS(5)+WF*PROT(N)/B
  LOCAL_CSSS(6)=LOCAL_CSSS(6)+WF*ISPR(1,L)/B
  ! CSSS(5,J)=CSSS(5,J)+WF*PROT(N)/B  !Fix Israel's bug
  ! CSSS(6,J)=CSSS(6,J)+WF*ISPR(1,L)/B
END IF
!remove!$omp end critical(reflect2)
!
XI=XB(J)
DX=DTR*PV(1,N)
DZ=0.D00
IF (IFX > 0) DY=DTR*PV(2,N)
IF (IFX == 2) DZ=DTR*PV(3,N)
IF (IFX == 0) X=XI+DX
IF (IFX > 0) CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
!write(*,*) n,xi,x-xi,x,dtm,dtr
!pause
!
RETURN
!
END SUBROUTINE REFLECT
