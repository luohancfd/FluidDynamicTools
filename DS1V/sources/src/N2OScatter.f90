!
!*****************************************************************************
SUBROUTINE N2OScatter(LS,MS,BMAX, VRI,VR,VRC,VRCP,IDT)
!
! -- calculate new scattering angles based on non VHS model
! -- only for N2+O collision
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
INTEGER :: LS,MS,IDT
REAL(KIND=8) :: B,VR,VRC(3),VRCP(3),RANF,C,D,VRI
REAL(8) :: COF(6), BMAX, ET, RMASS, EVIBEV,C1,C2,CHI,CCHI,SCHI
REAL(8) :: EPSI,CEPSI,SEPSI

COF(1) = 1.783071942310433
COF(2) = 0.016980219000487/1000.0D0
COF(3) = 2.846002511624816
COF(4) = 0.294962630947434/1000.0D0
COF(5) = 2.189954664950628
COF(6) = 0.078701684994745


EVIBEV =  0.063577602025999  !Erv(v=0,J=0) for N2

CALL ZGF(RANF,IDT)
B = DSQRT(RANF)*BMAX
C1 = COF(3)*DEXP(-COF(4)*VRI) +COF(5)*DEXP(-COF(6)*EVIBEV)
C2 = COF(1) + COF(2)*VRI
CHI = C2*(1.0d0 - 1.0d0/(C1+DLOG(2.0D0))*(B-DLOG(DCOSH(B-C1))))
! Scattering angle should only be a function of precollision conditions
! this two condition shouldn't occur, for safety we set it here
IF (CHI > PI) CHI = PI
IF (CHI < 0) CHI = 0.0D0
CCHI = DCOS(CHI); SCHI = DSIN(CHI)

CALL ZGF(RANF,IDT)
EPSI = RANF*2.0D0*PI
CEPSI = DCOS(EPSI)
SEPSI = DSIN(EPSI)

D=DSQRT(VRC(2)**2+VRC(3)**2)
VRCP(1) = CCHI*VRC(1) + SCHI*SEPSI*D
VRCP(2) = CCHI*VRC(2) + SCHI*(VR*VRC(3)*CEPSI-VRC(1)*VRC(2)*SEPSI)/D
VRCP(3) = CCHI*VRC(3) - SCHI*(VR*VRC(2)*CEPSI+VRC(1)*VRC(3)*SEPSI)/D

RETURN
END SUBROUTINE N2OScatter
