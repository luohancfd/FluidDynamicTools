!
!***************************************************************************
!
SUBROUTINE REAL_AIR
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=5
MMRM=1
MMVM=1

!
IF (JCD == 0) THEN
  MNRE=23
  MTBP=4
  MEX=0
  MMEX=0
END IF
IF (JCD == 1) THEN
  MNRE=0
  MTBP=0
  MEX=4
  MMEX=1
END IF
!
MNSR=0
CALL ALLOCATE_GAS
!--species 1 is oxygen
SP(1,1)=4.07D-10
SP(2,1)=273.D00
SP(3,1)=0.77D00
SP(4,1)=1.d00
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1               ! the number of vibrational modes
SPVM(1,1,1)=2256.D00          ! the characteristic vibrational temperature
SPVM(2,1,1)=18000.D00  !90000.D00        ! a constant Zv, or the reference Zv
SPVM(3,1,1)=2256.D00        ! -1 for a constant Zv, or the reference temperature
SPVM(4,1,1)=59500.D00
ISPVM(1,1,1)=3
ISPVM(2,1,1)=3
!--species 2 is nitrogen
SP(1,2)=4.17D-10
SP(2,2)=273.D00
SP(3,2)=0.74D00
SP(4,2)=1.D00
SP(5,2)=4.65D-26
ISPR(1,2)=2
ISPR(2,2)=0
SPR(1,2)=5.D00
ISPV(2)=1
SPVM(1,1,2)=3371.D00
SPVM(2,1,2)=52000.D00     !260000.D00
SPVM(3,1,2)=3371.D00
SPVM(4,1,2)=113500.D00
ISPVM(1,1,2)=4
ISPVM(2,1,2)=4
!--species 3 is atomic oxygen
SP(1,3)=3.D-10
SP(2,3)=273.D00
SP(3,3)=0.8D00
SP(4,3)=1.D00
SP(5,3)=2.656D-26
ISPR(1,3)=0
ISPV(3)=0
!--species 4 is atomic nitrogen
SP(1,4)=3.D-10
SP(2,4)=273.D00
SP(3,4)=0.8D00
SP(4,4)=1.0D00
SP(5,4)=2.325D-26
ISPR(1,4)=0
ISPV(4)=0
!--species 5 is NO
SP(1,5)=4.2D-10
SP(2,5)=273.D00
SP(3,5)=0.79D00
SP(4,5)=1.0D00
SP(5,5)=4.98D-26
ISPR(1,5)=2
ISPR(2,5)=0
SPR(1,5)=5.D00
ISPV(5)=1
SPVM(1,1,5)=2719.D00
SPVM(2,1,5)=14000.D00   !70000.D00
SPVM(3,1,5)=2719.D00
SPVM(4,1,5)=75500.D00
ISPVM(1,1,5)=3
ISPVM(2,1,5)=4
!--following data is required if JCD=1 (new reaction model)
IF (JCD == 1) THEN
!--set the recombination data for the molecule pairs
  ISPRC=0    !--data os zero unless explicitly set
  SPRC=0.D00
  ISPRC(3,3)=1    !--O+O -> O2  recombined species code for an O+O recombination
  SPRC(1,3,3,1)=0.04D00
  SPRC(2,3,3,1)=-1.3D00
  SPRC(1,3,3,2)=0.07D00
  SPRC(2,3,3,2)=-1.2D00
  SPRC(1,3,3,3)=0.08D00
  SPRC(2,3,3,3)=-1.2D00
  SPRC(1,3,3,4)=0.09D00
  SPRC(2,3,3,4)=-1.2D00
  SPRC(1,3,3,5)=0.065D00
  SPRC(2,3,3,5)=-1.2D00
  ISPRC(4,4)=2
  SPRC(1,4,4,1)=0.15D00
  SPRC(2,4,4,1)=-2.05D00
  SPRC(1,4,4,2)=0.09D00
  SPRC(2,4,4,2)=-2.1D00
  SPRC(1,4,4,3)=0.16D00
  SPRC(2,4,4,3)=-2.0D00
  SPRC(1,4,4,4)=0.17D00
  SPRC(2,4,4,4)=-2.0D00
  SPRC(1,4,4,5)=0.17D00
  SPRC(2,4,4,5)=-2.1D00
  ISPRC(3,4)=5
  SPRC(1,3,4,1)=0.3D00
  SPRC(2,3,4,1)=-1.9D00
  SPRC(1,3,4,2)=0.4D00
  SPRC(2,3,4,2)=-2.0D00
  SPRC(1,3,4,3)=0.3D00
  SPRC(2,3,4,3)=-1.75D00
  SPRC(1,3,4,4)=0.3D00
  SPRC(2,3,4,4)=-1.75D00
  SPRC(1,3,4,5)=0.15D00
  SPRC(2,3,4,5)=-1.9D00
  ISPRC(4,3)=5
  SPRC(1,4,3,1)=0.3D00
  SPRC(2,4,3,1)=-1.9D00
  SPRC(1,4,3,2)=0.4D00
  SPRC(2,4,3,2)=-2.0D00
  SPRC(1,4,3,3)=0.3D00
  SPRC(2,4,3,3)=-1.75D00
  SPRC(1,4,3,4)=0.3D00
  SPRC(2,4,3,4)=-1.75D00
  SPRC(1,4,3,5)=0.15D00
  SPRC(2,4,3,5)=-1.9D00
!--set the exchange reaction data
  SPEX=0.D00
  ISPEX=0
  NSPEX=0
  NSPEX(2,3)=1
  NSPEX(3,2)=1
  NSPEX(5,4)=1
  NSPEX(4,5)=1
  NSPEX(5,3)=1
  NSPEX(3,5)=1
  NSPEX(1,4)=1
  NSPEX(4,1)=1
!--N2+O->NO+N
  ISPEX(1,0,2,3)=2
  ISPEX(1,1,2,3)=5
  ISPEX(1,2,2,3)=4
  ISPEX(1,3,2,3)=1
  ISPEX(1,4,2,3)=1
  SPEX(1,1,2,3)=0.15D00
  SPEX(2,1,2,3)=0.D00
  SPEX(3,1,2,3)=-5.175D-19
  NEX(1,2,3)=1
  ISPEX(1,0,3,2)=2
  ISPEX(1,1,3,2)=5
  ISPEX(1,2,3,2)=4
  ISPEX(1,3,3,2)=1
  ISPEX(1,4,3,2)=1
  SPEX(1,1,3,2)=0.15D00
  SPEX(2,1,3,2)=0.D00
  SPEX(3,1,3,2)=-5.175D-19
  NEX(1,3,2)=1
!--NO+N->N2+0
  ISPEX(1,0,5,4)=5
  ISPEX(1,1,5,4)=2
  ISPEX(1,2,5,4)=3
  ISPEX(1,3,5,4)=1
  ISPEX(1,4,5,4)=1
  SPEX(1,1,5,4)=0.033D00
  SPEX(2,1,5,4)=0.8D00
  SPEX(3,1,5,4)=5.175D-19
  NEX(1,5,4)=2
  ISPEX(1,0,4,5)=5
  ISPEX(1,1,4,5)=2
  ISPEX(1,2,4,5)=3
  ISPEX(1,3,4,5)=1
  ISPEX(1,4,4,5)=1
  SPEX(1,1,4,5)=0.033D00
  SPEX(2,1,4,5)=0.8D00
  SPEX(3,1,4,5)=5.175D-19
  NEX(1,4,5)=2
!--NO+0->O2+N
  ISPEX(1,0,5,3)=5
  ISPEX(1,1,5,3)=1
  ISPEX(1,2,5,3)=4
  ISPEX(1,3,5,3)=1
  ISPEX(1,4,5,3)=1
  SPEX(1,1,5,3)=0.05D00
  SPEX(2,1,5,3)=0.7D00
  SPEX(3,1,5,3)=-2.719D-19
  NEX(1,5,3)=3
  ISPEX(1,0,3,5)=5
  ISPEX(1,1,3,5)=1
  ISPEX(1,2,3,5)=4
  ISPEX(1,3,3,5)=1
  ISPEX(1,4,3,5)=1
  SPEX(1,1,3,5)=0.05D00
  SPEX(2,1,3,5)=0.7D00
  SPEX(3,1,3,5)=-2.719D-19
  NEX(1,3,5)=3
!--O2+N->NO+O
  ISPEX(1,0,1,4)=1
  ISPEX(1,1,1,4)=5
  ISPEX(1,2,1,4)=3
  ISPEX(1,3,1,4)=1
  ISPEX(1,4,1,4)=1
  SPEX(1,1,1,4)=0.D00
  SPEX(2,1,1,4)=0.D00
  SPEX(3,1,1,4)=2.719D-19
  NEX(1,1,4)=4
  ISPEX(1,0,4,1)=1
  ISPEX(1,1,4,1)=5
  ISPEX(1,2,4,1)=3
  ISPEX(1,3,4,1)=1
  ISPEX(1,4,4,1)=1
  SPEX(1,1,4,1)=0.D00
  SPEX(2,1,4,1)=0.D00
  SPEX(3,1,4,1)=2.719D-19
  NEX(1,4,1)=4
!
END IF
IF (JCD == 0) THEN
!--the following data is required if JCD=0 (old reaction model)
! REACTION 1 IS O2+N->2O+N
  LE(1)=1
  ME(1)=4
  KP(1)=3
  LP(1)=3
  MP(1)=4
  CI(1)=1.
  AE(1)=8.197D-19
  AC(1)=5.993D-12
  BC(1)=-1.
  ER(1)=-8.197D-19
!--REACTION 2 IS O2+NO>2O+NO
  LE(2)=1
  ME(2)=5
  KP(2)=3
  LP(2)=3
  MP(2)=5
  CI(2)=1.
  AE(2)=8.197D-19
  AC(2)=5.993D-12
  BC(2)=-1.
  ER(2)=-8.197D-19
!--REACTION 3 IS O2+N2>2O+N2
LE(3)=1
ME(3)=2
KP(3)=3
LP(3)=3
MP(3)=2
CI(3)=1.5
AE(3)=8.197D-19
AC(3)=1.198D-11
BC(3)=-1.
ER(3)=-8.197D-19
!--REACTION 4 IS 2O2>2O+O2
LE(4)=1
ME(4)=1
KP(4)=3
LP(4)=3
MP(4)=1
CI(4)=1.5
AE(4)=8.197D-19
AC(4)=5.393D-11
BC(4)=-1.
ER(4)=-8.197D-19
!--REACTION 5 IS O2+O>3O
LE(5)=1
ME(5)=3
KP(5)=3
LP(5)=3
MP(5)=3
CI(5)=1.
AE(5)=8.197D-19
AC(5)=1.498D-10
BC(5)=-1.
ER(5)=-8.197D-19
!--REACTION 6 IS N2+O>2N+O
LE(6)=2
ME(6)=3
KP(6)=4
LP(6)=4
MP(6)=3
CI(6)=0.5
AE(6)=1.561D-18
AC(6)=3.187D-13
BC(6)=-0.5
ER(6)=-1.561D-18
!--REACTION 7 IS N2+O2>2N+O2
LE(7)=2
ME(7)=1
KP(7)=4
LP(7)=4
MP(7)=1
CI(7)=0.5
AE(7)=1.561D-18
AC(7)=3.187D-13
BC(7)=-0.5
ER(7)=-1.561D-18
!--REACTION 8 IS N2+NO>2N+NO
LE(8)=2
ME(8)=5
KP(8)=4
LP(8)=4
MP(8)=5
CI(8)=0.5
AE(8)=1.561D-18
AC(8)=3.187D-13
BC(8)=-0.5
ER(8)=-1.561D-18
!--REACTION 9 IS 2N2>2N+N2
LE(9)=2
ME(9)=2
KP(9)=4
LP(9)=4
MP(9)=2
CI(9)=1.
AE(9)=1.561D-18
AC(9)=7.968D-13
BC(9)=-0.5
ER(9)=-1.561D-18
!--REACTION 10 IS N2+N>3N
LE(10)=2
ME(10)=4
KP(10)=4
LP(10)=4
MP(10)=4
CI(10)=1.
AE(10)=1.561D-18
AC(10)=6.9E12
BC(10)=-1.5
ER(10)=-1.561D-18
!--REACTION 11 IS NO+N2>N+O+N2
LE(11)=5
ME(11)=2
KP(11)=4
LP(11)=3
MP(11)=2
CI(11)=1.
AE(11)=1.043D-18
AC(11)=6.59D-10
BC(11)=-1.5
ER(11)=-1.043D-18
!--REACTION 12 IS NO+O2>N+O+O2
LE(12)=5
ME(12)=1
KP(12)=4
LP(12)=3
MP(12)=1
CI(12)=1.
AE(12)=1.043D-18
AC(12)=6.59D-10
BC(12)=-1.5
ER(12)=-1.043D-18
!--REACTION 13 IS NO+NO>N+O+NO
LE(13)=5
ME(13)=5
KP(13)=4
LP(13)=3
MP(13)=5
CI(13)=1.
AE(13)=1.043D-18
AC(13)=1.318D-8
BC(13)=-1.5
ER(13)=-1.043D-18
!--REACTION 14 IS NO+O>N+O+O
LE(14)=5
ME(14)=3
KP(14)=4
LP(14)=3
MP(14)=3
CI(14)=1.
AE(14)=1.043D-18
AC(14)=1.318D-8
BC(14)=-1.5
ER(14)=-1.043D-18
!--REACTION 15 IS NO+N>2N+O
LE(15)=5
ME(15)=4
KP(15)=4
LP(15)=3
MP(15)=4
CI(15)=1.
AE(15)=1.043D-18
AC(15)=1.318D-8
BC(15)=-1.5
ER(15)=-1.043D-18
!--REACTION 16 IS NO+O>O2+N
LE(16)=5
ME(16)=3
KP(16)=0
LP(16)=1
MP(16)=4
CI(16)=0.
AE(16)=2.719D-19
AC(16)=5.279D-21
BC(16)=1.
ER(16)=-2.719D-19
!--REACTION 17 IS N2+O>NO+N
LE(17)=2
ME(17)=3
KP(17)=0
LP(17)=5
MP(17)=4
CI(17)=0.
AE(17)=5.175D-19
AC(17)=1.120D-16
BC(17)=0.
ER(17)=-5.175D-19
!--REACTION 18 IS O2+N>NO+O   !Exothermic Exchange
LE(18)=1
ME(18)=4
KP(18)=0
LP(18)=5
MP(18)=3
CI(18)=0.
AE(18)=4.968D-20
AC(18)=1.598D-18
BC(18)=0.5
ER(18)=2.719D-19
!--REACTION 19 IS NO+N>N2+O   !Exothermic Exchange
LE(19)=5
ME(19)=4
KP(19)=0
LP(19)=2
MP(19)=3
CI(19)=0.
AE(19)=0.
AC(19)=2.49D-17
BC(19)=0.
ER(19)=5.175D-19
!--REACTION 20 IS O+O+M1>O2+M1
LE(20)=3
ME(20)=3
KP(20)=-1
LP(20)=1
MP(20)=-1
CI(20)=0.
AE(20)=0.
AC(20)=8.297D-45
BC(20)=-0.5
ER(20)=8.197D-19
!--REACTION 21 IS N+N+M2>N2+M2
LE(21)=4
ME(21)=4
KP(21)=-1
LP(21)=2
MP(21)=-2
CI(21)=0.
AE(21)=0.
AC(21)=3.0051D-44
BC(21)=-0.5
ER(21)=1.561D-18
!--REACTION 22 IS N+N+N>N2+N
LE(22)=4
ME(22)=4
KP(22)=-1
LP(22)=2
MP(22)=-3
CI(22)=0.
AE(22)=0.
AC(22)=6.3962D-40
BC(22)=-1.5
ER(22)=1.5637D-18
!--REACTION 23 IS N+O+M3>NO+M3
LE(23)=4
ME(23)=3
KP(23)=-1
LP(23)=5
MP(23)=-4
CI(23)=0.
AE(23)=0.
AC(23)=2.7846D-40
BC(23)=-1.5
ER(23)=1.043D-18
!
  THBP(1,1)=9.
  THBP(1,2)=2.
  THBP(1,3)=25.
  THBP(1,4)=1.
  THBP(1,5)=1.
  THBP(2,1)=1.
  THBP(2,2)=2.5
  THBP(2,3)=1.
  THBP(2,4)=0.
  THBP(2,5)=1.
  THBP(3,1)=0.
  THBP(3,2)=0.
  THBP(3,3)=0.
  THBP(3,4)=1.
  THBP(3,5)=0.
  THBP(4,1)=1.
  THBP(4,2)=1.
  THBP(4,3)=20.
  THBP(4,4)=20.
  THBP(4,5)=20.
END IF
RETURN
END SUBROUTINE REAL_AIR
