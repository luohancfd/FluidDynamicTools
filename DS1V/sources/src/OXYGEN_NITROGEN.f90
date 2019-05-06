!
!***************************************************************************
!
SUBROUTINE OXYGEN_NITROGEN
!
USE GAS
USE CALC
USE MOLECS
USE EXPCOL,only : EXPCOL_INIT
!
IMPLICIT NONE
REAL(KIND=8) :: VDC=-1.d0, ZV=1.d12 !1.d12
INTEGER :: I,J,ZCHECK,IRATE=0,JRATE=0
CHARACTER (LEN=16) :: FILENAME
!--VDC>=0 for VDC model; <0 for tce model
!--IRATE=0 uncorrected rates; 1 corrected rates for recombination and exchange reactions
!--JRATE= same as IRATE but for dissociation reactions
!
MSP=8
MMRM=3
MMVM=3
!
IF (JCD == 1) THEN
  MNRE=0
  MTBP=0
  MEX=16
  MMEX=3
END IF
IF (JCD == 0) THEN
  IF(IRM == 0) MNRE=0     !Reactions are not included
  IF(IRM == 1) MNRE=37    !Shatalov
  IF(IRM == 2) MNRE=34    !Davidenko
  IF(IRM >= 100) MNRE=1   !testing only one reaction
  IF(IRM >= 200) MNRE=2   !testing only two reactions
  IF(IRM >= 400) MNRE=4   !testing only four reactions
  MTBP=0
  MEX=0
  MMEX=0
END IF
!
MNSR=0
!
CALL ALLOCATE_GAS
IQCTVT(1,4) = .true.
IQCTVT(4,1) = .true.
IQCTVT(3,4) = .true.
IQCTVT(4,3) = .true.
!
!--species 1 is hydrogen H2
!SP(1,1)=2.92D-10        !reference diameter  !vss/vhs  2.88/2.92
!SP(2,1)=273.D00         !reference temperature
!SP(3,1)=0.67D00         !viscosity index
!SP(4,1)=1.d00           !reciprocal of VSS scattering parameter; vss 1/1.35 (with BUG!)
!SP(5,1)=3.34D-27        !molecular mass
!ISPR(1,1)=2             !number of rotational degrees of freedom of species L
!ISPR(2,1)=0             !0,1 for constant, polynomial rotational relaxations collision number (Zrot)
!SPR(1,1)=5.d0           !=100 from Haas1995 !constant Zrot or constant in 2nd order polynomial in temperature
!ISPV(1)=1               !the number of vibrational modes
!SPVM(1,1,1)=6159.D00    !the characteristic vibrational temperature
!SPVM(2,1,1)=ZV  !18000.D00   !constant Zv, or reference Zv for mode K --estimate
!SPVM(3,1,1)=-1. !6159.D00    !-1 for constant Zv, or reference temperature
!SPVM(4,1,1)=52438.76D00 !charactheristic dissociation temperature; based on heats of formation
!ISPVM(1,1,1)=2          !the species code of the 1st dissociation product
!ISPVM(2,1,1)=2          !the species code of the 2nd dissociation product
!--species 2 is atomic hydrogen H
!SP(1,2)=2.33D-10         !SMILE/estimate 2.33/2.5
!SP(2,2)=273.D00
!SP(3,2)=0.75D00          !SMILE/estimate 0.75/0.80
!SP(4,2)=1.D00
!SP(5,2)=1.67D-27
!ISPR(1,2)=0
!ISPV(2)=0
!--species 1 is nitrogen N2
SP(1,1)=4.17d-10        !same as in DS1V
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0             !1 for polinomial Zr, 0 for constant Zr
SPR(1,1)=1.d0 !5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=ZV  !52000.D00 !260000.D00
SPVM(3,1,1)=-1.  !3371.D00
SPVM(4,1,1)=113500.D00
ISPVM(1,1,1)=2
ISPVM(2,1,1)=2
!--species 2 is atomic nitrogen N
SP(1,2)=3.d-10          !same as inDS1V
SP(2,2)=273.D00
SP(3,2)=0.8D00
SP(4,2)=1.0D00
SP(5,2)=2.325D-26
ISPR(1,2)=0
ISPV(2)=0
!--species 3 is oxygen O2
!SP(1,3)=3.985d-10       !wysong2014/vss/vhs 3.985/4.01/4.07
!SP(2,3)=273.D00         !t_ref
!SP(3,3)=0.71D00         !wysong2014

SP(1,3)=4.1515d-10     !Han correct by fitting esposito's data
SP(2,3)=273.D0
SP(3,3)=0.7318d0

SP(4,3)=1.d0            !vss 1./1.4d0 !there is some bug with vss model
SP(5,3)=5.312D-26       !vss correction in collision subroutine must be tested
ISPR(1,3)=2
ISPR(2,3)=0             !1 for polinomial Zr, 0 for constant Zr
SPR(1,3)=1.0 !5.d0           !constant Zr
ISPV(3)=1               !the number of vibrational modes
SPVM(1,1,3)=2273.54     !obtained from nist 2256.D00    !the characteristic vibrational temperature
SPVM(2,1,3)=ZV  !18000.D00   !constant Zv, or the reference Zv
SPVM(3,1,3)=-1. !2256.d0     !-1 for a constant Zv, or the reference temperature
SPVM(4,1,3)=59971.4D00  !the characteristic dissociation temperature
ISPVM(1,1,3)=4
ISPVM(2,1,3)=4
!--species 4 is atomic oxygen O
SP(1,4)=3.458d-10        !wysong2014
SP(2,4)=273.D00          !t_ref
SP(3,4)=0.76D00          !wysong2014
SP(4,4)=1.D00
SP(5,4)=2.656D-26
ISPR(1,4)=0
SPR(1,4)=5.d0
ISPV(4)=0
!--species 5 is hydroxy OH
!SP(1,5)=3.50d-10          !--SMILE/estimate 3.50/4.10
!SP(2,5)=273.D00
!SP(3,5)=0.75D00         !-estimate
!SP(4,5)=1.0D00
!SP(5,5)=2.823D-26
!ISPR(1,5)=2
!ISPR(2,5)=0
!SPR(1,5)=5.D00
!ISPV(5)=1
!SPVM(1,1,5)=5360.D00
!SPVM(2,1,5)=ZV  !18000.D00   !--estimate
!SPVM(3,1,5)=-1. !5360.D00    !--estimate
!SPVM(4,1,5)=51497.18D00
!ISPVM(1,1,5)=2
!ISPVM(2,1,5)=4
!--species 5 is nitric oxide NO
SP(1,5)=4.2D-10        !same as in DS1V
SP(2,5)=273.D00
SP(3,5)=0.79D00
SP(4,5)=1.0D00
SP(5,5)=4.981D-26
ISPR(1,5)=2
ISPR(2,5)=0
SPR(1,5)=5.D00
ISPV(5)=1
SPVM(1,1,5)=2719.D00
SPVM(2,1,5)=ZV  !14000.D00   !70000.D00
SPVM(3,1,5)=-1. !2719.D00
SPVM(4,1,5)=75500.D00
ISPVM(1,1,5)=2
ISPVM(2,1,5)=4
!--species 6 is water vapor H2O
SP(1,6)=4.5d-10         !--estimate
SP(2,6)=273.D00
SP(3,6)=1.0D00         !SMILE/estimate 1.0/0.75
SP(4,6)=1.0D00
SP(5,6)=2.99D-26
ISPR(1,6)=3
ISPR(2,6)=0
SPR(1,6)=10.d0        !--estimate Alexeenko (2003)
ISPV(6)=3
SPVM(1,1,6)=5261.D00           !--symmetric stretch mode
SPVM(2,1,6)=10.d0 !ZV !10.d0   !--estimate Alexeenko (2003); for same species collisions
SPVM(3,1,6)=-2.                !-2 indicates that Zv is different for same species collisions
SPVM(4,1,6)=60043.83D00
SPVM(1,2,6)=2294.D00  !--bend mode
SPVM(2,2,6)=250.d0 !ZV !250.d0    !--estimate Alexeenko (2003); for different species collisions
SPVM(3,2,6)=SPVM(3,1,6)
SPVM(4,2,6)=60043.83D00
SPVM(1,3,6)=5432.D00   !--asymmetric stretch mode
SPVM(2,3,6)=SPVM(2,2,6) !ZV !250.d0
SPVM(3,3,6)=SPVM(3,2,6)
SPVM(4,3,6)=60043.83D00
ISPVM(1,1,6)=2
ISPVM(2,1,6)=5
ISPVM(1,2,6)=2
ISPVM(2,2,6)=5
ISPVM(1,3,6)=2
ISPVM(2,3,6)=5
!--species 7 is hydroperoxy HO2
SP(1,7)=5.5d-10       !--estimate
SP(2,7)=273.D00
SP(3,7)=0.75D00      !-estimate
SP(4,7)=1.0D00
SP(5,7)=5.479D-26
ISPR(1,7)=2    !--assumes that HO2 is linear
ISPR(2,7)=0
SPR(1,7)=10.D00   !estimate
ISPV(7)=3
SPVM(1,1,7)=4950.D00
SPVM(2,1,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,1,7)=-1. ! 2500.D00    !--estimate
SPVM(4,1,7)=24988.08D00
SPVM(1,2,7)=2000.D00
SPVM(2,2,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,2,7)=-1. !2500.D00    !--estimate
SPVM(4,2,7)=24988.08D00
SPVM(1,3,7)=1580.D00
SPVM(2,3,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,3,7)=-1. !2500.D00    !--estimate
SPVM(4,3,7)=24988.08D00
ISPVM(1,1,7)=2
ISPVM(2,1,7)=3
ISPVM(1,2,7)=2
ISPVM(2,2,7)=3
ISPVM(1,3,7)=2
ISPVM(2,3,7)=3
!--species 8 is hydroperoxy CO2
SP(1,8)=5.687766d-10       !SMILE
SP(2,8)=273.d0
SP(3,8)=0.93d0
SP(4,8)=1.0d0
SP(5,8)=7.309999d-26
ISPR(1,8)=2    !--assumes that CO2 is linear
ISPR(2,8)=0
SPR(1,8)=10.d0 !estimate
ISPV(8)=3
SPVM(1,1,8)=945.d0
SPVM(2,1,8)=ZV !250.d0    !--estimate; for different species collisions
SPVM(3,1,8)=-1.       !-2 indicates that Zv is different for same species collisions
SPVM(4,1,8)=64015.d0
SPVM(1,2,8)=1903.d0
SPVM(2,2,8)=ZV !10.d0     !--for same species collisions
SPVM(3,2,8)=-1.
SPVM(4,2,8)=64015.d0
SPVM(1,3,8)=3329.d0
SPVM(2,3,8)=ZV !10.d0
SPVM(3,3,8)=-1.
SPVM(4,3,8)=64015.d0
ISPVM(1,1,8)=2
ISPVM(2,1,8)=3
ISPVM(1,2,8)=2
ISPVM(2,2,8)=3
ISPVM(1,3,8)=2
ISPVM(2,3,8)=3
!--Species 8 is argon
!SP(1,8)=4.17D-10  !vss/vhs 4.11/4.17
!SP(2,8)=273.15
!SP(3,8)=0.81
!SP(4,8)=1.        !vss 1./1.4
!SP(5,8)=6.63D-26
!ISPR(1,8)=0
!ISPV(8)=0
!
!ISPV(:)=0 !test mechanism without vibrational modes
!
IF (JCD ==1) THEN    !--Q-K model
!all the lines were deletd --isebasti
END IF
!
!--the following data is required if JCD=0 (TCE reaction model)
!
IF (JCD == 0 .AND. IRM == 2) THEN
!--rate data is based on Davidenko2006_SystematicStudySupersonicCombustion
!--------------------------------------------------------------
!--REACTION 1 IS H+H+H2>2H2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
  IREV(J)=7
!--REACTION 2 IS H+H+H>H2+H
  J=2
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
  IREV(J)=8
!--REACTION 3 IS H+H+O2>H2+O2
  J=3
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
  IREV(J)=9
!--REACTION 4 IS H+H+O>H2+O
  J=4
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
  IREV(J)=10
!--REACTION 5 IS H+H+OH>H2+OH
  J=5
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
  IREV(J)=11
!--REACTION 6 IS H+H+H2O>H2+H2O
  J=6
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
  IREV(J)=12
!--------------------------------------------------------------
! REACTION 7 IS H2+H2->2H+H2
  J=7
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
!--REACTION 8 IS H2+H>2H+H
  J=8
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=2
!--REACTION 9 IS H2+O2>2H+O2
  J=9
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=3
!--REACTION 10 IS H2+O>2H+O
  J=10
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=4
!--REACTION 11 IS H2+OH>2H+OH
  J=11
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=5
!--REACTION 12 IS H2+H2O>2H+H2O
  J=12
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=6
!--------------------------------------------------------------
!--REACTION 13 IS OH+H+H2>H2O+H2
  J=13
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
  IREV(J)=19
!--REACTION 14 IS OH+H+H>H2O+H
  J=14
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
  IREV(J)=20
!--REACTION 15 IS OH+H+O2>H2O+O2
  J=15
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
  IREV(J)=21
!--REACTION 16 IS OH+H+O>H2O+O
  J=16
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
  IREV(J)=22
!--REACTION 17 IS OH+H+OH>H2O+OH
  J=17
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
  IREV(J)=23
!--REACTION 18 IS OH+H+H2O>H2O+H2O
  J=18
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
  IREV(J)=24
!--------------------------------------------------------------
!--REACTION 19 IS H2O+H2>OH+H+H2
  J=19
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=13
!--REACTION 20 IS H2O+H>OH+H+H
  J=20
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=14
!--REACTION 21 IS H2O+O2>OH+H+O2
  J=21
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=15
!--REACTION 22 IS H2O+O>OH+H+O
  J=22
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=16
!--REACTION 23 IS H2O+OH>OH+H+OH
  J=23
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=17
!--REACTION 24 IS H2O+H2O>OH+H+H2O
  J=24
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=18
!--------------------------------------------------------------
!--REACTION 25 IS OH+OH>H2O+O
  J=25
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
  IREV(J)=26
!--REACTION 26 IS H2O+O>OH+OH
  J=26
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
  IREV(J)=25
!--------------------------------------------------------------
!--REACTION 27 IS H2+O>OH+H
  J=27
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
  IREV(J)=28
!--REACTION 28 IS OH+H>H2+O
  J=28
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
  IREV(J)=27
!--------------------------------------------------------------
!--REACTION 29 IS H2+OH>H20+H
  J=29
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
  IREV(J)=30
!--REACTION 30 IS H2O+H>H2+OH
  J=30
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ  !Dav 9300.d0  !Shat 9030.d0  !Conaire 9252.d0
  AC(J)=7.964d8/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  BC(J)=1.528d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  IF (IRATE == 1) BC(J)=1.939508d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
  IREV(J)=29
!--------------------------------------------------------------
!--REACTION 31 IS O2+H>OH+O
  J=31
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
  IREV(J)=32
!--REACTION 32 IS OH+O>O2+H
  J=32
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
  IREV(J)=31
!--------------------------------------------------------------
!--REACTION 33 IS H2+O2>OH+OH
  J=33
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
  IREV(J)=34
!--REACTION 34 IS OH+OH>H2+O2
  J=34
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
  IREV(J)=33
END IF
!
!--------------------------------------------------------------
!Testing individual reactions
J=1
!
IF (JCD == 0 .AND. IRM == 101) THEN
!--REACTION D1 IS H+H+H2>2H2  (Davidenko)
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 102) THEN
!--REACTION D2 IS H+H+H>H2+H
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 103) THEN
!--REACTION D3 IS H+H+O2>H2+O2
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 104) THEN
!--REACTION D4 IS H+H+O>H2+O
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 105) THEN
!--REACTION D5 IS H+H+OH>H2+OH
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 106) THEN
!--REACTION D6 IS H+H+H2O>H2+H2O
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 107) THEN
! REACTION D7 IS H2+H2->2H+H2
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 108) THEN
!--REACTION D8 IS H2+H>2H+H
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 109) THEN
!--REACTION D9 IS H2+O2>2H+O2
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 110) THEN
!--REACTION D10 IS H2+O>2H+O
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 111) THEN
!--REACTION D11 IS H2+OH>2H+OH
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 112) THEN
!--REACTION D12 IS H2+H2O>2H+H2O
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 113) THEN
!--REACTION D13 IS OH+H+H2>H2O+H2
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 114) THEN
!--REACTION D14 IS OH+H+H>H2O+H
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 115) THEN
!--REACTION D15 IS OH+H+O2>H2O+O2
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 116) THEN
!--REACTION D16 IS OH+H+O>H2O+O
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 117) THEN
!--REACTION D17 IS OH+H+OH>H2O+OH
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 118) THEN
!--REACTION D18 IS OH+H+H2O>H2O+H2O
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 119) THEN
!--REACTION D19 IS H2O+H2>OH+H+H2
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 120) THEN
!--REACTION D20 IS H2O+H>OH+H+H
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 121) THEN
!--REACTION D21 IS H2O+O2>OH+H+O2
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 122) THEN
!--REACTION D22 IS H2O+O>OH+H+O
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 123) THEN
!--REACTION D23 IS H2O+OH>OH+H+OH
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 124) THEN
!--REACTION D24 IS H2O+H2O>OH+H+H2O
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 125) THEN
!--REACTION D25 IS OH+OH>H20+O
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
END IF
!
IF (JCD == 0 .AND. IRM == 126) THEN
!--REACTION D26 IS H2O+O>OH+OH
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
END IF
!
IF (JCD == 0 .AND. IRM == 127) THEN
!--REACTION D27 IS H2+O>OH+H
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
END IF
!
IF (JCD == 0 .AND. IRM == 128) THEN
!--REACTION D28 IS OH+H>H2+O
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
END IF
!
IF (JCD == 0 .AND. IRM == 129) THEN
!--REACTION D29 IS H2+OH>H20+H
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
END IF
!
IF (JCD == 0 .AND. IRM == 130) THEN
!--REACTION D30 IS H2O+H>H2+OH
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ
  AC(J)=7.964d8/(AVOG*1000.d0)
  BC(J)=1.528d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.939508d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
END IF
!
IF (JCD == 0 .AND. IRM == 131) THEN
!--REACTION D31 IS O2+H>OH+O
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
END IF
!
IF (JCD == 0 .AND. IRM == 132) THEN
!--REACTION D32 IS OH+O>O2+H
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
END IF
!
IF (JCD == 0 .AND. IRM == 133) THEN
!--REACTION D33 IS H2+O2>OH+OH
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
END IF
!
IF (JCD == 0 .AND. IRM == 134) THEN
!--REACTION D34 IS OH+OH>H2+O2
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
END IF
!
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 140) THEN
!--REACTION 140 IS N2+N>N+N+N (from QCT data)
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=113200.d0*BOLTZ         !from AHO model
  AC(J)=6.6831d-6*1.0d-6        ! fit from Jaffe data
  BC(J)=-0.6996d0
  ER(J)=-AE(J)
END IF
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 141) THEN
!--REACTION 140 IS N2+N2>N2+N+N (from Bender data)
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=117000.d0*BOLTZ         !from AHO model
  AC(J)=4.5d-6*1.0d-6        ! fit from Jaffe data
  BC(J)=-0.675
  ER(J)=-AE(J)
END IF
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 150) THEN
!--REACTION 150 IS O2+O>O+O+O (from QCT data)
  LE(J)=3
  ME(J)=4
  KP(J)=4
  LP(J)=4
  MP(J)=4
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT          !from AHO model
  AC(J)=2.500d18/(AVOG*1000.d0)  !AC and BC are fit from QCT data
  BC(J)=-0.565d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN
      AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-0.185509d0
    ELSE
      AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-0.185509d0
    END IF
  END IF
  ER(J)=-5.21275d0*EVOLT !based on heat of formation !-5.21275d0*EVOLT !from AHO model
END IF
!
IF (JCD == 0 .AND. IRM == 160) THEN
!--REACTION 160 IS O2+O2>O+O+O2 (fit from park1990 but using the DSMC value of activation energy)
  LE(J)=3
  ME(J)=3
  KP(J)=4
  LP(J)=4
  MP(J)=3
  CI(J)=VDC
  !AE(J)=5.21275d0*EVOLT           !from AHO model
  !AC(J)=1.9337d22/(AVOG*1000.d0)  !AC and BC from park1990
  !BC(J)=-1.7334d0
  AE(J) = 59380.d0*BOLTZ+(5.21275d0-5.1153d0)*EVOLT
   !Han use Byron's rate but shift to add zeropoint energy
  AC(J) = 1.9392865d21/(AVOG*1000.d0)
  BC(J) = -1.5d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN
      AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-1.1919d0
    ELSE
      AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-1.1919d0
    END IF
  END IF
  !ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 161) THEN
!--REACTION 160 IS O2+N2>O+O+N2 (fit from park1990 but using the DSMC value of activation energy)
  LE(J)=3
  ME(J)=1
  KP(J)=4
  LP(J)=4
  MP(J)=1
  CI(J)=VDC
  AE(J) = 5.21275d0*EVOLT
  AC(J) = 1.586878535d15/(AVOG*1000.d0)
  BC(J) = 4.255d-3
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 162) THEN
!--REACTION 160 IS N2+O2>N+N+O2
  LE(J) = 1
  ME(J) = 3
  KP(J) = 2
  LP(J) = 2
  MP(J) = 3
  CI(J) = VDC
  !  Wray: doi:10.2514/5.2018-0240 fit good with Andrienko's QCT
  ! here we refit to match dissociation energy
  AC(J) = 1.1447510d18/(AVOG*1000.d0)
  BC(J) = -0.6823d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 163) THEN
!--REACTION 160 IS NO+N2>N+O+N2
  LE(J) = 5
  ME(J) = 1
  KP(J) = 2
  LP(J) = 4
  MP(J) = 1
  CI(J) = VDC
  ! Park's rate, match well with Andrienko's QCT
  ! 5.0e15, 0 , 75500K
  ! here we refit to match dissociation energy
  AC(J) = 1.08705682d16/(AVOG*1000.d0)
  BC(J) = -0.0766715d0
  AE(J) = 6.55879*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 164) THEN
!--REACTION 160 IS N2+NO>N+N+NO
  LE(J) = 1
  ME(J) = 5
  KP(J) = 2
  LP(J) = 2
  MP(J) = 5
  CI(J) = VDC
  ! Park uses N2+N2 as N2+NO
  ! Wray suggest N2+Ar as N2+NO
  ! here we refit Wray's data to match dissociation energy
  AC(J) = 9.0558360d17/(AVOG*1000.d0)
  BC(J) = -0.662299d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 170) THEN
!--REACTION 170 IS N2+O>2N+O (from QCT data)
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
 ! AE(J)=9.8216d0*EVOLT            !from AHO model
 ! AC(J)=8.9330d16/(AVOG*1000.d0)  !AC and BC are fit from QCT data
 ! BC(J)=-0.3842d0
  !IF (JRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)  !to be corrected
  !IF (JRATE == 1) BC(J)=0.7079989d0                 !to be corrected
 ! ER(J)=-9.8216d0*EVOLT !from AHO model
  AE(J) = 113950.d0*BOLTZ
  AC(J) = 8.934d-13
  BC(J) = -0.4807
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 180) THEN
!--REACTION 180 IS N2+O>NO+N (from QCT data)
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=VDC
  AE(J)=3.2628d0*EVOLT           !from AHO model
  AC(J)=2.4444d11/(AVOG*1000.d0) !AC and BC are fit from QCT data
  BC(J)=0.7071d0
  !IF (JRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00) !to be corrected
  !IF (JRATE == 1) BC(J)=0.7079989d0                !to be corrected
  ER(J)=-3.2628d0*EVOLT !from AHO model
END IF
!--------------------------------------------------------------
!Testing reaction pairs
!
IF (JCD == 0 .AND. IRM == 201) THEN
!--REACTION 1 IS H+H+H2>2H2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
  IREV(J)=2
! REACTION 7 IS H2+H2->2H+H2
  J=2
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 202) THEN
!--REACTION 2 IS H+H+H>H2+H
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 8 IS H2+H>2H+H
  J=2
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 203) THEN
!--REACTION 3 IS H+H+O2>H2+O2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 9 IS H2+O2>2H+O2
  J=2
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 204) THEN
!--REACTION 4 IS H+H+O>H2+O
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 10 IS H2+O>2H+O
  J=2
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 205) THEN
!--REACTION 5 IS H+H+OH>H2+OH
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 11 IS H2+OH>2H+OH
  J=2
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 206) THEN
!--REACTION 6 IS H+H+H2O>H2+H2O
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 12 IS H2+H2O>2H+H2O
  J=2
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 213) THEN
!--REACTION 13 IS OH+H+H2>H2O+H2
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 19 IS H2O+H2>OH+H+H2
  J=2
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 214) THEN
!--REACTION 14 IS OH+H+H>H2O+H
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 20 IS H2O+H>OH+H+H
  J=2
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 215) THEN
!--REACTION 15 IS OH+H+O2>H2O+O2
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 21 IS H2O+O2>OH+H+O2
  J=2
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 216) THEN
!--REACTION 16 IS OH+H+O>H2O+O
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 22 IS H2O+O>OH+H+O
  J=2
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 217) THEN
!--REACTION 17 IS OH+H+OH>H2O+OH
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 23 IS H2O+OH>OH+H+OH
  J=2
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 218) THEN
!--REACTION 18 IS OH+H+H2O>H2O+H2O
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 24 IS H2O+H2O>OH+H+H2O
  J=2
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 225) THEN
!--REACTION 25 IS OH+OH>H20+O
  J=1
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
  IREV(J)=2
!--REACTION 26 IS H2O+O>OH+OH
  J=2
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 227) THEN
!--REACTION 27 IS H2+O>OH+H
  J=1
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
  IREV(J)=2
!--REACTION 28 IS OH+H>H2+O
  J=2
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 229) THEN
!--REACTION 29 IS H2+OH>H20+H
  J=1
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
  IREV(J)=2
!--REACTION 30 IS H2O+H>H2+OH
  J=2
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ  !Dav 9300.d0  !Shat 9030.d0  !Conaire 9252.d0
  AC(J)=7.964d8/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  BC(J)=1.528d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  IF (IRATE == 1) BC(J)=1.939508d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 231) THEN
!--REACTION 31 IS O2+H>OH+O
  J=1
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
  IREV(J)=2
!--REACTION 32 IS OH+O>O2+H
  J=2
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 233) THEN
!--REACTION 33 IS H2+O2>OH+OH
  J=1
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
  IREV(J)=2
!--REACTION 34 IS OH+OH>H2+O2
  J=2
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 250) THEN
!--REACTION 150 IS O2+O>O+O+O (from QCT data)
  J=1
  LE(J)=3
  ME(J)=4
  KP(J)=4
  LP(J)=4
  MP(J)=4
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT          !from AHO model
  AC(J)=2.500d18/(AVOG*1000.d0)  !AC and BC are fit from QCT data
  BC(J)=-0.565d0
  IF (IRATE == 1) AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
  IF (IRATE == 1) BC(J)=-0.185509d0
  ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  IREV(J)=2
!--REACTION 160 IS O2+O2>O+O+O2 (fit from park1990 but using the DSMC value of activation energy)
  J=2
  LE(J)=3
  ME(J)=3
  KP(J)=4
  LP(J)=4
  MP(J)=3
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT           !from AHO model
  AC(J)=1.9337d22/(AVOG*1000.d0)  !AC and BC from park1990
  BC(J)=-1.7334d0
  IF (IRATE == 1) AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
  IF (IRATE == 1) BC(J)=-1.1919d0
  ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 251) THEN
!--REACTION 1 is O2+N2>O+O+N2
  LE(1) = 3
  ME(1) = 1
  KP(1) = 4
  LP(1) = 4
  MP(1) = 1
  CI(1) = VDC
  ! Park 2E21cc/mol/sec -1.5 59700
  ! Boyd 8.132E-9 cc/sec -0.131 59380
  ! fit1 2.7212E14 cc/mol/sec 0.22006 59380
  ! fit2 1.5869E15 cc/mol/sec 0.004255  60490
  AC(1) = 1.586878535d15/(AVOG*1000.d0)
  BC(1) = 4.255d-3
  AE(1) = 5.21275d0*EVOLT
  ER(1) = -AE(1)
!--REACTION 2 is N2+O2->N+N+O2
  LE(2) = 1
  ME(2) = 3
  KP(2) = 2
  LP(2) = 2
  MP(2) = 3
  CI(2) = VDC
  !  Wray: doi:10.2514/5.2018-0240
  ! here we refit to match dissociation energy
  AC(2) = 1.1447510d18/(AVOG*1000.d0)
  BC(2) = -0.6823
  AE(2) = 9.82163d0*EVOLT
  ER(2) = -AE(2)
END IF
!
IF (JCD == 0 .AND. IRM == 252) THEN
  !--REACTION 160 IS NO+N2>N+O+N2
  J = 1
  LE(J) = 5
  ME(J) = 1
  KP(J) = 2
  LP(J) = 4
  MP(J) = 1
  CI(J) = VDC
  ! Park's rate, match well with Andrienko's QCT
  ! 5.0e15, 0 , 75500K
  ! here we refit to match dissociation energy
  AC(J) = 1.08705682d16/(AVOG*1000.d0)
  BC(J) = -0.0766715d0
  AE(J) = 6.55879*EVOLT
  ER(J) = -AE(J)
  !--REACTION 160 IS N2+NO>N+N+NO
  J = 2
  LE(J) = 1
  ME(J) = 5
  KP(J) = 2
  LP(J) = 2
  MP(J) = 5
  CI(J) = VDC
  ! Park uses N2+N2 as N2+NO
  ! Wray suggest N2+Ar as N2+NO
  ! here we refit Wray's data to match dissociation energy
  AC(J) = 9.0558360d17/(AVOG*1000.d0)
  BC(J) = -0.662299d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!--------------------------------------------------------------
!
IF (IPRS == 1) THEN !--read pre-collision vibrational pdfs
  DO I=1,MNRE
    DO J=1,ITMAX
      WRITE(FILENAME,700) I,J
      700 FORMAT('DS1VIBF1_',i2.2,'_',i2.2)
      800 CONTINUE
      OPEN (7,FILE=TRIM(FILENAME)//'.BIN',FORM='UNFORMATTED',ERR=800)
      READ (7) FPTEMP(J),FEVIB(I,J,:,:),FPVIB(I,J,:,:,:),ZCHECK
      CLOSE(7)
      IF (ZCHECK /= 1234567) THEN
        WRITE (9,*) 'Failed reading DS1FVIB for IKA-IT-TEMP; Check:',I,J,FPTEMP(J),ZCHECK
        STOP
      ELSE
        WRITE (9,*) 'Success readding DS1FVIB for IKA-IT-TEMP; Check:',I,J,FPTEMP(J),ZCHECK
      END IF
    END DO
  END DO
END IF
!--------------------------------------------------------------
IF (QCTMODEL > 0) THEN
  VIBEN=1.d3            !initializing variable
!
  J=1                   !N2 vibrational levels from: Han's Luo (Purdue AAE) Master Thesis (2016)
  IVMODEL(J,1)=1        !use pre-computed vibrational levels for this species
  IVMODEL(J,2)=59        !max possible vibrational level
  VIBEN(0,J)=-9.758052d0
  VIBEN(1,J)=-9.467844d0
  VIBEN(2,J)=-9.181767d0
  VIBEN(3,J)=-8.899808d0
  VIBEN(4,J)=-8.621953d0
  VIBEN(5,J)=-8.348193d0
  VIBEN(6,J)=-8.078518d0
  VIBEN(7,J)=-7.812921d0
  VIBEN(8,J)=-7.551397d0
  VIBEN(9,J)=-7.293940d0
  VIBEN(10,J)=-7.040548d0
  VIBEN(11,J)=-6.791219d0
  VIBEN(12,J)=-6.545954d0
  VIBEN(13,J)=-6.304754d0
  VIBEN(14,J)=-6.067621d0
  VIBEN(15,J)=-5.834561d0
  VIBEN(16,J)=-5.605578d0
  VIBEN(17,J)=-5.380680d0
  VIBEN(18,J)=-5.159875d0
  VIBEN(19,J)=-4.943173d0
  VIBEN(20,J)=-4.730586d0
  VIBEN(21,J)=-4.522125d0
  VIBEN(22,J)=-4.317806d0
  VIBEN(23,J)=-4.117644d0
  VIBEN(24,J)=-3.921656d0
  VIBEN(25,J)=-3.729860d0
  VIBEN(26,J)=-3.542278d0
  VIBEN(27,J)=-3.358931d0
  VIBEN(28,J)=-3.179843d0
  VIBEN(29,J)=-3.005039d0
  VIBEN(30,J)=-2.834546d0
  VIBEN(31,J)=-2.668394d0
  VIBEN(32,J)=-2.506615d0
  VIBEN(33,J)=-2.349241d0
  VIBEN(34,J)=-2.196308d0
  VIBEN(35,J)=-2.047855d0
  VIBEN(36,J)=-1.903922d0
  VIBEN(37,J)=-1.764553d0
  VIBEN(38,J)=-1.629794d0
  VIBEN(39,J)=-1.499695d0
  VIBEN(40,J)=-1.374309d0
  VIBEN(41,J)=-1.253693d0
  VIBEN(42,J)=-1.137909d0
  VIBEN(43,J)=-1.027022d0
  VIBEN(44,J)=-9.211040d-1
  VIBEN(45,J)=-8.202310d-1
  VIBEN(46,J)=-7.244866d-1
  VIBEN(47,J)=-6.339614d-1
  VIBEN(48,J)=-5.487544d-1
  VIBEN(49,J)=-4.689744d-1
  VIBEN(50,J)=-3.947413d-1
  VIBEN(51,J)=-3.261886d-1
  VIBEN(52,J)=-2.634655d-1
  VIBEN(53,J)=-2.067411d-1
  VIBEN(54,J)=-1.562088d-1
  VIBEN(55,J)=-1.120932d-1
  VIBEN(56,J)=-7.466090d-2
  VIBEN(57,J)=-4.423639d-2
  VIBEN(58,J)=-2.123195d-2
  VIBEN(59,J)=-6.208831d-3
  VIBEN(:,J)=9.821630d0+VIBEN(:,J)  !add the zero point energy
!
  J=3                   !O2 vibrational levels from: Esposito et al, Chemical Physics 351 (2008)
  IVMODEL(J,1)=1        !use pre-computed vibrational levels for this species
  IVMODEL(J,2)=46       !max possible vibrational level
  VIBEN(0,J)=-5.1153d0  !in eV
  VIBEN(1,J)=-4.9221d0
  VIBEN(2,J)=-4.7315d0
  VIBEN(3,J)=-4.5433d0
  VIBEN(4,J)=-4.3578d0
  VIBEN(5,J)=-4.1749d0
  VIBEN(6,J)=-3.9948d0
  VIBEN(7,J)=-3.8175d0
  VIBEN(8,J)=-3.6431d0
  VIBEN(9,J)=-3.4717d0
  VIBEN(10,J)=-3.3034d0
  VIBEN(11,J)=-3.1381d0
  VIBEN(12,J)=-2.9760d0
  VIBEN(13,J)=-2.8172d0
  VIBEN(14,J)=-2.6617d0
  VIBEN(15,J)=-2.5096d0
  VIBEN(16,J)=-2.3610d0
  VIBEN(17,J)=-2.2160d0
  VIBEN(18,J)=-2.0746d0
  VIBEN(19,J)=-1.9369d0
  VIBEN(20,J)=-1.8030d0
  VIBEN(21,J)=-1.6729d0
  VIBEN(22,J)=-1.5469d0
  VIBEN(23,J)=-1.4249d0
  VIBEN(24,J)=-1.3071d0
  VIBEN(25,J)=-1.1935d0
  VIBEN(26,J)=-1.0842d0
  VIBEN(27,J)=-0.97939d0
  VIBEN(28,J)=-0.87913d0
  VIBEN(29,J)=-0.78354d0
  VIBEN(30,J)=-0.69274d0
  VIBEN(31,J)=-0.60685d0
  VIBEN(32,J)=-0.52601d0
  VIBEN(33,J)=-0.45036d0
  VIBEN(34,J)=-0.38004d0
  VIBEN(35,J)=-0.31522d0
  VIBEN(36,J)=-0.25607d0
  VIBEN(37,J)=-0.20277d0
  VIBEN(38,J)=-0.15550d0
  VIBEN(39,J)=-0.11445d0
  VIBEN(40,J)=-0.079818d0
  VIBEN(41,J)=-0.051751d0
  VIBEN(42,J)=-0.030323d0
  VIBEN(43,J)=-0.015397d0
  VIBEN(44,J)=-0.0063748d0
  VIBEN(45,J)=-0.0019261d0
  VIBEN(46,J)=-0.00029275d0
  VIBEN(:,J)=5.21275d0+VIBEN(:,J)  !add the zero point energy
!
  J=5                    !NO vibrational levels from: Han's Luo (Purdue AAE) Master Thesis (2016)
  IVMODEL(J,2)=52        !max possible vibrational level
  IVMODEL(J,1)=1         !use pre-computed vibrational levels for this species
  VIBEN(0,J)=-6.495503d0 !in eV
  VIBEN(1,J)=-6.261055d0
  VIBEN(2,J)=-6.030857d0
  VIBEN(3,J)=-5.804906d0
  VIBEN(4,J)=-5.583201d0
  VIBEN(5,J)=-5.365739d0
  VIBEN(6,J)=-5.152521d0
  VIBEN(7,J)=-4.943543d0
  VIBEN(8,J)=-4.738807d0
  VIBEN(9,J)=-4.538309d0
  VIBEN(10,J)=-4.342050d0
  VIBEN(11,J)=-4.150030d0
  VIBEN(12,J)=-3.962248d0
  VIBEN(13,J)=-3.778704d0
  VIBEN(14,J)=-3.599398d0
  VIBEN(15,J)=-3.424332d0
  VIBEN(16,J)=-3.253505d0
  VIBEN(17,J)=-3.086920d0
  VIBEN(18,J)=-2.924577d0
  VIBEN(19,J)=-2.766480d0
  VIBEN(20,J)=-2.612630d0
  VIBEN(21,J)=-2.463031d0
  VIBEN(22,J)=-2.317685d0
  VIBEN(23,J)=-2.176598d0
  VIBEN(24,J)=-2.039773d0
  VIBEN(25,J)=-1.907217d0
  VIBEN(26,J)=-1.778935d0
  VIBEN(27,J)=-1.654934d0
  VIBEN(28,J)=-1.535221d0
  VIBEN(29,J)=-1.419806d0
  VIBEN(30,J)=-1.308699d0
  VIBEN(31,J)=-1.201909d0
  VIBEN(32,J)=-1.099450d0
  VIBEN(33,J)=-1.001335d0
  VIBEN(34,J)=-9.075788d-1
  VIBEN(35,J)=-8.181986d-1
  VIBEN(36,J)=-7.332134d-1
  VIBEN(37,J)=-6.526442d-1
  VIBEN(38,J)=-5.765146d-1
  VIBEN(39,J)=-5.048512d-1
  VIBEN(40,J)=-4.376836d-1
  VIBEN(41,J)=-3.750456d-1
  VIBEN(42,J)=-3.169752d-1
  VIBEN(43,J)=-2.635155d-1
  VIBEN(44,J)=-2.147164d-1
  VIBEN(45,J)=-1.706349d-1
  VIBEN(46,J)=-1.313377d-1
  VIBEN(47,J)=-9.690349d-2
  VIBEN(48,J)=-6.742670d-2
  VIBEN(49,J)=-4.302306d-2
  VIBEN(50,J)=-2.383927d-2
  VIBEN(51,J)=-1.007100d-2
  VIBEN(52,J)=-2.005209d-3
  VIBEN(:,J)=6.55879d0+VIBEN(:,J)   !add the NO dissociation energy
!
  VIBEN=VIBEN*EVOLT  !convert to Joules
END IF
!--------------------------------------------------------------
!
RETURN
END SUBROUTINE OXYGEN_NITROGEN
