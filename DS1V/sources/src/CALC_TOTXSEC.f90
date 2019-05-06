!
!*****************************************************************************
!
SUBROUTINE CALC_TOTXSEC(LS,MS,VR,VRR,ET0,EVIB,TOTXSEC,BMAX,CVR)
!
!--calculate the total cross sections, i.e. collision cross sections
! EVIB is the vibrational energy of N2 molecule for N2-O nonVHS model
! TOTXSEC in angstrom^2
! CVR,VR,VRR, EVIB in si
! ET0 in eV
USE GAS, only : INONVHS, SPM
USE CALC, only: BOLTZ, EVOLT, PI
USE EXPCOL, only: EXPCOL_TOTXSEC
IMPLICIT NONE
INTEGER,intent(in) :: LS,MS
REAL(8),intent(in) :: VRR, EVIB, VR, ET0
REAL(8) :: EVIBEV
REAL(8),intent(out) ::TOTXSEC, CVR, BMAX
REAL(8),PARAMETER :: CTOT(8) = (/46.039504155886434, 0.051004420857885, -0.584551057667247, &
                              &  -0.002969283806997, 0.281618756125794, 0.030181202283512, &
                              &  -0.436592532266083, 0.152224780739684/)
INTEGER :: IERROR

! ET = 0.5d0*SPM(1,LS,MS)*VRR/EVOLT  ! collisional energy in eV

IF (INONVHS(LS, MS) == 2) THEN
  CALL EXPCOL_TOTXSEC(LS, MS, ET0, TOTXSEC, IERROR)
ELSE IF (INONVHS(LS,MS) == 0 .or. (INONVHS(LS,MS) == 1 .and. EVIB < 0.0D0)) THEN
  TOTXSEC = SPM(2,LS,MS)*((2.D00*BOLTZ*SPM(5,LS,MS)/(SPM(1,LS,MS)*VRR))**(SPM(3,LS,MS)-0.5D00))*SPM(6,LS,MS)*1.0d20
ELSE
  EVIBEV = EVIB/EVOLT  ! convert to eV
  IF (EVIBEV >= 6.9d0  ) EVIBEV = 6.9d0  ! extrapolate
  TOTXSEC = DLOG(ET0) - (CTOT(8) + EVIBEV*(CTOT(7) + EVIBEV*CTOT(6)))
  TOTXSEC = (CTOT(3)+CTOT(4)*EVIBEV*EVIBEV)*DTANH(CTOT(5)*TOTXSEC) + CTOT(2)*EVIBEV
  TOTXSEC = CTOT(1)*DEXP(TOTXSEC)
END IF
BMAX = DSQRT(TOTXSEC/PI)  !angstrom
CVR = TOTXSEC/1.0d20*VR

END SUBROUTINE CALC_TOTXSEC
