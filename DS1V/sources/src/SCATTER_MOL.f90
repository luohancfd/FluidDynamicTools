!
!*****************************************************************************
!
SUBROUTINE SCATTER_MOL(LS,MS,BMAX,BR,ET0,VRI,VR,VRC,VRCP,IDT)
!-- a wrapper of different scattering model
USE GAS, only: INONVHS
USE EXPCOL, only: EXPCOL_Scatter
IMPLICIT NONE
INTEGER,INTENT(IN) :: LS, MS
REAL(8),INTENT(IN) :: BMAX,BR,ET0,VRI, VR, VRC(3)
REAL(8) :: VRCP(3), B
INTEGER :: IDT
! VRI original relative speed
! ET0: initial collisional energy in eV
! VR, VRC new speed without scattering
! all the following functions should only change VRCP
IF (INONVHS(LS,MS) == 0 ) THEN
  CALL VHSS(LS, MS, BR, VR, VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 1) THEN
  B = DSQRT(BR) * BMAX
  CALL N2OScatter(LS, MS,B,VRI, VR,VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 2) THEN
  B = DSQRT(BR) * BMAX
  CALL EXPCOL_Scatter(LS,MS,B,ET0,VR,VRC,VRCP,IDT)
END IF

END SUBROUTINE SCATTER_MOL
