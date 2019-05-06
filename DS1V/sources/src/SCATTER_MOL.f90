!
!*****************************************************************************
!
SUBROUTINE SCATTER_MOL(LS,MS,BMAX,ET0,VRI,VR,VRC,VRCP,IDT)
!-- a wrapper of different scattering model
USE GAS, only: INONVHS
USE EXPCOL, only: EXPCOL_Scatter
IMPLICIT NONE
INTEGER,INTENT(IN) :: LS, MS
REAL(8),INTENT(IN) :: BMAX,ET0,VRI, VR, VRC(3)
REAL(8) :: VRCP(3)
INTEGER :: IDT
! VRI original relative speed
! VR, VRC new speed without scattering
! all the following functions should only change VRCP
IF (INONVHS(LS,MS) == 0 ) THEN
  CALL VHSS(LS, MS, VR, VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 1) THEN
  CALL N2OScatter(LS, MS, BMAX,VRI, VR,VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 2) THEN
  CALL EXPCOL_Scatter(LS,MS,BMAX,ET0,VR,VRC,VRCP,IDT)
END IF

END SUBROUTINE SCATTER_MOL
