!
!*****************************************************************************
!
SUBROUTINE INELASTIC_VT(LS, MS, EV, EC, DOF, P)
!
!--calculate P/P_max for VT process
! energy should be given in si
! DOF is the total dof of all mode in EC other than vibrational mode, only need for LB
! EXPCOL type collision model only handles VT
!
USE GAS,ONLY : INONVHS
USE EXPCOL,ONLY : EXPCOL_VT
IMPLICIT NONE
INTEGER, INTENT(IN) :: LS, MS
REAL(KIND=8), INTENT(IN) :: EV, EC, DOF
REAL(KIND=8) :: P

IF (INONVHS(LS,MS) == 2) THEN
  P = EXPCOL_VT(LS, MS, EV, EC)
ELSE
  ! use regular LB model
  P = (1.d0 - EV/EC)**(DOF*0.5d0 - 1.0d0)  ! bird eq 5.61
END IF
!
END SUBROUTINE INELASTIC_VT
