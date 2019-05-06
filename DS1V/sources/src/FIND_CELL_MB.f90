!
!*****************************************************************************
!
SUBROUTINE FIND_CELL_MB(X,NCC,NSC,TIM)
!
USE MOLECS
USE GEOM
USE CALC
!
IMPLICIT NONE
!
INTEGER :: NSC,NCC,ND
REAL(KIND=8) :: X,A,B,TIM
!
!--NCC collision cell number
!--NSC sampling cell number
!--X location
!--ND division number
!--DSC the ratio of the sub-division width to the division width
!--TIM the time
!
A=(XB(2)+VELOB*TIM-XB(1))/DFLOAT(NDIV)      !--new DDIV
ND=(X-XB(1))/A+0.99999999999999D00
B=XB(1)+DFLOAT(ND-1)*A
!
!the division is a level 0 sampling cell
NSC=-JDIV(0,ND)
NCC=NCIS*(X-B)/A+0.99999999999999D00
NCC=NCC+ICELL(NSC)
RETURN
WRITE (9,*) 'No cell for molecule at x=',X
STOP
END SUBROUTINE FIND_CELL_MB
