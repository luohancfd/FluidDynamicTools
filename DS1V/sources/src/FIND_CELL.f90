!
!*****************************************************************************
!
SUBROUTINE FIND_CELL(X,NCC,NSC)
!
USE MOLECS
USE GEOM
USE CALC
!
IMPLICIT NONE
!
INTEGER :: N,L,M,NSC,NCC,ND
REAL(KIND=8) :: X,FRAC,DSC
!
!--NCC collision cell number
!--NSC sampling cell number
!--X location
!--ND division number
!--DSC the ratio of the sub-division width to the division width
!
ND=(X-XB(1))/DDIV+0.99999999999999D00
!
IF (JDIV(0,ND) < 0) THEN    !the division is a level 0 (no sub-division) sampling cell
  NSC=-JDIV(0,ND)
!  IF (IFX == 0)
  NCC=NCIS*(X-CELL(2,NSC))/(CELL(3,NSC)-CELL(2,NSC))+0.9999999999999999D00
  NCC=NCC+ICELL(NSC)
!  IF (NCC == 0) NCC=1
  RETURN
ELSE  !--the molecule is in a subdivided division
  FRAC=(X-XB(1))/DDIV-DFLOAT(ND-1)
  M=ND
  DO N=1,ILEVEL
    DSC=1.D00/DFLOAT(N+1)
    DO L=1,2  !over the two level 1 subdivisions
      IF (((L == 1).AND.(FRAC < DSC)).OR.((L == 2).AND.(FRAC >= DSC))) THEN
        M=JDIV(N-1,M)+L  !the address in JDIV
        IF (JDIV(N,M) < 0) THEN
          NSC=-JDIV(N,M)
          NCC=NCIS*(X-CELL(2,NSC))/(CELL(3,NSC)-CELL(2,NSC))+0.999999999999999D00
          IF (NCC == 0) NCC=1
          NCC=NCC+ICELL(NSC)
          RETURN
        END IF
      END IF
    END DO
    FRAC=FRAC-DSC
  END DO
END IF
WRITE (9,*) 'No cell for molecule at x=',X
STOP
END SUBROUTINE FIND_CELL
