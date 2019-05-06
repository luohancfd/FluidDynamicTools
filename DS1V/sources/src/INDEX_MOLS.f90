!
!*****************************************************************************
!
SUBROUTINE INDEX_MOLS
!
!--index the molecules to the collision cells
!
USE MOLECS
USE CALC
USE GEOM
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: N,M,K
!
!--N,M,K working integer
!
ICCELL(2,:)=0    !--isebasti: replace original loop by this
!
IF (NM /= 0) THEN
  DO N=1,NM      !--isebasti: no openmp; tests show best efficiency for 1 thread
    M=IPCELL(N)
    ICCELL(2,M)=ICCELL(2,M)+1
  END DO
!
  M=0
  DO N=1,NCCELLS
    ICCELL(1,N)=M
    M=M+ICCELL(2,N)
  END DO
!
  ICCELL(2,:)=0  !--isebasti: included
!
  DO N=1,NM
    M=IPCELL(N)
    ICCELL(2,M)=ICCELL(2,M)+1
    K=ICCELL(1,M)+ICCELL(2,M)
    ICREF(K)=N
  END DO
!
END IF
!
RETURN
!
END SUBROUTINE INDEX_MOLS
