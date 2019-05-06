!
!*****************************************************************************
!
SUBROUTINE REMOVE_MOL(N)
!
!--remove molecule N and replaces it by NM
USE MOLECS
USE CALC
USE GEOM
USE GAS
!
IMPLICIT NONE
!
INTEGER :: N,NC,M,K

!--N the molecule number
!--M,K working integer
!
IF (N /= NM) THEN
  PX(N)=PX(NM)
  PV(1:3,N)=PV(1:3,NM)
  IF (MMRM > 0) PROT(N)=PROT(NM)
  IF (MMVM > 0) IPVIB(:,N)=IPVIB(:,NM)
  IPCELL(N)=IPCELL(NM)  !ABS(IPCELL(NM))  !--isebasti: removed ABS
  IPSP(N)=IPSP(NM)
  IPCP(N)=IPCP(NM)
  PTIM(N)=PTIM(NM)
END IF
NM=NM-1
!
RETURN
!
END SUBROUTINE REMOVE_MOL
