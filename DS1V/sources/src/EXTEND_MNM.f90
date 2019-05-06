!
!*****************************************************************************
!
SUBROUTINE EXTEND_MNM(FAC)
!
!--the maximum number of molecules is increased by a specified factor
!--the existing molecules are copied TO disk storage
!
USE MOLECS
USE CALC
USE GAS
!
IMPLICIT NONE
!
INTEGER :: M,N,MNMN
REAL(8) :: FAC
!
!--M,N working integers
!--MNMN extended value of MNM
!--FAC the factor for the extension
MNMN=FAC*MNM
WRITE (*,*) 'Maximum number of molecules is to be extended from',MNM,' to',MNMN
WRITE (*,*) '( if the additional memory is available !! )'
OPEN (7,FILE='EXTMOLS.SCR',FORM='UNFORMATTED') !--isebasti: replace binary by unformatted
WRITE (*,*) 'Start write to disk storage'
DO N=1,MNM
  IF (MMVM > 0) THEN
    WRITE (7) PX(N),PTIM(N),PROT(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N),(IPVIB(M,N),M=0,MMVM)
  ELSE
    IF (MMRM > 0) THEN
      WRITE (7) PX(N),PTIM(N),PROT(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N)
    ELSE
      WRITE (7) PX(N),PTIM(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N)
    END IF
  END IF
END  DO
WRITE (*,*) 'Disk write completed'
CLOSE (7)
IF (MMVM > 0) THEN
  DEALLOCATE (PX,PTIM,PROT,PV,IPSP,IPCELL,ICREF,IPCP,IPVIB,STAT=ERROR)
ELSE
  IF (MMRM > 0) THEN
    DEALLOCATE (PX,PTIM,PROT,PV,IPSP,IPCELL,ICREF,IPCP,STAT=ERROR)
  ELSE
    DEALLOCATE (PX,PTIM,PV,IPSP,IPCELL,ICREF,IPCP,STAT=ERROR)
  END IF
END IF
IF (ERROR /= 0) THEN
  WRITE (*,*)'PROGRAM COULD NOT DEALLOCATE MOLECULES',ERROR
!  STOP
END IF
!
WRITE (*,*) 'Molecule arrays have been deallocated'
!
IF (MMVM > 0) THEN
  ALLOCATE (PX(MNMN),PTIM(MNMN),PROT(MNMN),PV(3,MNMN),IPSP(MNMN),IPCELL(MNMN),&
             ICREF(MNMN),IPCP(MNMN),IPVIB(0:MMVM,MNMN),STAT=ERROR)
ELSE
  IF (MMRM > 0) THEN
    ALLOCATE (PX(MNMN),PTIM(MNMN),PROT(MNMN),PV(3,MNMN),IPSP(MNMN),IPCELL(MNMN),ICREF(MNMN),IPCP(MNMN),STAT=ERROR)
  ELSE
    ALLOCATE (PX(MNMN),PTIM(MNMN),PV(3,MNMN),IPSP(MNMN),IPCELL(MNMN),ICREF(MNMN),IPCP(MNMN),STAT=ERROR)
  END IF
END IF
IF (ERROR /= 0) THEN
  WRITE (*,*)'PROGRAM COULD NOT ALLOCATE SPACE FOR EXTEND_MNM',ERROR
!  STOP
END IF
!
PX=0.; PTIM=0.; PV=0.; IPSP=0; IPCELL=0; ICREF=0; IPCP=0
IF (MMRM > 0) PROT=0.
IF (MMVM > 0) IPVIB=0
!--restore the original molecules
OPEN (7,FILE='EXTMOLS.SCR',FORM='UNFORMATTED') !--isebasti: replace binary by unformatted
WRITE (*,*) 'Start read back from disk storage'
DO N=1,MNM
  IF (MMVM > 0) THEN
    READ (7) PX(N),PTIM(N),PROT(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N),(IPVIB(M,N),M=0,MMVM)
  ELSE
    IF (MMRM > 0) THEN
      READ (7) PX(N),PTIM(N),PROT(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N)
    ELSE
      READ (7) PX(N),PTIM(N),(PV(M,N),M=1,3),IPSP(N),IPCELL(N),ICREF(N),IPCP(N)
    END IF
  END IF
END DO
WRITE (*,*) 'Disk read completed'
CLOSE (7,STATUS='DELETE')
!
MNM=MNMN
!
RETURN
END SUBROUTINE EXTEND_MNM
