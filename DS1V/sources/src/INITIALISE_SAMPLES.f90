!
!*****************************************************************************
!
SUBROUTINE INITIALISE_SAMPLES
!
!--start a new sample
!
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE MOLECS
USE MFDSMC,only: IMF, IMFS, MF_INIT_SAMPLE
!
IMPLICIT NONE
!
!
NSAMP=0.
TISAMP=FTIME
NMISAMP=NM
ENERS=0.d0
COLLS=0.D00 ; WCOLLS=0.D00 ; CLSEP=0.D00
TCOL=0.D00 ; CSCR=0.D00
TDISS=0.D00   !--isebasti: uncommented
TRECOMB=0.D00 !--isebasti: uncommented
TFOREX=0.D00  !--isebasti: uncommented
TREVEX=0.D00  !--isebasti: uncommented
TREACG=0
TREACL=0
TNEX=0.D00    !--isebasti: uncommented
!
 CS=0. ; CSS=0. ; CSSS=0. ; CST=0.d0; BINS=0.d0; BIN=0.d0 !--isebasti: CST,BINS,BIN included
 CST(0,:)=1.d0; BINS(0,:,:)=1.d0; BIN(0,:)=1.d0           !--isebasti: to avoid dividing by zero
! Han: I don't think the following line makes sense
! CCELL(4,:)=SQRT(2.D00*BOLTZ*VAR(8,:)/SP(5,3))*SPM(2,3,3) !--isebasti: included



!
REAC=0.       !--isebasti: uncommented
SREAC=0.      !--isebasti: uncommented
!
NDISSOC=0   !used in qk
NRECOMB=0   !used in qk
NDISSL=0    !used in qk
NDROT=0     !bin counter
NDVIB=0     !bin
NDVIB(:,0,:,0)=1.d0 !to avoid dividing by zero
EVREM = 0.0d0
!
IF (MNRE > 0) NPVIB=0 !bin counter
IF (MNRE > 0) NEVIB=0 !bin counter
IF (MNRE > 0  .and. IMF .ne. 0 .and. IREAC == 2 .and. IMFS == 1) THEN
   CALL MF_INIT_SAMPLE()
 END IF
!
END SUBROUTINE INITIALISE_SAMPLES
