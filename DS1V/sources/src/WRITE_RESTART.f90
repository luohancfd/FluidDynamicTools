!
!*****************************************************************************
!
SUBROUTINE WRITE_RESTART
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE MFDSMC, only:IMF,IMFS,NMFpair,NMFANG, IMFpair, IMFdia,&
  MFRMASS, NMFER0, NMFET0,NMFEV0,NMFER,&
  NMFET,NMFEV,NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
!
IMPLICIT NONE
!
INTEGER :: ZCHECK
!
ZCHECK=1234567
!
101 CONTINUE
OPEN (7,FILE='PARAMETERS.DAT',FORM='UNFORMATTED',ERR=101) !--isebasti: replace binary by unformatted
WRITE(7) NCCELLS,NCELLS,MMRM,MMVM,MNM,MNRE,MNSR,MSP,MTBP,ILEVEL,MDIV,IRECOM,MMEX,MEX,ISF,NBINS,NSNAP,ITMAX !--isebasti: included ISF,NBINS,NSAP,ITMAX
WRITE(7) IMF, IMFS, NMFpair, IMFdia, nonVHS !--han: included these for MF model
CLOSE(7)
!
102 CONTINUE
OPEN (7,FILE='RESTART.DAT',FORM='UNFORMATTED',ERR=102) !--isebasti: replace binary by unformatted
WRITE (7)AC,AE,AVDTM,BC,BOLTZ,EVOLT,CCELL,CELL,CI,CLSEP,COLLS, &
         CPDTM,CR,CS,CSH,CSS,CSSS,CTM,CXSS,CST,BINS,BIN,PDFS,PDF,NSPDF,NDROT,NDVIB,DDIV,DPI,DTM,DTSAMP,DTOUT, &
         ENTMASS,ENTR,ER,ERROR,ERS,FDEN,FMA,FND,FNUM,FRACSAM,FSP,FP,FPR,FREM,FSPEC, &
         FTMP,FTIME,FRTMP,FVTMP,GASCODE,ICCELL,ICELL,ICREF,IFX,IMTS,IPCELL,IPCP, &
         IPSP,IPVIB,IRCD,IREA,IREM,ISECS,ISF,ISPEX,ISPR,ISPRC,ISPV,ISPVM,ISRCD,ITYPE,IVB,IWF,JCD, &
         JDIV,JREA,KP,LE,LIS,LP,LRS,ME,MOLSC,MP,MVER,NCCELLS,NCELLS, &
         NCIS,NDIV,NEX,NM,NMISAMP,NNC,NOUT,NREA,NRSP,NSAMP,NSPEX,NVER,OUTRAT,SAMPRAT,PI,PROT,PTIM,PV,PX, &
         REAC,RGFS,RMAS,REA,SP,SPEX,SPI,SPM,SPR,SPRC,SPV,SPVM,IDL,CPER,DPER,ENERS,ENERS0,IREV,SREAC,NPVIB,FPVIB,FPTEMP,IRM, &
         TCOL,CSCR,TDISS,TFOREX,TRECOMB,TREVEX,THBP,TISAMP,TPOUT,TREF,TLIM,TOTCOL,PCOLLS,TOTDUP,TOTMOV,     &
         TREACG,TREACL,TOUT,TPDTM,TREF,TSAMP,TSURF,VAR,VARS,VARSP,VELOB,VFX,VFY,UVFX,UFND,UFTMP,UVMP,VMP, &
         VMPM,VNMAX,VSURF,WCOLLS,WFM,XB,XREM,XVELS,YVELS,TNEX,NEVIB,FEVIB,QCTMODEL,IVMODEL,VIBEN,IGS,AMEG,EVREM,ZCHECK
!--isebasti: CST,BINS,BIN,PDFS,PDF,UVFX,IDL,CPER,DPER,ENERS,ENERS0,IREV,NPVIB,FPVIB,FPTEMP,UVFX,UFND,UFTMP,UVMP included

! Restart variable related to MFmodel
IF (MNRE >0 .AND. IMF .ne. 0) THEN
  WRITE(7) NMFANG, MFRMASS,IMFpair, ZCHECK
  ! the following dump is huge, thus we write it to a individual file
  IF (IMFS == 1 .AND. NMFpair > 0)THEN
103 CONTINUE
    OPEN (11,FILE='RESTART_MF_MODEL.DAT',FORM='UNFORMATTED',ERR=103)
    WRITE(11) NMFER0,NMFET0,NMFEV0,NMFER,NMFET,NMFEV,&
      NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR,ZCHECK
    CLOSE(8)
  END IF
END IF
WRITE(7) INONVHS, IQCTVT, ZCHECK
CLOSE(7)
!
WRITE (9,*) 'Restart files written'
!
RETURN
END SUBROUTINE WRITE_RESTART
