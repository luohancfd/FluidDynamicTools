!
!***********************************************************************
!*****************************MAIN PROGRAM******************************
!***********************************************************************
!
!--Version x.xx (must be set in NVER,MVER)
!
!--the version changes when there are major enhancements to the program
!--the first subscript must change if there is any alteration to the data
!--the second subscript is the release number
!
PROGRAM DS1
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE OMP_LIB
USE MFDSMC,only : IMF, IMFS, MF_CLEAN_AHO, IMFdia
!
IMPLICIT NONE
!
INTEGER :: IRUN,NSEED,N,I,J,IRETREM,NTHREADS, ITHREAD
REAL(KIND=8) :: WCLOCK(5)
INTEGER(KIND=8) :: COUNT0, COUNT1, COUNT_RATE
REAL(8) :: CALC_TIME
LOGICAL :: FILE_EXIST

NTHREADS = 1
!$OMP PARALLEL PRIVATE(ITHREAD)
!$   ITHREAD = omp_get_thread_num()
!$   IF (ITHREAD == 0) THEN
!$       NTHREADS = omp_get_num_threads()
!$   END IF
!$OMP END PARALLEL
!
NVER=1          !--for changes to basic architecture
MVER=22         !--significant changes, but must change whenever the data in DS1VD.DAT changes
!
!$ WCLOCK(:)=omp_get_wtime()
!
!--set constants
PI=3.1415926535897932D00
DPI=6.283185307179586D00
SPI=1.772453850905516D00
BOLTZ=1.38064852d-23
EVOLT=1.602176565d-19    !--isebasti: multiply to convert from electron to Joules
AVOG=6.022169D26
IJSR=123456789           !--isebasti: base seed for Ziggurate Random# Generator
!
!--set internal switches
IRUN=0       !initializing (do not change)
IADAPT=0     !initializing (do not change)
ELACK=0.d0   !initializing (do not change)
ITCV=1       !0 ternary collision volume depends on the macroscopic temperature; 1 collision temperature (used only for QK)
IEAA=1       !0 activation energy adjustments depend on the macroscopic temperature; 1 collision temperature (used only for QK)
IZV=0        !0 vibrational relaxation collision number depends on the Ttrans (PREFERABLE); 1 collision temperature
IPDF=1       !>0 sample PDFs for one particular cell (NSPDF)
NSPDF=1      !value can be modified in SAMPLE_PDF
NSCELLS=1    !number of sampling cells for vibrational pdfs
IENERS=1     !>0 sample total energy in domain (must be extended for GASCODE /= 8)
NBINS=10000 !60 !number of bins in velocity/speed PDF sampling
NSNAP=10000  !number of particles in snapshot files
ITHP=0       !>0 compute thermophoretic forces
IUBC=0       !>0 update boundary conditions
IPRS=2       !post-reaction vibrational levels are populated using 0 LB, 1 pre-reaction pdf sampling, 2 modified LB
ITMAX=4      !number of temperature intervlas for pre-reaction sampling (it might be redefined after reading restart files)
IRELAX=0     !=0 for isothermal relaxation tests: 0 moving/indexing are prohibited, 1 allowed
!=== Change of Zr and Zv
IRELAX_ROT=0 !=1 do nothing with Zr   =0 make Zr equal to 1E20
IRELAX_VIB=1 !=1 do nothing with Zv   =0 make Zv equal to 1E20
! This two variables are NOT saved
!==============
IREAC=0      !0 reaction on and rate sampling off (standard case)
!             1 reaction on and rate sampling on (samples only NSPDF cell, requires ISF>0)
!             2 reaction off and rate sampling on (called by run_ireac_ds1v.sh, samples only NSPDF cell, requires ISF=0 & NSEED>= 0)
QCTMODEL=1   !0 for TCE+LB+SHO vibrational levels
             !1 for TCE/MF+LB+AHO
             !2 for TCE/MF+MEQCT+AHO (MEQCT for O2+O and N2+O), IMF should also be turned on
             !3 for TCE+SSD+SSE+MEQCT+AHO (QCT SSD for O2+O, N2+O and SSE for N2+O)
IMFS=0     ! 0 for not sampling 1 for sampling
! WARNING: All above variables are hard coded, pay attenton before use it

! ----- WARNING: The following variables are overwriten by reading input ------------------
IMF=0      ! 0: do not use MF model
           ! 1: use MF+SHO
           ! 2: use MF+AHO(QCT ladder)
           ! 3: use MF+AHO(Morse Potential)
           ! 4: use MF+AHO(Morse Potential)+Effective translational energy
IMFdia=1   ! 0: for collision of same molecules, dissociate the one with higher vibrational energy
           ! 1: for collision of same molecules, dissociate the one gives lower threshold energy
nonVHS=2
!-- nonVHS 0 for VHS/VSS model
!          1 for VHS/VSS model+QCT N2O model
!          2 exponential model, fallback to 0 if parameters not found
!          3 special fix for MF-DSMC O2/O case
!-------------------------------------------------------------------------------------------
!
!--variables for vibratioal sampling
!
! TAG: SHOCK SAMPLE
! When you start sampling of vibrational levels, you should set appropriate value
! for NSCELLS, I, J, NSVEC
ALLOCATE (NSVEC(NSCELLS))
NSVEC=NSCELLS/2.+0.9999999d0
IF (NSCELLS > 1) THEN
  I=400              !user defined; values for shockwave sampling (nonreactive case)
  DO J=1,NSCELLS
    NSVEC(J)=I+J
  END DO
  NSVEC(1)=250       !user defined; override initial pdf sampling cell
  NSVEC(NSCELLS)=750 !user defined; override final pdf sampling cell
END IF
!
OPEN (9,FILE='DIAG.TXT',ACCESS='APPEND')
WRITE (9,*,IOSTAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'Stop the DS1.EXE that is already running and try again'
  STOP
ELSE
  WRITE (9,*) 'File DIAG.TXT has been opened'
END IF
!
DO WHILE ((IRUN < 1).OR.(IRUN > 3))
  WRITE (*,*) 'DS1 Version',NVER,'.',MVER
  WRITE (*,*) 'Enter 1 to continue the current sample or'
  WRITE (*,*) 'enter 2 to continue with a new sample, or'
  WRITE (*,*) 'enter 3 to start a new run :-'
  READ (*,*)  IRUN
END DO
!
IF (IRUN == 1) WRITE (9,*) 'Continuing the current run and sample'
IF (IRUN == 2) WRITE (9,*) 'Continuing the current run with a new sample'
IF (IRUN == 3) WRITE (9,*) 'Starting a new run'
!
IF (IRUN <= 2) THEN
  WRITE (*,*) 'Enter a seed (integer) for Random# Generator (a negative value calls run_multiple_cases.sh)'
  READ (*,*) NSEED
  IJSR=IJSR+100*NSEED
  CALL SET_ZIGGURAT
  WRITE (*,*) 'Setting Ziggurat Random# Generator: done'
  CALL READ_RESTART
  !ISF=1                    !to restart a case that is already ISF=0 as an unsteady-case
  IF (NSEED >= 9999) THEN
    ISF = 0
    IF (NSEED > 9999) THEN
      TPOUT = NSEED - 9999
      OUTRAT = TPOUT
      NSEED = 9999
    END IF
  END IF
  !IF (NSEED == 9999) ISF=0 !special code for shockwave sampling
  !CALL OXYGEN_NITROGEN     !special code to update reaction rates
  !CALL INIT_REAC           !special code to update reaction rates
  WRITE(9,*) ' ISF = ', ISF
  WRITE(9,*) ' NSEED = ',NSEED
  TSAMP=FTIME
  TOUT=FTIME
END IF
!
IF ((IRUN == 2).AND.(IVB == 0)) THEN
  WRITE (*,*) 'Enter 0 to continue with current cells;'
  WRITE (*,*) '      1 to adapt cells based on reference number density;'
  WRITE (*,*) '      2 to adapt cells based on mininum number density;'
  READ(*,*) IADAPT
!
  IF (IADAPT > 0) THEN
    WRITE (*,*) 'Adapting cells after beginning'
    CALL ADAPT_CELLS
    CALL INDEX_MOLS
    CALL WRITE_RESTART
  ELSE
    WRITE (*,*) 'Continuing with existing cells'
  END IF
!
  IF (IREM == 2) THEN
    WRITE (*,*) 'Enter 1 to reset the removal details, or 0 to continue with current details:-'
    READ(*,*) IRETREM
    IF (IRETREM == 1) THEN
      FREM=-1.D00
      DO WHILE ((FREM < -0.0001).OR.(FREM > 2.))
        WRITE (*,*) 'Enter the fraction of entering molecules that are removed:-'
        READ (*,*) FREM
        WRITE (*,999) FREM
      END DO
      WRITE (9,999) FREM
999   FORMAT (' The ratio of removed to entering molecules is ',G15.5)
      IF (FREM > 1.D-10) THEN
      XREM=XB(1)-1.
      DO WHILE ((XREM < XB(1)-0.0001).OR.(XREM > XB(2)+0.0001))
        WRITE (*,*) 'Enter x coordinate of the upstream removal limit:-'
        READ (*,*) XREM
        WRITE (*,998) XREM,XB(2)
      END DO
      WRITE (9,998) XREM,XB(2)
998   FORMAT (' The molecules are removed from ',G15.5,' to',G15.5)
      END IF
    END IF
  END IF
!
  CALL INITIALISE_SAMPLES
!
END IF
!
IF (IRUN == 3) THEN
  WRITE (*,*) 'Enter a seed (integer) for Random# Generator (a negative value calls run_multiple_cases.sh)'
  READ (*,*) NSEED
  IJSR=IJSR+100*NSEED
  CALL SET_ZIGGURAT
  WRITE (*,*) 'Setting Ziggurat Random# Generator: done'
  WRITE (*,*) 'Enter 0 to continue with current cells;'
  WRITE (*,*) '      1 to adapt cells based on reference number density;'
  WRITE (*,*) '      2 to adapt cells based on mininum number density;'
  READ (*,*) IADAPT
  CALL READ_DATA
  IF (IREAC == 0 .and. IRELAX == 0) THEN !extracting MeQct Zv(T) values
    WRITE (*,*) 'Enter temperature value for reaction rate sampling'
    READ (*,*) FTMP(1)
    READ (*,*) FVTMP(1)
    FTMP=FTMP(1)
    FRTMP=FTMP(1)
    FVTMP = FVTMP(1)
    !FVTMP=500.d0 !1000.d0/500.d0 !extracting MeQct Zv(T)
    TSURF=FTMP
    !FND=101325.d0/(BOLTZ*FTMP(1))
    !IF (AE(1) <= 0.) FND=2.5d25*(FTMP/300.d0)  !101325.d0/(BOLTZ*FTMP(1))
  END IF
  IF (IREAC == 2 .and. IRELAX == 0) THEN
    WRITE (*,*) 'Enter temperature value for reaction rate sampling'
    READ (*,*) FTMP(1)
    READ (*,*) FVTMP(1)
!    FTMP(1)=5.d3 !testing MeQct Nonequilibrium Rates
    FTMP=FTMP(1)
    FRTMP=FTMP(1)
    FVTMP=FVTMP(1)
    TSURF = FTMP
    FTMP0 = FTMP(1)
    FVTMP0 = FVTMP(1)
    ! FND=101325.d0/(BOLTZ*FTMP(1))
    !  IF (AE(1) <= 0.) FND=2.5d25*(FTMP/300.d0)  !101325.d0/(BOLTZ*FTMP(1))
  END IF
  CALL SET_INITIAL_STATE
END IF
!
!$ write(*,*) 'Pre-processing wallclock time (s):  ', omp_get_wtime()-WCLOCK(1)
!
TOTCOLI=TOTCOL
TOTMOVI=TOTMOV
!
IF (NOUT == -1) THEN
  DO I=1,1
    FTIME=FTIME+DTM
    CALL MOLECULES_MOVE
    IF (ITYPE(1)==0.OR.ITYPE(2)==0) CALL MOLECULES_ENTER
    CALL INDEX_MOLS
    IF (IPDF > 0) CALL SAMPLE_PDF
    CALL COLLISIONS
    CALL SAMPLE_FLOW
  END DO
  CALL OUTPUT_RESULTS
!
  IF ((IRUN == 3) .AND. (IADAPT > 0)) THEN
    WRITE(*,*) 'Adapting cells at beginning'
    CALL ADAPT_CELLS
    CALL INDEX_MOLS
    CALL INITIALISE_SAMPLES
    DO I=1,10
      FTIME=FTIME+DTM
      CALL MOLECULES_MOVE
      IF (ITYPE(1)==0.OR.ITYPE(2)==0) CALL MOLECULES_ENTER
      CALL INDEX_MOLS
      IF (IPDF > 0) CALL SAMPLE_PDF
      CALL COLLISIONS
      CALL SAMPLE_FLOW
    END DO
    NOUT=NOUT-1
    CALL OUTPUT_RESULTS
    WRITE(*,*) 'Adapting cells at beginning: done'
  END IF
!
  CALL WRITE_RESTART
  CALL INITIALISE_SAMPLES
  WRITE(*,*) 'Writing outputs for initial condition: done'
END IF
!
IF (NSEED < 0) THEN
  IF (IREAC == 2) THEN
    WRITE (*,*) 'Set up does not support multiple copies'
    WRITE (9,*) 'Set up does not support multiple copies'
    STOP
  END IF
  WRITE (*,*) 'Creating multiple cases (different seeds) from RESTART file'
  WRITE (9,*) 'Creating multiple cases (different seeds) from RESTART file'
  CALL SYSTEM ("./run_multiple_cases.sh")
  STOP
END IF
CLOSE(9)
!flush diag.txt
!
!------------------------------------------------------------------------------
!
INQUIRE(FILE="RunningTime.DAT", EXIST=FILE_EXIST)
IF (.NOT. FILE_EXIST) THEN
  OPEN(119, FILE="RunningTime.DAT", STATUS="REPLACE")
  WRITE(119,"(10(A,1X,I2,2X))") "# IRELAX:",IRELAX,"IREAC:",IREAC,"nonVHS:",nonVHS,"GASCODE:",GASCODE
  WRITE(119,"(10(A,1X,I2,2X))") "# QCTMODEL:",QCTMODEL,"IMF:",IMF,"IMFdia:",IMFdia,"IMFS:",IMFS
  WRITE(119,"(A,G14.6)") '# SAMPRAT = ', SAMPRAT
  WRITE(119,"(A,G14.6)") '# OUTRAT = ', OUTRAT
  WRITE(119,"(A,I5)") '# Nthreads = ', NTHREADS
  WRITE(119,"(10(A,1X))") 'VARIABLES = ','"NOUT",', '"CALC_TIME (s)",', '"FTIME (s)",', &
                    &   '"DTM (s)",', '"NMOL,"', '"NSAMP",', '"TISAMP",','"TPOUT"'
ELSE
  OPEN(119, FILE="RunningTime.DAT", POSITION="APPEND")
  WRITE(119, "(A)") "# Running is restarted: Sample Run"
  WRITE(119,"(A,G14.6)") '# SAMPRAT = ', SAMPRAT
  WRITE(119,"(A,G14.6)") '# OUTRAT = ', OUTRAT
  WRITE(119,"(A,G14.6)") '# TPOUT = ', TPOUT
  WRITE(119,"(A,G14.6)") '# ISF,IPDF = ', ISF,IPDF
END IF
CLOSE(119)
CALL SYSTEM_CLOCK(COUNT=COUNT0, COUNT_RATE=COUNT_RATE)


DO WHILE (FTIME < TLIM)
  OPEN (9,FILE='DIAG.TXT',ACCESS='APPEND')
  IF (TPOUT >= 1000) THEN
    OPEN(91, FILE="TPOUT_LOG.TXT",ACCESS="APPEND")
  END IF
!
!$ WCLOCK(2)=omp_get_wtime()
!
  DO N=1,INT(TPOUT)
    !
    IF (MOD(N,50) == 0 .and. TPOUT >= 1000) THEN
      WRITE(91,*) " TPOUT = ", N, "/", TPOUT
      CALL FLUSH(91)
    END IF
    DO I=1,INT(SAMPRAT)
      FTIME=FTIME+DTM
!      !$ WCLOCK(1)=omp_get_wtime()
      IF(IREAC <= 1 .AND. IRELAX > 0) CALL MOLECULES_MOVE
!      !$ WRITE(*,*) 'MOVE wallclock time (s):    ', omp_get_wtime()-WCLOCK(1)
!
!      !$ WCLOCK(1)=omp_get_wtime()
      IF ((ITYPE(1) == 0).OR.(ITYPE(2) == 0)) CALL MOLECULES_ENTER
!      !$ WRITE(*,*) 'ENTER wallclock time (s):   ', omp_get_wtime()-WCLOCK(1)
!
!      !$ WCLOCK(1)=omp_get_wtime()
      IF(IREAC <= 1 .AND. IRELAX > 0) CALL INDEX_MOLS
!      !$ WRITE(*,*) 'INDEX wallclock time (s):   ', omp_get_wtime()-WCLOCK(1)
!
!      !$ WCLOCK(1)=omp_get_wtime()
      CALL COLLISIONS
!      !$ WRITE(*,*) 'COLL wallclock time (s):    ', omp_get_wtime()-WCLOCK(1)
    END DO
!
!    !$ WCLOCK(1)=omp_get_wtime()
    CALL SAMPLE_FLOW
!    !$ WRITE(*,*) 'SAMPLE wallclock time (s):   ', omp_get_wtime()-WCLOCK(1)
!
!    !$ WCLOCK(1)=omp_get_wtime()
    IF (IPDF > 0) THEN
      IF (ISF == 0) THEN
        CALL SAMPLE_PDF
      ELSE IF ((ISF == 1).AND.(DBLE(N/TPOUT) >= (1.d0-FRACSAM))) THEN
        CALL SAMPLE_PDF
      ELSE IF (IPDF == 2 .AND. TPOUT < 10) THEN
        CALL SAMPLE_PDF
      ELSE IF ((ISF == 2).AND.(N >= INT(TPOUT-10))) THEN
        CALL SAMPLE_PDF
      END IF
    END IF
!    !$ WRITE(*,*) 'PDF wallclock time (s):      ', omp_get_wtime()-WCLOCK(1)
!
!    !$ WCLOCK(1)=omp_get_wtime()
    IF (ITHP > 0) CALL THERMOPHORETIC
!    !$ WRITE(*,*) 'THERMO wallclock time (s):   ', omp_get_wtime()-WCLOCK(1)
!    I need to create a variable with updated thermoph forces so that there is no problem on reseting CST samples
!
!    !$ WCLOCK(1)=omp_get_wtime()
    IF (N == INT(TPOUT*FRACSAM)*(N/INT(TPOUT*FRACSAM))) THEN
!      write(*,*) 'TPOUT,N,N/TPOUT',TPOUT,N,N/TPOUT
      !IF (ISF  > 0) CALL UPDATE_MP !need to be updated
      IF (IUBC > 0) CALL UPDATE_BC
      IF ((ISF > 0).AND.(N/TPOUT <= (1.d0-FRACSAM))) CALL INITIALISE_SAMPLES
    END IF
!    !$ WRITE(*,*) 'UPDATE_MP wallclock time (s):', omp_get_wtime()-WCLOCK(1)
!
    IF ((ISF == 2).AND.(N == INT(TPOUT-10))) CALL INITIALISE_SAMPLES
    CALL FLUSH(6)
!
  END DO
!
  !$ WCLOCK(1)=omp_get_wtime()
  CALL OUTPUT_RESULTS
  CALL SYSTEM_CLOCK(COUNT=COUNT1)
  CALC_TIME = dble(COUNT1-COUNT0) / dble(COUNT_RATE)

  ! io
  CALL FLUSH(9)
  OPEN(119, FILE="RunningTime.DAT", POSITION="APPEND")
  WRITE(119,118) NOUT, CALC_TIME, FTIME, DTM, NM, NSAMP, TISAMP, TPOUT
  CLOSE(119)
118  FORMAT(I10,2X, 10(G14.6,2X))

  IF (IREAC .ne. 2 .AND. IMFS .ne. 1) THEN
    ! WARNING: do not write restart for ireac = 2, file is too huge
    CALL WRITE_RESTART
  END IF
  IF (ISF > 0) CALL INITIALISE_SAMPLES
  IF (IREAC == 2) THEN
    IF ((NOUT < 2).OR.(REAC(1) == 0)) CALL INITIALISE_SAMPLES
   !IF ((NOUT > 100).OR.(REAC(1) > 1.d9)) STOP
    IF (REAC(1) > 1.d9) STOP
  END IF
  IF ((IREM == 2).AND.(XREM < XB(1)).AND.(VAR(3,INT(NCELLS*.666))/FND(1) > 2.d0)) THEN
    WRITE (9,*) 'Shockwave reached ~center of domain'
    WRITE (*,*) 'Shockwave reached ~center of domain'
    STOP
  END IF
  !$ WRITE(*,*) 'OUTPUT wallclock time (s):  ', omp_get_wtime()-WCLOCK(1)
  !$ WRITE(*,*) 'LOOP wallclock time (s):    ', omp_get_wtime()-WCLOCK(2)
  !$ WRITE(9,*) 'LOOP wallclock time (s):    ', omp_get_wtime()-WCLOCK(2)
  WRITE(*,*)
  WRITE(9,*)
  CLOSE(9)
  IF (TPOUT >= 1000) THEN
    CLOSE(91)
  END IF
!
END DO
CALL MF_CLEAN_AHO()
!
STOP
END PROGRAM DS1
