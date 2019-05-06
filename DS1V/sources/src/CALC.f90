!
!********************************************************************
!
MODULE CALC
!
!--declares the variables associated with the calculation
!
IMPLICIT NONE
!
INTEGER :: NVER,MVER,MOLSC,JCD,ISF,ISECS,IGS,IREM,IRECOM,NNC,IMTS,ERROR,ITCV,IEAA,IZV,&
           IFI,IRM,IADAPT,IJSR,IPDF,ITHP,IUBC,IREAC,IRELAX,IPRS,QCTMODEL !--isebasti: IFI,IRM,IADAPT,IJSR,IPDF,ITHP,IUBC,IREAC,IPRS included
REAL(KIND=8) :: FREM,XREM,FTIME,TLIM,PI,SPI,DPI,BOLTZ,FNUM,DTM,TREF,TSAMP,TOUT,AMEG,SAMPRAT,OUTRAT,TOTCOLI,TOTMOVI,&  !--isebasti: RANF removed
                DTSAMP,DTOUT,TPOUT,FRACSAM,TOTMOV,TOTCOL,PCOLLS,ENTMASS,CPDTM,TPDTM,TDISS,TRECOMB,TFOREX,TREVEX,TOTDUP,&
                AVOG,ELACK,EVOLT !--isebasti: ELACK included
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: VNMAX
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: TCOL,CSCR
INTEGER, ALLOCATABLE, DIMENSION(:) :: ISEED     !--isebasti: included
INTEGER :: nonVHS  !--han add
!
!--NVER.MVER the version number
!--AMEG the initial number of megabytes to be used by the program
!--MOLSC the target number of molecules per sampling cell
!--FREM fraction of molecule removal
!--XREM the coordinate at which removal commences
!--FTIME the flow time
!--TLIM the time at which the calculation stops
!--FNUM the number of real molecules represented by each simulated molecule
!--CPDTM the maximum number of collisions per time step (standard 0.2)
!--TPDTM the maximum number of sampling cell transit times of the flow per time step
!--TOTMOV total molecule moves
!--TOTCOL total collisions
!--PCOLLS total collisions in NSPDF cell
!--TOTDUP total duplicated collisions
!--TDISS dissociations since sample reset
!--TFOREX forward exchange reactions since sample reset
!--TREVEX backward exchange reactions since sample reset
!--TRECOMB recombinations since sample reset
!--ENTMASS the current entry mass of which a fraction FREM is to be removed
!--VNMAX(L) the maximum normal velocity component of species L
!--JCD  0 if chemical reactions are based on rate equations (GASCODE=6 only), 1 if they based on the Q-K model
!--TCOL species dependent collision counter
!--ISF 0,1 for steady, unsteady flow sampling
!--IRECOM 0 to neglect recombination reactions, 1 to include them
!--ISECS 0,1 for no secondary stream,a secondary stream that applies for positive values of x
!--IREM data item to set type of molecule removal
!--NNC 0 for normal collisions, 1 for nearest neighbor collisions
!--IMTS 0 for uniform move time steps, 1 for time steps that vary over the cells
!--IGS 0 for initial gas, 1 for stream(s) or reference gas, 2 for stream(s) or reference gas with TDF temperature distribution
!--ITCV 0 to base ternary collision volume on macroscopic temperature, 1 on translational collision temperature
!--IEAA 0 to base activation energy adjustments on macroscopic temperature, 1 on translational collision temperature
!--N.B. the coefficent a for tern. coll. vols. and act. en. adj. requre a transformation if ITCV or IEAA = 1
!--IZV  0 to base vibrational relaxation collision numbers on macroscopic temperature, 1 on quantized collision temperature
!--IFI 0 no forced ignition or >0 for forced ignition
!--IRM to select suboptions concerning GASCODE choice
!--IADAPT 0 for no cell adapting, 1 for adapting cells with VAR(3,N)/FND(1)>DCRIT,
!         2 for adapting cell with VAR(3,N)/MIN(VAR3,:)>DCRIT
!--IDT thread id; =0 in serial calculation but =0,1,2,...,mnth (maximum number of threads) in parallel cases (LOCAL VARIABLE)
!--ISEED(0:mnth) vector containing current number generator seed correspoding to each thread
!--IJSR initial seed that produces ISEED values
!--IUBC 0 for no BCs updating, 1 for updating BCs based on target velocities
!       2 for updating BCs based on target inlet pressures
!-- nonVHS 0 for VHS/VSS model
!          1 for VHS/VSS model+QCT N2O model
!          2 exponential model, fallback to 0 if parameters not found
!          3 special fix for MF-DSMC O2/O case
!
END MODULE CALC
