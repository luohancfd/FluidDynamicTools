!*************************************************************************************  
! 
MODULE POST  
!  
!--declares variables associated with post processing
!  
IMPLICIT NONE  
!
INTEGER :: IRUN,FRUN,NRUNS,CTIME,NSS  
INTEGER, ALLOCATABLE :: NMS(:),ISNAP(:,:)
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: PSNAP(:,:),VSNAP(:,:)
CHARACTER (LEN=4) :: E,NRUN
CHARACTER (LEN=256) :: DIR,DIRCASE,DIROUT
!
!--IRUN initial run to be averaged  
!--FRUN final run to be averaged
!--NRUNS total number of runs to be averaged  
!
!--auxiliar variables
!  
INTEGER :: IVNSAMP,IVNM,IVCTIME,IVNDISSOC,IVNRECOMB,IVNTSAMP
INTEGER, ALLOCATABLE :: IVNMS(:),IVIPSP(:,:)
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IVTREACG,IVTREACL
REAL(KIND=8):: DVTISAMP,DVFTIME,DVTOTMOV,DVTOTCOL,DVTOTDUP,DVTOTCOLI,DVTOTMOVI,DVSTARTIME,&
               DVOVTEMP,DVOVTTEMP,DVOVRTEMP,DVOVVTEMP,DVDTM,DVUVFX(2)
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: DVSTEMP,DVTRANSTEMP,DVROTTEMP,DVVIBTEMP,DVREAC 
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: DVTCOL,DVVARS,DVVAR,DVCST,DVPDF,DVBIN,DVPV,DVPX 
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: DVVARSP,DVPDFS,DVBINS
!
END MODULE POST  
!
!          BELOW THIS LINE, JUST COPY MODULES FROM DS1_openmp.f90!  
!*************************************************************************************  
! 
MODULE MOLECS  
!  
!--declares the variables associated with the molecules  
!  
IMPLICIT NONE  
!  
INTEGER, ALLOCATABLE, DIMENSION(:) :: IPCELL,IPSP,ICREF,IPCP  
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPVIB  
INTEGER :: NM,MNM  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: PV  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: PX,PTIM,PROT  
!  
!--PX(N) position (x coordinate) of molecule N  
!--PTIM(N) molecule time  
!--PROT(N) rotational energy  
!--PV(1-3,N) u,v,w velocity components  
!--IPSP(N) the molecular species  
!--IPCELL(N) the collision cell number  
!--ICREF the cross-reference array (molecule numbers in order of collision cells)  
!--IPCP(N) the code number of the last collision partner of molecule   
!--IPVIB(K,N) level of vibrational mode K of molecule N  
!--NM number of molecules  
!--MNM the maximum number of molecules  
!  
END MODULE MOLECS  
!  
!*************************************************************************************  
!  
MODULE GEOM  
!  
!--declares the variables associated with the flowfield geometry and cell structure  
!  
IMPLICIT NONE  
!  
INTEGER :: NCELLS,NCCELLS,NCIS,NDIV,MDIV,ILEVEL,IFX,JFX,IVB,IWF      
INTEGER, DIMENSION(2) :: ITYPE              
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICELL  
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICCELL,JDIV  
REAL(KIND=8) :: DDIV,XS,VELOB,WFM,AWF  
REAL(KIND=8), DIMENSION(2) :: XB  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: CELL,CCELL  
!  
!--XB(1), XB(2) the minimum, maximum x coordinate  
!--DDIV the width of a division  
!--ITYPE(K) the tpe of boundary at the minimum x (K=1) and maximum x (K=2) boundaries  
!          0 for a stream boundary  
!          1 for a plane of symmetry  
!          2 for a solid surface  
!          3 for a vacuum  
!--NCELLS the number of sampling cells  
!--NCCELLS the number of collision cells  
!--NCIS the number of collision cells in a sampling cell  
!  MDIV the maximum number of sampling cell divisions at any level of subdivision  
!--IVB 0,1 for stationary, moving outer boundary  
!--IWF 0 for no radial weighting factors, 1 for radial weighting factors  
!--WFM, set in data as the maximum weighting factor -1, then divided by the maximum radius  
!--AWF overall ratio of real to weighted molecules   
!--VELOB the speed of the outer boundary  
!--JDIV(N,M) (-cell number) or (start address -1 in JDIV(N+1,M), where M is MDIV  
!--JFX  IFX+1  
!--CELL(M,N) information on sampling cell N  
!    M=1 x coordinate  
!    M=2 minimum x coordinate  
!    M=3 maximum x cooedinate   
!    M=4 volume  
!--ICELL(N) number of collision cells preceding those in sampling cell N  
!--CCELL(M,N) information on collision cell N  
!    M=1 volume  
!    M=2 remainder in collision counting  
!    M=3 half the allowed time step  
!    M=4 maximum value of product of cross-section and relative velocity   
!    M=5 collision cell time  
!--ICCELL(M,N) integer information on collision cell N  
!    M=1 the (start address -1) in ICREF of molecules in collision cell N  
!    M=2 the number of molecules in collision cell N  
!    M=3 the sampling cell in which the collision cell lies  
!  
END MODULE GEOM  
!  
!********************************************************************  
!  
MODULE GAS  
!  
!--declares the variables associated with the molecular species and the stream definition  
!  
IMPLICIT NONE  
!  
REAL(KIND=8) :: RMAS,CXSS,RGFS,VMPM,FDEN,FPR,FMA,FP,CTM  
REAL(KIND=8), DIMENSION(2) :: FND,FTMP,FRTMP,FVTMP,VFX,VFY,TSURF,FSPEC,VSURF,UVFX,UFND,UFTMP !--isebasti: UVFX,UFND,UFTMP,UVMP,FRTMP included  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CI,AE,AC,BC,ER,ERS,CR,TNEX,PSF  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: FSP,SP,SPR,SPV,REA,THBP,VMP,UVMP  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: SPM,SPVM,ENTR  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: SPEX,SPRC  
INTEGER :: MSP,MMVM,MMRM,MNRE,MNSR,MTBP,GASCODE,MMEX,MEX  
INTEGER, ALLOCATABLE, DIMENSION(:) :: ISP,ISPV,LE,ME,KP,LP,MP  
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISPR,IREA,NREA,NRSP,LIS,LRS,ISRCD,ISPRC,TREACG,TREACL,NSPEX  
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ISPVM,NEX,IRCD,JREA  
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ISPEX  
!  
!--MSP the number of molecular species  
!--MMVM the maximum number of vibrational modes of any species  
!--MEX number of exchange or chain reactions  
!--MMEX the maximum number of exchange reactions involving the same precollision pair of molecules  
!--MMRM 0 if gass is completely monatomic, 1 if some species have rotation  
!--MNRE the number of gas phase chemical reactions (old scheme)  
!--MNSR the number oF surface reactions  
!--MTBP the number of third body tables  
!--SP(1,L) the reference diameter of species L         
!--SP(2,L) the reference temperature of species L  
!--SP(3,L) the viscosity-temperature power law of species L  
!--SP(4,L) the reciprocal of the VSS scattering parameter  
!--SP(5,L) molecular mass of species L  
!--ISPR(1,L) number of rotational degrees of freedom of species L  
!--ISPR(2,L) 0,1 for constant, polynomial rotational relaxation collision number  
!--SPR(1,L) constant rotational relaxation collision number of species L  
!           or the constant in a second order polynomial in temperature  
!--SPR(2,L) the coefficient of temperature in the polynomial  
!--SPR(3,L) the coefficient of temperature squared in the polynomial  
!--SPM(1,L,M) the reduced mass for species L,M  
!--SPM(2,L,M) the reference collision cross-section for species L,M  
!--SPM(3,L,M) the mean value of the viscosity-temperature power law  
!--SPM(4,L,M) the reference diameter for L,M collisions  
!--SPM(5,L,M) the reference temperature for species L,M  
!--SPM(6,L,M) reciprocal of the gamma function of (5/2-w) for species L,M  
!--SPM(7,L,M) rotational relaxation collision number for species L,M, or const in polynomial  
!--SPM(8,L,M) reciprocal of VSS scattering parameter  
!--ISPV(L) the number of vibrational modes  
!--SPVM(1,K,L) the characteristic vibrational temperature  
!--SPVM(2,K,L) constant Zv, or reference Zv for mode K  
!--SPVM(3,K,L) -1. for constant Zv, or reference temperature  
!--SPVM(4,K,L) the characteristic dissociation temperature  
!--ISPVM(1,K,L) the species code of the first dissociation product  
!--ISPVM(2,K,L) the species code of the second dissociation product  
!--ISPRC(L,M) the species of the recombined molecule from species L and M  
!--SPRC(1,L,M) the constant a in the ternary collision volume  
!--SPRC(2,L,M) the temperature exponent b in the ternary collision volume  
!--NSPEX(L,M) the number of exchange reactios with L,M as the pre-collision species  
!--in the following variables, J is the reaction number (1 to NSPEX(L,M))  
!--ISPEX(J,0,L,M) the species that splits in an exchange reaction  
!--ISPEX(J,1,L,M) the post collision species of the molecule that splits in an L,M exchange reaction (all ISPEX are set to 0 if no exchange reaction)  
!--ISPEX(J,2,L,M) the post collision species of either the atom or the molecule that does not split in an L,M exchange reaction  
!--ISPEX(J,3,L,M) vibrational mode that is associated with this reaction  
!--ISPEX(J,4,L,M) degeneracy of this reaction  
!--SPEX(1,J,L,M) the constant a in the activation energy  
!--SPEX(2,J,L,M) the temperature exponent b in the activation energy  
!--SPEX(3,J,L,M)  for the heat of reaction  
!--TNEX(N) total number of exchange reaction N   
!--NEX(N,L,M) the number of exchange or chain reactions in L,M collisions  
!  
!--The reaction information for the Q-K reaction model (JCD=1) is in the above  
!--The reaction information for the traditional reaction model (JCD =0) follows   
!  
!--IRCD(N,M,L) the code number of the Nth reaction between species M,L   
!--IREA(1,N),IREA(2,N) the codes of the pre-reaction species  
!--JREA(1,N,K) the code of the Kth post reaction species from first mol.  
!--JREA(2,N,K) similarly for second, but -ve indicates 3rd body table  
!--NREA(1,N),NREA(2,N) the number of post-reaction molecules from each  
!--NRSP(L,M) the number of chemical reactions for L-M precollision  
!--REA(1,N) number of contributing internal degrees of freedom  
!--REA(2,N) the activation energy for reaction N  
!--REA(3,N) the constant in the reaction rate  
!--REA(4,N) temperature exponent in the rate equation  
!--REA(5,N) heat of reaction (+,- for exothermic,endothermic)  
!--REA(6,N) non-energy terms in the steric factor expressions (6.10,14)  
!--LE(N) the species code of the first molecule  
!--ME(N) the species code of the second molecule  
!--for three post-collision particles  
!----KP(N) the first post-reaction species code  
!----LP(N) the second post-reaction species code  
!----MP(N) the third post-reaction species code  
!--for one post-collision particle  
!----KP(N)=-1  
!----LP(N) the species code of the molecule  
!----MP(N) if positive the code of the third body in a recombination  
!          or, if negative, the code of the third-body table  
!--for two post-collision particles  
!----KP(N)=0  
!----LP(N) the first post-reaction species code  
!----MP(N) the second post-reaction species code  
!----CI(N) the number of internal degrees of freedom that contribute  
!----AE(N) the activation energy for the reaction  
!----AC(N) the pre-exponential parameter  
!----BC(N) the temperature exponent in the rate equation  
!----ER(N) the energy of the reaction  
!--THBP(N,L) third body efficiency in table N of species L   
!--RMAS reduced mass for single species case  
!--CXSS reference cross-section for single species case  
!--RGFS reciprocal of gamma function for single species case  
!--for the following, J=1 for the reference gas and/or the minimum x boundary, J=2 for the secondary sream at maximum x boundary  
!--FND(J) stream or reference gas number density  
!--FTMP(J) stream temperature  
!--FRTMP(J) the rotational temperature in the freestream  
!--FVTMP(J) the vibrational temperature in the freestream  
!--VFX(J)  the x velocity components of the stream  
!--VFY(J) the y velocity component in the stream  
!--FSP(N,J)) fraction of species N in the stream  
!--FMA stream Mach number  
!--VMP(N,J) most probable molecular velocity of species N at FTMP(J)  
!--VMPM the maximum value of VMP in stream 1  
!--ENTR(M,L,K) entry/removal information for species L at K=1 for 1, K=2 for XB(2)  
!    M=1 number per unit time  
!    M=2 remainder  
!    M=3 speed ratio  
!    M=4 first constant  
!    M=5 second constant  
!    M=6 the maxinum normal velocity component in the removal zone (> XREM)  
!--LIS(1,N) the species code of the first incident molecule  
!--LIS(2,N) the species code of the second incident molecule (0 if none)  
!--LRS(1,N) the species code of the first reflected molecule  
!--LRS(2,N) the species code of the second reflected molecule (0 if none)  
!--LRS(3,N) the species code of the third reflected molecule (0 if none)  
!--LRS(4,N) the species code of the fourth reflected molecule (0 if none)  
!--LRS(5,N) the species code of the fifth reflected molecule (0 if none)  
!--LRS(6,N) the species code of the sixth reflected molecule (0 if none)  
!--ERS(N) the energy of the reaction (+ve for recombination, -ve for dissociation)   
!--NSRSP(L) number of surface reactions that involve species L as incident molecule  
!--ISRCD(N,L) code number of Nth surface reaction with species L as incident molecule    
!--CTM approximate mean collision time in stream (exact for simple gas)  
!--FP approximate mean free path  
!--FDEN stream 1 density    
!--FPR stream 1 pressure  
!--FMA stream 1 Mach number  
!--RMAS reduced mass for single species case  
!--CXSS reference cross-section for single species case  
!--RGFS reciprocal of gamma function for single species case  
!--CR(L) collision rate of species L  
!--TREACG(N,L) the total number of species L gained from reaction type N=1 for dissociation, 2 for recombination, 3 for forward exchange, 4 for reverse exchange  
!--TREACL(N,L) the total number of species L lost from reaction type N=1 for dissociation, 2 for recombination, 3 for forward exchange, 4 for reverse exchange  
!         
END MODULE GAS  
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
           IFI,IRM,IADAPT,IJSR,IPDF,ITHP,IUBC !--isebasti: IFI,IRM,IADAPT,IJSR,IPDF,ITHP,IUBC included  
REAL (KIND=8) :: FREM,XREM,FTIME,TLIM,PI,SPI,DPI,BOLTZ,FNUM,DTM,TREF,TSAMP,TOUT,AMEG,SAMPRAT,OUTRAT,TOTCOLI,TOTMOVI,&  !--isebasti: RANF removed
                 DTSAMP,DTOUT,TPOUT,FRACSAM,TOTMOV,TOTCOL,ENTMASS,CPDTM,TPDTM,TDISS,TRECOMB,TFOREX,TREVEX,TOTDUP,AVOG  
REAL (KIND=8), ALLOCATABLE, DIMENSION(:) :: VNMAX  
REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:) :: TCOL
INTEGER, ALLOCATABLE, DIMENSION(:) :: ISEED !--isebasti: included  
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
!
!  
END MODULE CALC  
!  
!********************************************************************  
!  
MODULE OUTPUT  
!  
!--declares the variables associated with the sampling and output  
!  
IMPLICIT NONE  
!  
INTEGER :: NSAMP,NMISAMP,NOUT,NDISSOC,NRECOMB,IDISTS,NTSAMP,NBINS,NSNAP  !--isebasti: included NBINS,NSNAP 
INTEGER, DIMENSION(0:100) :: NDISSL  
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: NDSAMPLES  
INTEGER(8), ALLOCATABLE, DIMENSION(:,:) :: NDTRANS,NDROT,NDVIB  
REAL(KIND=8):: TISAMP,XVELS,YVELS,AVDTM,OVTEMP,OVTTEMP,OVRTEMP,OVVTEMP,STARTTIME,TDISTS,STARTIME,DBINV(3),DBINC(3) !--isebasti: DBINV,C included 
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: COLLS,WCOLLS,CLSEP,REAC,SREAC,STEMP,TRANSTEMP,ROTTEMP,VIBTEMP  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: VAR,VARS,CSSS,SUMVIB,CST,PDF,BIN !--isebasti:CST,PDF,BIN included  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CS,VARSP,VIBFRAC,PDFS,BINS !--isebasti:PDFS,BINS included  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: CSS  
!  
!--NSAMP the number of samples  
!--TISAMP the time at which the sampling was last reset  
!--MNISAMP the number of molecules at the last reset  
!--AVDTM the average value of DTM in the cells  
!--NOUT the number of output intervals  
!--COLLS(N) total number of collisions in sampling cell N  
!--WCOLLS(N) total weighted collisins in N  
!--CLSEP(N) sum of collision pair separation in cell N  
!--CS(0,N,L) sampled number of species L in cell N  
!--CS(1,N,L) sampled weighted number of species L in cell N  
!----all the following CS are weighted sums   
!--CS(2,N,L), CS(3,N,L), CS(4,N,L) sampled sum of u, v, w  
!--CS(5,N,L), CS(6,N,L), CS(7,N,L) sampled sum of u*u, v*v, w*w  
!--CS(8,N,L) sampled sum of rotational energy of species L in cell N  
!--CS(8+K,N,L) sampled sum of vibrational energy of species L in cell N  
!              K is the mode  
!  
!--in CSS, M=1 for incident molecules and M=2 for reflected molecules   
!--J=1 for surface at x=XB(1), 2 for surface at x=XB(2)  
!  
!--CSS(0,J,L,M) number sum of molecules of species L  
!--CSS(1,J,L,M) weighted number sum of molecules of species L  
!----all the following CSS are weighted    
!--CSS(2,J,L,M) normal momentum sum to surface   
!--CSS(3,J,L,M) y momentum sum to surface   
!--CSS(4,J,L,M) z momentum sum to surface   
!--CSS(5,J,L,M) tranlational energy sum to surface   
!--CSS(6,J,L,M) rotational energy sum to surface  
!--CSS(7,J,L,M) vibrational energy sum to the surface   
!--CSS(8,J,L,M) reaction energy sum to the surface   
!  
!--CSSS(1,J) weighted sum (over incident AND reflected molecules) of 1/normal vel. component  
!----all the following CSSS are weighted    
!--CSSS(2,J) similar sum of molecular mass / normal vel. component  
!--CSSS(3,J) similar sum of molecular mass * parallel vel. component / normal vel. component  
!--CSSS(4,J) similar sum of molecular mass * speed squared / normal vel. component  
!--CSSS(5,J) similar sum of rotational energy / normal vel. component  
!--CSSS(6,J) similar sum of rotational degrees of freedom /normal velocity component  
!  
!--REAC(N) the number of type N reactions  
!--SREAC(N) the number of type N surface reactions  
!  
!--VAR(M,N) the flowfield properties in cell N  
!--M=1 the x coordinate   
!--M=2 sample size  
!--M=3 number density  
!--M=4 density  
!--M=5 u velocity component  
!--M=6 v velocity component  
!--M=7 w velocity component  
!--M=8 translational temperature  
!--M=9 rotational temperature  
!--M=10 vibrational temperature  
!--M=11 temperature  
!--M=12 Mach number  
!--M=13 molecules per cell  
!--M=14 mean collision time / rate  
!--M=15 mean free path  
!--M=16 ratio (mean collisional separation) / (mean free path)   
!--M=17 flow speed  
!--M=18 scalar pressure nkT  
!--M=19 x component of translational temperature TTX  
!--M=20 y component of translational temperature TTY  
!--M=21 z component of translational temperature TTZ  
!  
!--VARSP(M,N,L) the flowfield properties for species L in cell N  
!--M=0 the sample size  
!--M=1 the fraction  
!--M=2 the temperature component in the x direction  
!--M=3 the temperature component in the y direction  
!--M=4 the temperature component in the z direction  
!--M=5 the translational temperature  
!--M=6 the rotational temperature  
!--M=7 the vibrational temperature  
!--M=8 the temperature  
!--M=9 the x component of the diffusion velocity  
!--M=10 the y component of the diffusion velocity  
!--M=11 the z component of the diffusion velocity  
!   
!--VARS(N,M) surface property N on interval L of surface M  
!  
!--N=1 the incident sample  
!--N=2 the reflected sample  
!--N=3 the incident number flux  
!--N=4 the reflected number flux  
!--N=5 the incident pressure  
!--N=6 the reflected pressure  
!--N=7 the incident parallel shear tress  
!--N=8 the reflected parallel shear stress  
!--N=9 the incident normal-to-plane shear stress  
!--N=10 the reflected normal shear stress  
!--N=11 the incident translational heat flux  
!--N=12 the reflected translational heat fluc  
!--N=13 the incident rotational heat flux  
!--N=14 the reflected rotational heat flux  
!--N=15 the incident vibrational heat flux  
!--N=16 the reflected vibrational heat flux  
!--N=17 the incident heat flux from surface reactions  
!--N=18 the reflected heat flux from surface reactions  
!--N=19 slip velocity      
!--N=20 temperature slip     
!--N=21 rotational temperature slip     
!--N=22 the net pressure  
!--N=23 the net parallel in-plane shear  
!--N=24 the net parallel normal-to-plane shear  
!--N=25 the net translational energy flux  
!--N=26 the net rotational heat flux  
!--N=27 the net vibrational heat flux  
!--N=28 the heat flux from reactions  
!--N=29 total incident heat transfer  
!--N=30 total reflected heat transfer  
!--N=31 net heat transfer   
!--N=32 surface temperature   --not implemented  
!--N=32+K the percentage of species K  
!  
!--COLLS(N) the number of collisions in sampling cell N  
!--WCOLLS(N) weighted number 
!--CLSEP(N) the total collision partner separation distance in sampling cell N  
!  
!--VIBFRAC(L,K,M) the sum of species L mode K in level M  
!--SUMVIB(L,K) the total sample in VIBFRAC  
!  
!--THE following variables apply in the sampling of distribution functions  
!--(some are especially for the dissociation of oxygen  
!  
!--IDISTS 0 before dist sampling, 1 when bulding temperature samples, 2 when building distributions  
!--STARTTIME the time at which sampling starts  
!--TDISTS the nominal time to start distribution function sampling  
!--NDISSOC the number of dissociations  
!--NRECOMB the number of recombinations  
!--NDISSL(L) the number of dissociations from level   
!--NTSAMP the number of temperature samples  
!--OVTEMP,OVTTEMP,OVRTEMP,OVVTEMP the overall temperature, trans. temp, rot. temp,vib, temp   
!--STEMP(L) the temperature of species L  
!--TRANSTEMP(L) the translational temperature of species N  
!--ROTTEMP(L) rotational temperature of species N  
!--VIBTEMP(L) vibrational temperature of species N  
!--NDSAMPLES(L) number in sample for species L  
!--NDTRANS(L,N) number of species L in speed range N  
!--NDROT(L,N) number of species L in rotational energy range N  
!--NDVIB(L,N) number of species L in vibrational level N  
!--STARTIME the time at which the distribution sampling started  
!  
END MODULE OUTPUT
!  
!       ABOVE THIS LINE, JUST COPY MODULES FROM DS1_openmp.f90
!*********************************************************************** 
!*****************************MAIN PROGRAM******************************   
!***********************************************************************  
!  
PROGRAM DS1_POST  
!  
USE MOLECS  
USE GEOM  
USE GAS  
USE CALC  
USE OUTPUT
USE OMP_LIB
USE POST
!  
IMPLICIT NONE  
! 
INTEGER :: N
REAL(KIND=8) :: WCLOCK(5)  
!  
!--set/initialize data
IJSR=123456789
IDISTS=0
!!$ WCLOCK(1)=omp_get_wtime()  
!
!WRITE(*,*) 'Enter output interval to be averaged'
READ(*,*) NOUT
!
!WRITE(*,*) 'Enter case directory'
READ(*,6) DIR
!
!WRITE(*,*) 'Enter output directory'
READ(*,6) DIROUT
!
!WRITE(*,*) 'Enter first run to be considered'
READ(*,*) IRUN
!
!WRITE(*,*) 'Enter last run to be considered'
READ(*,*) FRUN
!
CALL NUMCHAR4 (NOUT,E)
NRUNS=(FRUN-IRUN)+1
!
CALL SET_ZIGGURAT
!
NSS=0 
DO N=IRUN,FRUN
  WRITE(NRUN,5) N
  DIRCASE=TRIM(DIR)//TRIM(NRUN)
  IF(N == IRUN) CALL INIT_AUXVARS
  CALL READ_FILES
  CALL SUM_DATA
END DO
!
CALL AVERAGE_DATA
CALL OUTPUT_RESULTS
!
!!$ WRITE(*,*) 'POST wallclock time (s):    ', omp_get_wtime()-WCLOCK(1)
!
5 FORMAT(I3.3)
6 FORMAT(A)
STOP  
END PROGRAM DS1_POST  
!  
!*****************************************************************************  
!
SUBROUTINE OUTPUT_RESULTS  
!  
!--calculate the surface and flowfield properties  
!--generate TECPLOT files for displaying these properties  
!--calculate collisiion rates and flow transit times and reset time intervals  
!--add molecules to any flow plane molecule output files  
!  
USE MOLECS  
USE CALC  
USE GEOM  
USE GAS  
USE OUTPUT
USE POST  
!   
IMPLICIT NONE  
!  
INTEGER :: IJ,J,JJ,K,L,M,N,NN,NMCR,IDT=0     !--isebasti: included IDT, removed CTIME  
INTEGER(KIND=8) :: NNN   
REAL :: AO,BO,CO,AS,AT  
REAL(KIND=8) :: A,B,C,SDTM,SMCR,DOF,AVW,UU,SVDF,VDOFM,TVIBM,VEL,DTMI,TT,EVIBM,RANF  !--isebasti: included RANF  
REAL(KIND=8), DIMENSION(0:12) :: SUM   
REAL(KIND=8), DIMENSION(0:8,2) :: SUMS  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TVIB,VDOF,PPA  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: TV,THCOL  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: DF  
!INTEGER, ALLOCATABLE, DIMENSION(:) :: NMS   !--isebasti: commented IDT
CHARACTER (LEN=16) :: FILENAME  
!CHARACTER (LEN=4) :: E                      !--isebasti: commented E
!  
!--CTIME  computer time (microseconds)  
!--SUMS(N,L) sum over species of CSS(N,J,L,M) for surface properties  
!   
!--For flowfield properties,where <> indicates sampled sum  
!--SUM(0) the molecular number sum over all species  
!--SUM(1) the weighted number sum over all species  
!--SUM(2) the weighted sum of molecular masses   
!--SUM(3),(4),(5) the weighted sum over species of m*<u>,<v>,<w>  
!--SUM(6) the weighted sum over species of m*(<u**2>+<v**2>+<w**2>)  
!--SUM(7) the weighted sum over species of <u**2>+<v**2>+<w**2>  
!--SUM(8) the weighted sum of rotational energy  
!--SUM(9) the weighted sum of rotational degrees of freedom  
!--SUM(10) the weighted sum over species of m*<u**2>  
!--SUM(11) the weighted sum over species of m*<v**2>  
!--SUM(12) sum over species of m*<w**2>  
!--UU velocity squared  
!--DOF degrees of freedom  
!--AVW the average value of the viscosity-temperature exponent  
!--DVEL velocity difference  
!--TVEL thermal speed  
!--SMCR sum of mcs/mfp over cells  
!--NMCR number in the sum  
!--VDOFM effective vibrational degrees of freedom of mixture  
!--TVIB(L)  
!--VDOF(L)  
!--TV(K,L) the temperature of vibrational mode K of species L  
!--SVDF sum of vibrational degrees of freedom  
!--PPA particles per atom  
!--NMS number per species 
!--NSNAP number of particles considered in SNAPSHOT files 
!  
!---------------------------------------------------------------------------- 
!--allocate/set variables  
!  
ALLOCATE (TV(MMVM,MSP),TVIB(MSP),DF(NCELLS,MMVM,MSP),VDOF(MSP),PPA(MSP),&
          THCOL(MSP,MSP),STAT=ERROR)  !--isebasti: removed NMS(MSP),ISNAP(1,NSNAP),PSNAP(1,NSNAP),VSNAP(3,NSNAP)
IF (ERROR /= 0) THEN  
  WRITE (*,*)'PROGRAM COULD NOT ALLOCATE OUTPUT VARIABLES ',ERROR  
!  STOP  
END IF   
!
IF (NOUT > 9999) NOUT=NOUT-9999
CALL NUMCHAR4 (NOUT,E)    
WRITE (*,*) 'Generating files for output interval',NOUT
!
!----------------------------------------------------------------------------
!--write general/surface properties
!  
IF(ISF == 0) OPEN (3,FILE=TRIM(DIROUT)//'/DS1GEN.DAT')
IF(ISF == 1) OPEN (3,FILE=TRIM(DIROUT)//'/DS1GEN_'//E//'.DAT')
!
IF (IFX == 0) WRITE (3,*) 'DSMC program DS1 for a one-dimensional plane flow'  
IF (IFX == 1) WRITE (3,*) 'DSMC program DS1 for a cylindrical flow'  
IF (IFX == 2) WRITE (3,*) 'DSMC program DS1 for a spherical flow'  
WRITE (3,*)  
!  
WRITE (3,*) 'Interval',NOUT,'Time ',FTIME, ' with',NSAMP,' samples from',TISAMP  
!WRITE (*,*) 'TOTAL MOLECULES = ',NM  
!
!NMS=0  
!DO N=1,NM  
!  M=IPSP(N)  
!  NMS(M)=NMS(M)+1  
!END DO  
WRITE (3,*) 'Total simulated molecules =',NM  
DO N=1,MSP  
! WRITE (*,*) 'SPECIES ',N,' TOTAL = ',NMS(N)  
  WRITE (3,*) 'Species ',N,' total = ',NMS(N)  
END DO  
!
NNN=DINT(TOTMOV)    
WRITE (3,*) 'Total molecule moves   =',NNN  
NNN=DINT(TOTCOL)  
WRITE (3,*) 'Total collision events =',NNN  
NNN=DINT(TOTDUP)  
WRITE (3,*) 'Total duplicate collisions =',NNN  
!  
WRITE (3,*) 'Species dependent collision numbers in current sample'  
DO N=1,MSP  
  WRITE (3,901) TCOL(N,1:MSP)
END DO   
901 FORMAT(20G13.5)   
!  
!CTIME=MCLOCK()
WRITE (3,*) 'Computation time',FLOAT(CTIME)/1000.,'seconds'
WRITE (3,*) 'Collision events per second',(TOTCOL-TOTCOLI)*1000.D00/DFLOAT(CTIME)
WRITE (3,*) 'Molecule moves per second',(TOTMOV-TOTMOVI)*1000.D00/DFLOAT(CTIME)  
WRITE (3,*)  
!  
IF (IDISTS == 2) THEN  
  WRITE (3,*) 'Distribution function sampling started at',STARTIME  
  WRITE (3,*) 'Total dissociations',NDISSOC  
  WRITE (3,*) 'Total recombinations',NRECOMB  
  WRITE (3,*) 'Gas temperature',OVTEMP/DFLOAT(NTSAMP),' K'  
  WRITE (3,*) 'Gas translational temperature',OVTTEMP/DFLOAT(NTSAMP),' K'  
  WRITE (3,*) 'Gas rotational temperature',OVRTEMP/DFLOAT(NTSAMP),' K'  
  WRITE (3,*) 'Gas vibrational temperature',OVVTEMP/DFLOAT(NTSAMP),' K'  
  DO M=1,MSP  
    WRITE (3,*) 'Species',M,' overall temperature',STEMP(M)/DFLOAT(NTSAMP),' K'     
    WRITE (3,*) 'Species',M,' translational temperature',TRANSTEMP(M)/DFLOAT(NTSAMP),' K'         
    WRITE (3,*) 'Species',M,' rotational temperature',ROTTEMP(M)/DFLOAT(NTSAMP),' K'   
    WRITE (3,*) 'Species',M,' vibrational temperature',VIBTEMP(M)/DFLOAT(NTSAMP),' K'   
  END DO  
WRITE (3,*)      
END IF      
!
IF ((ITYPE(1) == 2).OR.(ITYPE(2) == 2)) WRITE (3,*) 'Surface quantities'
DO JJ=1,2  
  IF (ITYPE(JJ) == 2) THEN
    WRITE (3,*)  
    WRITE (3,*) 'Surface at',XB(JJ)  
    WRITE (3,*) 'Incident sample',VARS(0,JJ)   
    WRITE (3,*) 'Number flux',VARS(3,JJ),' /sq m/s'  
    WRITE (3,*) 'Inc pressure',VARS(5,JJ),' Refl pressure',VARS(6,JJ)  
    WRITE (3,*) 'Pressure', VARS(5,JJ)+VARS(6,JJ),' N/sq m'  
    WRITE (3,*) 'Inc y shear',VARS(7,JJ),' Refl y shear',VARS(8,JJ)  
    WRITE (3,*) 'Net y shear',VARS(7,JJ)-VARS(8,JJ),' N/sq m'  
    WRITE (3,*) 'Net z shear',VARS(9,JJ)-VARS(10,JJ),' N/sq m'  
    WRITE (3,*) 'Incident translational heat flux',VARS(11,JJ),' W/sq m'  
    IF (MMRM > 0) WRITE (3,*) 'Incident rotational heat flux',VARS(13,JJ),' W/sq m'  
    IF (MMVM > 0) WRITE (3,*) 'Incident vibrational heat flux',VARS(15,JJ),' W/sq m'  
    WRITE (3,*) 'Total incident heat flux',VARS(29,JJ),' W/sq m'  
    WRITE (3,*) 'Reflected translational heat flux',VARS(12,JJ),' W/sq m'  
    IF (MMRM > 0) WRITE (3,*) 'Reflected rotational heat flux',VARS(14,JJ),' W/sq m'  
    IF (MMVM > 0) WRITE (3,*) 'Reflected vibrational heat flux',VARS(16,JJ),' W/sq m'  
    WRITE (3,*) 'Total reflected heat flux',VARS(30,JJ),' W/sq m'  
    WRITE (3,*) 'Net heat flux',VARS(31,JJ),' W/sq m'  
    WRITE (3,*) 'Slip velocity (y direction)',VARS(19,JJ),' m/s'  
    WRITE (3,*) 'Translational temperature slip',VARS(20,JJ),' K'  
    IF (MMRM > 0) WRITE (3,*) 'Rotational temperature slip',VARS(21,JJ),' K'  
    IF (MSP > 1) THEN  
      DO L=1,MSP  
        WRITE (3,*) 'Species',L,' percentage',VARS(L+32,JJ)  
      END DO  
    END IF
    WRITE (3,*) 
  END IF
END DO
!
PPA=0  
DO N=1,NCELLS  
  DO M=1,MSP  
    PPA(M)=PPA(M)+VARSP(0,N,M)  
  END DO  
END DO
CLOSE(3)
!
!----------------------------------------------------------------------------
!--write number of reactions
!  
!IF(ISF == 0) OPEN (3,FILE=TRIM(DIROUT)//'/DS1REAC.DAT')
!IF(ISF == 1) OPEN (3,FILE=TRIM(DIROUT)//'/DS1REAC_'//E//'.DAT')
OPEN (3,FILE=TRIM(DIROUT)//'/DS1REAC.DAT',ACCESS='APPEND')
!
IF ((JCD == 1).AND.(MSP > 1)) THEN  
  WRITE (3,*) 'New integrated DSMC chemistry model with no experimentally based rates'  
  WRITE (3,*)  
  WRITE(3,*) 'GAINS FROM REACTIONS'  
  WRITE(3,*) '                          Dissoc.     Recomb. Endo. Exch.  Exo. Exch.'  
  DO M=1,MSP  
  WRITE (3,*) ' SPECIES',M,TREACG(1,M),TREACG(2,M),TREACG(3,M),TREACG(4,M)     
  END DO  
  WRITE (3,*)   
  WRITE(3,*) 'LOSSES FROM REACTIONS'  
  WRITE(3,*) '                          Dissoc.     Recomb. Endo. Exch.  Exo. Exch.'  
  DO M=1,MSP  
  WRITE (3,*) ' SPECIES',M,TREACL(1,M),TREACL(2,M),TREACL(3,M),TREACL(4,M)  
  END DO  
  WRITE (3,*)  
  WRITE (3,*) 'TOTALS'  
  DO M=1,MSP  
    WRITE (3,*) ' SPECIES',M,' GAINS',TREACG(1,M)+TREACG(2,M)+TREACG(3,M)+TREACG(4,M),&  
                ' LOSSES',TREACL(1,M)+TREACL(2,M)+TREACL(3,M)+TREACL(4,M)    
  END DO  
END IF
!  
IF ((JCD == 0).AND.(MSP > 1)) THEN  
!  WRITE (3,994) '# Traditional DSMC chemistry model based on rate equations'
!  WRITE (3,994) '# Macro. or Coll. Temps: ITCV, IEAA, IZV',ITCV,IEAA,IZV  
!  WRITE (3,994) '# MNRE:',MNRE
!  WRITE (3,994) '# Reaction    Number    Fraction'
  AS=0.
  DO M=1,MNRE
    AS=AS+REAC(M)
  END DO  
  IF (AS == 0.) AS=1.d0
!  DO M=1,MNRE  
    WRITE (3,993) FTIME*1.d6,REAC(1:MNRE)/AS
!  END DO   
END IF
CLOSE(3) 
!
!----------------------------------------------------------------------------
!--write flowfield overall properties
!  
IF(ISF == 0) OPEN (3,FILE=TRIM(DIROUT)//'/DS1FP00.DAT')
IF(ISF == 1) OPEN (3,FILE=TRIM(DIROUT)//'/DS1FP00_'//E//'.DAT')
!
WRITE (3,994) '# Flowfield Overall Properties '
WRITE (3,994) '# NCELLS:', NCELLS  
WRITE (3,994) '# NSAMP: ', NSAMP
WRITE (3,994) '# X-coord.     Cell   Sample     Number Dens.   Density   u velocity &  
            &  v velocity   w velocity   Trans. Temp.   Rot. Temp.   Vib. Temp. &  
            &  Temperature  Mach no.     Mols/cell    m.c.t        m.f.p     &  
            &  mcs/mfp        speed      Pressure        TTX          TTY          TTZ     &
            &  dtm/mct     <dx/mfp>      Fx          Fy            Fz         Qtransfer &
            &  Species Fractions'  
DO N=1,NCELLS  
  WRITE (3,996) VAR(1,N),N,VAR(2:21,N),DTM/VAR(14,N),(CELL(3,N)-CELL(2,N))/DFLOAT(NCIS)/VAR(15,N),&
                CST(1:4,N)/CST(0,N),0.01*VARSP(1,N,1:MSP)
END DO
CLOSE(3)
!
!----------------------------------------------------------------------------
!--write flowfield properties per species
!
DO L=1,MSP  
    WRITE(FILENAME,777) L
777 FORMAT('/DS1FP',i2.2)
IF(ISF == 0) OPEN (3,FILE=TRIM(DIROUT)//TRIM(FILENAME)//'.DAT')
IF(ISF == 1) OPEN (3,FILE=TRIM(DIROUT)//TRIM(FILENAME)//'_'//E//'.DAT')
!
  WRITE (3,994) '# Flowfield Properties Species:',L  
  WRITE (3,994) '# NCELLS:', NCELLS  
  WRITE (3,994) '# NSAMP: ', NSAMP
  WRITE (3,994) '# X-coord.     Cell   Sample    Percentage   Species TTx   Species TTy  Species TTz &  
             & Trans. Temp.  Rot. Temp.   Vib. Temp.   Spec. Temp  u Diff Vel   v Diff Vel   w Diff Vel.'   
  DO N=1,NCELLS  
    WRITE (3,996) VAR(1,N),N,VARSP(0,N,L),0.01*VARSP(1,N,L),VARSP(2:11,N,L)
  END DO
CLOSE(3)  
END DO
!
!----------------------------------------------------------------------------    
!--write composition of a reacting gas as a function of time  
!
OPEN (10,FILE=TRIM(DIROUT)//'/COMPOSITION.DAT',ACCESS='APPEND')  
  AS=NM  
  WRITE (10,993) FTIME*1.d6,NM,UVFX(:),VAR(11,1),NOUT,NMS(1:MSP)!/AS 
CLOSE (10)  
!
!----------------------------------------------------------------------------
!--write SAMPLE_PDF files  
!
IF(IPDF == 1) THEN  
!  
IF(ISF == 0) OPEN (7,FILE=TRIM(DIROUT)//'/PDFUVW.DAT',FORM='FORMATTED')
IF(ISF == 1) OPEN (7,FILE=TRIM(DIROUT)//'/PDFUVW_'//E//'.DAT',FORM='FORMATTED')    
   
WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
WRITE (7,994) '# Interval, Equilibrium pdf, U V W pdfs total, U V W pdfs spec1, U V W pdfs spec2 ... '
DO N=1,NBINS
  A=DBINV(1)+DBINV(3)*(DFLOAT(N)-0.5d0)                !normalized U/VMP
  B=(1.d0/SPI)/DEXP(A*A)                               !Boltzmann distribution
  WRITE (7,993) A,B,PDF(N,1:3),(PDFS(N,1:3,L),L=1,MSP) !pdfs
END DO  
CLOSE(7) 
!
OPEN (7,FILE='PDFC.DAT',FORM='FORMATTED')  
WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
WRITE (7,994) '# Interval, Equilibrium pdf, C pdf total, C pdf spec1, C pdf spec2 ... '
DO N=1,NBINS
  A=DBINC(1)+DBINC(3)*(DFLOAT(N)-0.5d0)                !normalized C/VMP
  B=((4.d0/SPI)*A*A)/DEXP(A*A)                         !Boltzmann distribution 
  WRITE (7,993) A,B,PDF(N,4),(PDFS(N,4,L),L=1,MSP)     !pdfs
END DO  
CLOSE(7)
!
END IF
!
!----------------------------------------------------------------------------
!--write SNAPSHOT files  
!  
OPEN (7,FILE=TRIM(DIROUT)//'/SNAPSHOT_'//E//'.xyz',FORM='FORMATTED')  
WRITE (7,993) NSS
WRITE (7,994) '# X Y Z Vx Species Radius - at time', FTIME*1.d6
DO N=1,NSS  
  CALL ZGF(RANF,IDT)
  A=0.2d0*RANF
  CALL ZGF(RANF,IDT)
  B=0.2d0*RANF
  L=ISNAP(1,N)
  WRITE (7,993) PSNAP(1,N)/(XB(2)-XB(1)),A,B,VSNAP(1,N),L,SP(1,L)*1.d7
END DO  
CLOSE(7)
!
!----------------------------------------------------------------------------
!--write vibrational distribution
!IF (MMVM > 0) THEN  
!  CALL CHECK_VIB_DIST   
!  IF (JCD == 1) CALL CHECK_DISS_DIST  
!END IF   
!  
!----------------------------------------------------------------------------  
!--write a special output file for vibrational temperature and temperature
!  versus collision number
!
!IF (ISF == 1) THEN
!  OPEN (10,FILE='RELAX.DAT',ACCESS='APPEND')  
!  AO=2.*TOTCOL/NM  
!  BO=0.  
!  CO=0.  
!  DO N=1,NCELLS  
!    BO=BO+VAR(11,N)  
!    CO=CO+VAR(10,N)  
!  END DO    
!  WRITE (10,*) AO,BO/NCELLS,CO/NCELLS  
!  CLOSE (10)  
!END IF
!
!----------------------------------------------------------------------------
!--I/O format
!
999 FORMAT (I5,50G13.5)  
998 FORMAT (G270.0)  
997 FORMAT (G175.0)
996 FORMAT (G13.5,I5,50G13.5)  
995 FORMAT (22G13.5)
994 FORMAT (A,22G13.5)
993 FORMAT (99G13.5)
!
!----------------------------------------------------------------------------
!--deallocate local variables  
!  
DEALLOCATE (TV,TVIB,DF,VDOF,THCOL,STAT=ERROR)  
IF (ERROR /= 0) THEN  
  WRITE (*,*)'PROGRAM COULD NOT DEALLOCATE OUTPUT VARIABLES',ERROR  
END IF  
!
!----------------------------------------------------------------------------
!  
RETURN  
!  
END SUBROUTINE OUTPUT_RESULTS   
!   
!**************************************************************************************
!  
SUBROUTINE INIT_AUXVARS  
!  
USE MOLECS  
USE GEOM  
USE GAS  
USE CALC  
USE OUTPUT
USE POST  
!  
IMPLICIT NONE  
!  
!
101 CONTINUE  
OPEN (7,FILE=TRIM(DIRCASE)//'/PARAMETERS.DAT',FORM='UNFORMATTED',ERR=101)
READ (7) NCCELLS,NCELLS,MMRM,MMVM,MNM,MNRE,MNSR,MSP,MTBP,ILEVEL,MDIV,IRECOM,MMEX,MEX,ISF,NBINS,NSNAP !--isebasti: included ISF,NBINS,NSAP  
CLOSE(7)  
!  
IF (MMVM > 0) THEN  
  ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),      &  
       IPVIB(MMVM,MNM),STAT=ERROR)  
ELSE  
  IF (MMRM > 0) THEN  
    ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),STAT=ERROR)  
  ELSE  
    ALLOCATE (PX(MNM),PTIM(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),STAT=ERROR)  
  END IF    
END IF  
IF (ERROR /= 0) THEN  
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR MOLECULE ARRAYS',ERROR  
ENDIF  
!  
ALLOCATE (JDIV(0:ILEVEL,MDIV),STAT=ERROR)  
IF (ERROR /= 0) THEN  
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR JDIV ARRAY',ERROR  
ENDIF  
!  
ALLOCATE (CELL(4,NCELLS),ICELL(NCELLS),CCELL(5,NCCELLS),ICCELL(3,NCCELLS),STAT=ERROR)  
IF (ERROR /= 0) THEN  
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR CELL ARRAYS',ERROR  
ENDIF  
!  
ALLOCATE (COLLS(NCELLS),WCOLLS(NCELLS),CLSEP(NCELLS),REAC(MNRE),SREAC(MNSR),VAR(21,NCELLS), &  
          VARSP(0:11,NCELLS,MSP),VARS(0:32+MSP,2),CS(0:8+MMVM,NCELLS,MSP),CSS(0:8,2,MSP,2), &
          CSSS(6,2),CST(0:4,NCELLS),BINS(0:NBINS,4,MSP),BIN(0:NBINS,4),&
          PDFS(0:NBINS,4,MSP),PDF(0:NBINS,4),STAT=ERROR) !--isebasti: CST,BINS,BIN.PDFS,PDF included 
IF (ERROR /= 0) THEN  
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR SAMPLING ARRAYS',ERROR  
ENDIF  
!  
IF (MMVM > 0) THEN  
  ALLOCATE (VIBFRAC(MSP,MMVM,0:150),SUMVIB(MSP,MMVM),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR RECOMBINATION ARRAYS',ERROR  
  ENDIF  
END IF  
!  
CALL ALLOCATE_GAS 
!
!--ABOVE THIS LINE EVERYTHING IS THE SAME AS IN DS1_openmp.f90
!
IF (IDISTS > 0) THEN    
  ALLOCATE (STEMP(MSP),TRANSTEMP(MSP),ROTTEMP(MSP),VIBTEMP(MSP),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE DISTRIBUTION FUNCTION VARIABLES',ERROR  
  END IF
END IF
!  
ALLOCATE (NMS(MSP),ISNAP(1,NSNAP),PSNAP(1,NSNAP),VSNAP(3,NSNAP),STAT=ERROR)  
IF (ERROR /= 0) THEN  
  WRITE (*,*)'PROGRAM COULD NOT ALLOCATE OUTPUT VARIABLES ',ERROR  
END IF  
!
!--allocate/initiate auxiliar variables (after reading the first file)
!
ALLOCATE (IVNMS(MSP),IVTREACG(4,MSP),IVTREACL(4,MSP),DVREAC(MNRE),DVTCOL(MSP,MSP),&
          DVSTEMP(MSP),DVTRANSTEMP(MSP),DVROTTEMP(MSP),DVVIBTEMP(MSP),&
          DVVARS(0:32+MSP,2),DVVAR(21,NCELLS),DVCST(0:4,NCELLS),DVPDF(0:NBINS,4),DVBIN(0:NBINS,4),&
          IVIPSP(1,NSNAP),DVPX(1,NSNAP),DVPV(3,NSNAP),DVVARSP(0:11,NCELLS,MSP),&
          DVPDFS(0:NBINS,4,MSP),DVBINS(0:NBINS,4,MSP),&
          STAT=ERROR)  
IF (ERROR /= 0) THEN  
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE AUXILIAR POST VARIABLES',ERROR  
ENDIF
!
IVNSAMP=0.
IVNM=0.
IVCTIME=0.
IVTREACG=0.
IVTREACL=0.
IVNMS=0.
!
DVTISAMP=0.
DVFTIME=0.	
DVTOTMOV=0.
DVTOTCOL=0.
DVTOTDUP=0.
DVTCOL=0.
DVTOTCOLI=0.
DVTOTMOVI=0.
IF (IDISTS > 0) THEN
  IVNTSAMP=0.
  IVNDISSOC=0.
  IVNRECOMB=0.
  DVSTARTIME=0.
  DVOVTEMP=0.
  DVOVTTEMP=0.
  DVOVRTEMP=0.
  DVOVVTEMP=0.
  DVSTEMP=0.
  DVTRANSTEMP=0.
  DVROTTEMP=0.
  DVVIBTEMP=0.
END IF
DVVARS=0.
DVVARSP=0.
DVREAC=0.
DVVAR=0.
DVDTM=0.
DVCST=0.
DVUVFX=0.
DVBIN=0.
DVBINS=0.
DVPDF=0.
DVPDFS=0.
!
RETURN  
END SUBROUTINE INIT_AUXVARS  
!  
!**************************************************************************************
!  
SUBROUTINE READ_FILES  
!  
USE MOLECS  
USE GEOM  
USE GAS  
USE CALC  
USE OUTPUT
USE POST  
!  
IMPLICIT NONE  
!  
INTEGER :: ZCHECK
!  
102 CONTINUE  
OPEN (7,FILE=TRIM(DIRCASE)//'/DS1OUT_'//E//'.BIN',FORM='UNFORMATTED',ERR=102)  
READ (7)IFX,NOUT,FTIME,NSAMP,TISAMP,NM,TOTMOV,TOTCOL,TOTDUP,TCOL,&
        CTIME,TOTCOLI,TOTMOVI,IDISTS,STARTIME,NDISSOC,NRECOMB,OVTEMP,NTSAMP,OVTTEMP,OVRTEMP,&
        OVVTEMP,STEMP,TRANSTEMP,ROTTEMP,VIBTEMP,ITYPE,XB,VARS,NCELLS,VARSP,&  !end of general properties
        JCD,TREACG,TREACL,ITCV,IEAA,IZV,REAC,VAR,DTM,CELL,NCIS,CST,UVFX,NMS,& !end of reactions,flow properties, and composition
        IPDF,BIN,BINS,DBINV,SPI,PDF,PDFS,DBINC,ISNAP,PSNAP,VSNAP,SP,ZCHECK    !end of sample pdf and snapshot
CLOSE(7) 
!  
IF (ZCHECK /= 1234567) THEN  
  WRITE (*,*) NM,' Molecules, Check integer =',ZCHECK  
  STOP  
ELSE  
!  WRITE (*,*) TRIM(DIRCASE)//'/DS1OUT_'//E//'.BIN file read, Check integer=',ZCHECK    
END IF 
!  
RETURN  
END SUBROUTINE READ_FILES  
!  
!**************************************************************************************
!  
SUBROUTINE SUM_DATA  
!  
USE MOLECS  
USE GEOM  
USE GAS  
USE CALC  
USE OUTPUT
USE POST  
!  
IMPLICIT NONE
!
INTEGER :: K,L,M,N,IDT=0
REAL(8) :: RANF  
!  
!
IVNSAMP=IVNSAMP+NSAMP
IVNM=IVNM+NM
IVCTIME=IVCTIME+CTIME
IVTREACG=IVTREACG+TREACG
IVTREACL=IVTREACL+TREACL
IVNMS=IVNMS+NMS
!
DVTISAMP=DVTISAMP+TISAMP
DVFTIME=DVFTIME+FTIME	
DVTOTMOV=DVTOTMOV+TOTMOV
DVTOTCOL=DVTOTCOL+TOTCOL
DVTOTDUP=DVTOTDUP+TOTDUP
DVTCOL=DVTCOL+TCOL
DVTOTCOLI=DVTOTCOLI+TOTCOLI
DVTOTMOVI=DVTOTMOVI+TOTMOVI
IF (IDISTS > 0) THEN
  IVNTSAMP=IVNTSAMP+NTSAMP
  IVNDISSOC=IVNDISSOC+NDISSOC
  IVNRECOMB=IVNRECOMB+NRECOMB
  DVSTARTIME=DVSTARTIME+STARTIME
  DVOVTEMP=DVOVTEMP+OVTEMP
  DVOVTTEMP=DVOVTTEMP+OVTTEMP
  DVOVRTEMP=DVOVRTEMP+OVRTEMP
  DVOVVTEMP=DVOVVTEMP+OVVTEMP
  DVSTEMP=DVSTEMP+STEMP
  DVTRANSTEMP=DVTRANSTEMP+TRANSTEMP
  DVROTTEMP=DVROTTEMP+ROTTEMP
  DVVIBTEMP=DVVIBTEMP+VIBTEMP
END IF
DVVARS=DVVARS+VARS
DVVARSP=DVVARSP+VARSP
DVREAC=DVREAC+REAC
DVVAR=DVVAR+VAR
DVDTM=DVDTM+DTM
DVCST=DVCST+CST
DVUVFX=DVUVFX+UVFX
DVBIN=DVBIN+BIN
DVBINS=DVBINS+BINS
DVPDF=DVPDF+PDF
DVPDFS=DVPDFS+PDFS
!
!--create vectors by sampling particles from all NRUNS
!
M=NSNAP/DFLOAT(NRUNS)
DO K=1,M
  CALL ZGF(RANF,IDT)
  N=INT(RANF*DFLOAT(NSNAP))+1
  NSS=NSS+1
  IVIPSP(1,NSS)=ISNAP(1,N)
  DVPV(:,NSS)=VSNAP(:,N)  
  DVPX(1,NSS)=PSNAP(1,N)
END DO
!  
RETURN  
END SUBROUTINE SUM_DATA  
!  
!**************************************************************************************
!  
SUBROUTINE AVERAGE_DATA  
!  
USE MOLECS  
USE GEOM  
USE GAS  
USE CALC  
USE OUTPUT
USE POST  
!  
IMPLICIT NONE
!
REAL(8) :: A
!  
!
A=DFLOAT(NRUNS)
!
NSAMP=IVNSAMP/A
NM=IVNM/A
CTIME=IVCTIME/A
TREACG=IVTREACG/A
IVTREACL=IVTREACL/A
IVNMS=IVNMS/A
!
TISAMP=DVTISAMP/A
FTIME=DVFTIME/A	
TOTMOV=DVTOTMOV/A
TOTCOL=DVTOTCOL/A
TOTDUP=DVTOTDUP/A
TCOL=DVTCOL/A
TOTCOLI=DVTOTCOLI/A
TOTMOVI=DVTOTMOVI/A
IF (IDISTS > 0) THEN
  NTSAMP=IVNTSAMP/A
  STARTIME=DVSTARTIME/A
  NDISSOC=IVNDISSOC/A
  NRECOMB=IVNRECOMB/A
  OVTEMP=DVOVTEMP/A
  OVTTEMP=DVOVTTEMP/A
  OVRTEMP=DVOVRTEMP/A
  OVVTEMP=DVOVVTEMP/A
  STEMP=DVSTEMP/A
  TRANSTEMP=DVTRANSTEMP/A
  ROTTEMP=DVROTTEMP/A
  VIBTEMP=DVVIBTEMP/A
END IF
VARS=DVVARS/A
VARSP=DVVARSP/A
REAC=DVREAC/A
VAR=DVVAR/A
DTM=DVDTM/A
CST=DVCST/A
UVFX=DVUVFX/A
BIN=DVBIN/A
BINS=DVBINS/A
PDF=DVPDF/A
PDFS=DVPDFS/A
!
!--rewrite original variables using sampled data
!
ISNAP=IVIPSP
VSNAP=DVPV
PSNAP=DVPX
!  
RETURN  
END SUBROUTINE AVERAGE_DATA  
!  
!***************************************************************************
!  
SUBROUTINE ALLOCATE_GAS  
!  
USE GAS  
USE CALC  
!  
IMPLICIT NONE  
!  
ALLOCATE (FSP(MSP,2),SP(5,MSP),SPR(3,MSP),SPM(8,MSP,MSP),ISPR(2,MSP),ISPV(MSP),ENTR(6,MSP,2),NRSP(MSP,MSP),      &  
          VMP(MSP,2),UVMP(MSP,2),VNMAX(MSP),CR(MSP),TCOL(MSP,MSP),ISPRC(MSP,MSP),SPRC(2,MSP,MSP,MSP),STAT=ERROR)  !--isebasti: UVMP included
!  
IF (ERROR /= 0) THEN  
  WRITE (*,*)'PROGRAM COULD NOT ALLOCATE SPECIES VARIABLES',ERROR  
END IF  
!  
ALLOCATE (NEX(MMEX,MSP,MSP),NSPEX(MSP,MSP),SPEX(3,MMEX,MSP,MSP),ISPEX(MMEX,0:4,MSP,MSP),TREACG(4,MSP),         &  
          PSF(MMEX),TREACL(4,MSP),TNEX(MEX),STAT=ERROR)  
!  
IF (ERROR /= 0) THEN  
  WRITE (*,*)'PROGRAM COULD NOT ALLOCATE Q-K REACTION VARIABLES',ERROR  
END IF  
!  
IF (MMVM > 0) THEN  
  
  ALLOCATE (SPVM(4,MMVM,MSP),ISPVM(2,MMVM,MSP),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE VIBRATION VARIABLES',ERROR  
  END IF  
END IF  
!  
IF (MNRE > 0) THEN  
  ALLOCATE (CI(MNRE),AE(MNRE),AC(MNRE),BC(MNRE),ER(MNRE),REA(6,MNRE),LE(MNRE),ME(MNRE),KP(MNRE),LP(MNRE),MP(MNRE),      &  
            IREA(2,MNRE),JREA(2,MNRE,2),IRCD(MNRE,MSP,MSP),NREA(2,MNRE),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE REACTION VARIABLES',ERROR  
  END IF  
END IF   
!  
IF (MTBP > 0) THEN  
  ALLOCATE (THBP(MTBP,MSP),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE THIRD BODY VARIABLES',ERROR  
  END IF  
END IF  
!  
IF (MNSR > 0) THEN  
  ALLOCATE (ERS(MNSR),LIS(2,MNSR),LRS(6,MNSR),ISRCD(MNSR,MSP),STAT=ERROR)  
  IF (ERROR /= 0) THEN  
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE SURFACE REACTION VARIABLES',ERROR  
  END IF  
END IF  
!  
RETURN  
END SUBROUTINE ALLOCATE_GAS  
!  
!*****************************************************************************  
!  
SUBROUTINE NUMCHAR4(NNN,E)  
!  
!--produces the character equivalent E of a 4 digit integer NNN   
!  
CHARACTER(LEN=1) :: A  
CHARACTER(LEN=1) :: B  
CHARACTER(LEN=1) :: C  
CHARACTER(LEN=1) :: D  
CHARACTER(LEN=4) :: E  
A='0' ; B='0' ; C='0' ; D='0'  
N=NNN  
IF (N.GT.999) THEN  
  L=N/1000  
  A=CHAR(48+L)  
  N=N-1000*L  
END IF  
IF (N.GT.99) THEN  
  L=N/100  
  B=CHAR(48+L)  
  N=N-100*L  
END IF  
IF (N.GT.9) THEN  
  L=N/10  
  C=CHAR(48+L)  
  N=N-10*L  
END IF  
D=CHAR(48+N)  
E=A//B//C//D  
!  
RETURN  
!  
END SUBROUTINE NUMCHAR4  
!   
!*****************************************************************************  
!
SUBROUTINE SET_ZIGGURAT
!
!--author: isebasti
!--calculations required by Ziggurate Method
!--call this subroutine with a different jsrl value for a different set of seeds
!
USE CALC
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: nth,idtl,ishr3  !local variables
!
nth=1                      !for serial compilation
!$omp parallel
!$omp master
!$    nth=omp_get_num_threads( )
!!$    write(*,*) 'Available threads     ', omp_get_max_threads( )
!!$    write(*,*) 'Threads in use        ', nth
!$omp end master
!$omp end parallel
!
allocate(iseed(0:nth-1))
!
do idtl=0,nth-1
call ishr(ishr3)
iseed(idtl)=ishr3         !set a seed for each thread
!!$    write(*,*) 'Thread',idtl,'Seed',iseed(idtl)
end do
!
RETURN
END SUBROUTINE SET_ZIGGURAT
!
!*****************************************************************************
!
SUBROUTINE ISHR(ISHR3L)
!
!--author: isebasti
!--calculations required by Ziggurate Method; initialize different seeds
!
USE CALC
!
IMPLICIT NONE
!
INTEGER :: jsrl,jsr_inputl,ishr3l  !local variables
!
jsrl=ijsr
jsr_inputl = jsrl
jsrl = ieor(jsrl,ishft(jsrl,13))
jsrl = ieor(jsrl,ishft(jsrl,-17))
jsrl = ieor(jsrl,ishft(jsrl,  5))
ishr3l = jsr_inputl+jsrl
ijsr=jsrl
!
RETURN
END SUBROUTINE ISHR
!
!*****************************************************************************
!
SUBROUTINE ZGF(RANFL,IDTL)
!
!--author: isebasti
!--generate a random fraction ]0,1[ using Ziggurate Method (Marsaglia & Tsang)
!--this openmp implemation is based on code available at
!--people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat_openmp/ziggurat_openmp.html
!
USE CALC
!
IMPLICIT NONE
!
REAL(KIND=8) :: ranfl
INTEGER :: jsrl,idtl,jsr_inputl   !local variables
!
jsrl=iseed(idtl)
jsr_inputl = jsrl
jsrl = ieor(jsrl,ishft(jsrl, 13))
jsrl = ieor(jsrl,ishft(jsrl,-17))
jsrl = ieor(jsrl,ishft(jsrl,  5))
ranfl = 0.5e0+0.2328306e-9*real(jsr_inputl+jsrl)
iseed(idtl) = jsrl
!
RETURN
END SUBROUTINE ZGF
!
!*****************************************************************************
!
