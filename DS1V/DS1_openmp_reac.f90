!******************************************************************************
!****DSMC program for one-dimensional plane, cylindrical or spherical flows****
!******************************************************************************
!
!--Version 1.22
!
!******************************************************************************
!
!--Development History
!
!--Version 1.1  writing commenced 16 November 2007, and was used for the RGD26 chemistry work
!--Version 1.2  added the full coding for cylindrical and spherical flows (from April 2009)
!--Version 1.3  added the moving wall and enhanced unsteady sampling options
!--Version 1.4  added radial weighting factors
!--Version 1.5  extended the new chemistry model to cope with O2-H2 sysyem
!--Version 1.6  modified model for exothermic reactions
!--Version 1.7  added velocity in y direction
!--Version 1.8  removed diagnostic subroutines - first version for web (POSTED Nov 25, 2010)
!--Version 1.9  revised volume for ternary collisions, only the relevant vibrational modes are considered in reactions,
!               an improved procedure when there are multiple reactions possible for the collision pair
!--Version 1.10 new Q-K procedure for exothermic exchange and chain reactions
!--Version 1.11 includes the newer Q-K procedures as described in the Phys Fluids paper
!--Version 1.12 gives priority to dissociation when considering reactions
!--Version 1.13 option to base, ternary collision volume, activation energy adjustments, and
!               vibrational relaxation numbers on either the macroscopic or collision temperatures
!--Version 1.14 includes the sampling of distribution functions
!--Version 1.15 TCE model is set as standard case, bug in reflected particles is fixed --isebasti (Nov 2013)
!--Version 1.16 added thermophoretic calculations --isebasti (Nov 2013)
!--Version 1.17 openmp version + diffusion flame set up --isebasti (May 2015)
!--Version 1.18 added Marat's QCT-VVT relaxation model for O2-O collisions --isebasti (Oct 2015)
!--Version 1.19 integrated both TCE+LB and QCT-based (Marat's) model in the same code (Mar 2016)
!--Version 1.20 special coding for shockwave sampling (Apr 2016)
!--Version 1.21 adding calibrated O2-O2 kinetic parameters (May 2016)
!--Version 1.22 adding Han's MEQCTVT model for N2-O collisions --isebasti (Sep 2016)
!******************************************************************************
!
!--Base SI units are employed throughout the program
!
!--Variables in data file DS1VD.DAT
!
!--NVER the n in version n.m that generated the data file
!--MVER the m in version n.m that generated the data file
!
!--m must change whenever the data file is altered
!
!--AMEG the approximate number of megabytes that are initially available to the calculation
!
!--IFX 0 for plane flows, 1 for cylindrical flows or 2 for spherical flows (the axis for cyl. flows and the center of spherical flows is at x=0)
!
!--the type of the boundaries of the flowfield is set by ITYPE, the number code is
!    0 for a boundary with the stream
!    1 plane of symmetry for IFX=0, or center for IFX>0 and XB(1)=0, or mandatory for outer boundary if IVB=1
!    2 for solid surface
!    3 for vacuum
!
!--XB(1) the minimum x coordinate of the flow (must be zero or positive if IFX>0)
!--ITYPE(1) sets type of XB(1)
!
!--XB(2) the maximum x coordinate of the flow
!--ITYPE(2) sets type of XB(2)
!
!--if IFX>0
!----IWF 0 for no weighting factors, 1 for weighting factors based on radius
!------the weighting factors are 1+WFM*(r/r_max)**IFX
!----if (IWF=1) WFM the (maximum weighting factor -1)
!--end if
!
!--GASCODE the code number of the gas
!    1 for a hard sphere gas
!    2 for argon
!    3 for ideal gas
!    4 for real oxygen
!    5 for ideal air
!    6 for real air @ 7.5 km/s
!    7 for helium-xenon mixture with argon molecular weight
!    8 for combustible hydrogen-oxygen
!
!--IGS 0 if the initial state is a vacuum, 1 if it is the stream or reference gas and (if present) the secondary stream
!      2 if it is the stream or reference gas with TDF temperature distribution
!--ISECS=0 if there is no secondary stream, 1 if there is a secondary stream that applies when x>0. for IFX=0
!          (XB(1) must then be negative and XB(2) must be positive),or x>XS for IFX>0.
!--if (IFX>0 and ISECS=1)
!----XS the initial boundary between the initial and secondary streams
!----ISECS must be zero if there are radial weighting factors
!--in a loop over K=1,2 (1 for XB(1), 2 for XB(2)
!----if K=1 or (K=2 and ISECS=1)
!------FND(K)  the number density of the stream or reference gas
!------FTMP(K) temperature of stream or reference gas
!------FRTMP(K) rotational temperature of stream or reference gas (if any species have rotational modes)
!------FVTMP(K) vibrational temperature of stream or reference gas (if any species have vibrational modes)
!------VFX(K) the velocity component in the x direction
!------VFY(K) the velocity component in the y direction
!------in a loop over the number of species (if there is more than one species)
!--------FSP(N,K)) the fraction of species N in the stream or reference gas
!------end if
!----end if
!----if ITYPE(K)=2
!------TSURF(K) the surface temperature
!------FSPEC(K) the fraction of specular reflection at this surface
!------VSURF(K) the y velocity component
!----end if
!--end loop over ends
!
!--if ITYPE(1)=0 (stream) and IFX=0
!----IREM 0 for no molecule removal (based on "stagnation streamline" theory)
!         1 for steady removal (FREM=1) set as data (e.g. for shock wave problem)
!         2 for removal set as a restart option
!----if IREM=1
!------XREM the coordinate at which the removal starts (uniform from XREM to XB(2))
!----end if
!--end if
!
!--if ITYPE(2)=1
!----IVB 0 if the maximum x boundary is stationary, 1 if it is a moving specularly reflecting piston
!------if IVB=1, ITYPE(2)=1
!----if (IVB = 1) VELOB  the velocity of the boundary
!--end if
!
!--MOLSC target number of molecules per sampling cell
!--NCIS number of collision cells in a sampling cell
!--CPDTM collisions in time step
!--collisions appropriate to CPDTM*(local mean collision time)are calculated in a collision cell
!  when its time falls half this time interval behind the flow time
!--TPDTM maximum flow transits of a sampling cell in a time step
!--a molecule is moved appropriate to the minimum of CPDTM*(local mean coll time) and TPDTM*(local transit time)
!  when its time falls half this time interval behind the flow time
!--NNC 0 to select collision partner randomly from collision cell, 1 for nearest neighbor collisions
!--IMTS 0 for move time steps that are the same for all cells, 1 if they are allowed to vary over the cells
!--SAMPRAT the number of time steps between samplings
!--OUTRAT the number of sampling steps bewteen the output of results
!--ISF 0 if the sample is reset only as a user restart option, 1 if it is reset at the end of each output interval (unsteady flow sampling)
!--FRACSAM fraction (at the end) of the output interval in which samples are made in a flow with unsteady sampling;
!  sample for implicit boudary conditions is reset every FRACSAM*OUTRAT samples; FRACSAM possible values=0.5,0.25,0.2,0.1,0.05,0.02,0.025,0.01...
!--JCD 0 if any reactions are based on rate equations (GASCODE=6 or 8 only), 1 if they are based on the Q-K reaction model
!--IRECOM 0 if recombination reactions can be neglected, 1 if they must be included
!
!************************************************************************************
!
!--Restart options in Version 1.1
!
!--if IREM=2
!----FREM fraction of molecule removal
!----XREM the coordinate at which removal commences
!--end if
!
!************************************************************************************
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
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: PV,VIBEN
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: PX,PTIM,PROT
!
!--PX(N) position (x coordinate) of molecule N
!--PTIM(N) molecule time
!--PROT(N) rotational energy
!--PV(1-3,N) u,v,w velocity components
!--VIBEN(IV,L) vibrational energy of level IV for species L
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
REAL(KIND=8), DIMENSION(2) :: FND,FTMP,FRTMP,FVTMP,VFX,VFY,TSURF,FSPEC,VSURF,UVFX,UFND,UFTMP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CI,AE,AC,BC,ER,ERS,CR,TNEX,PSF,DPER,FPTEMP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: FSP,SP,SPR,SPV,REA,THBP,VMP,UVMP,CPER
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: SPM,SPVM,ENTR
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: SPEX,SPRC,FEVIB
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FPVIB
INTEGER :: MSP,MMVM,MMRM,MNRE,MNSR,MTBP,GASCODE,MMEX,MEX,ITMAX
INTEGER, ALLOCATABLE, DIMENSION(:) :: ISP,ISPV,LE,ME,KP,LP,MP,IREV
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISPR,IREA,NREA,NRSP,LIS,LRS,ISRCD,ISPRC,TREACG,TREACL,NSPEX,IDL
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ISPVM,NEX,IRCD,JREA,NEVIB
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ISPEX
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: NPVIB  !--isebasti: UVFX,UFND,UFTMP,UVMP,FRTMP,IDL,CPER,NPVIB,FVIB,FPTEM,ITMAX included
!REAL(8)  :: NCANGLE(4,180),NCRANGLE(4,180)
REAL(KIND=8) :: FTMP0,FVTMP0
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IVMODEL !--isebasti: included
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INONVHS(:,:) !-- Han add, for nonVHS model

!--NCANGLE sample of angle
!--NCRANGLE sample of angle which react
!--NMFVT count of mf Et and Ev
!--NMFVTR count of MF Et and Ev react
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
!--IDL(K,L) the dissociation vibrational level of species L for mode K
!--CPER(K,L) the numerator constant in post reaction energy redistribution rule for reaction K, molecule L
!--CPER(K,L) the denominator constant
!--NPVIB(J,IKA,M,KV,I) sum of pre/post-collisional (J=1/2) vibrational levels I, mode KV, molecule M involved in reaction IKA
!--FPVIB(IKA,IT,M,KV,I) sampled probability of molecule M reacting from vibrational levels I, mode KV, in reaction IKA at Temperature bin IT
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
!--IREV(N) the index of the corresponding reverse reaction
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
!
END MODULE CALC
!
!********************************************************************
!
MODULE EXPCOL
!
! number of polynomials
INTEGER,PARAMETER :: NT = 10
REAL(8) :: gauss_w(2*NT), gauss_x(2*NT), gx(NT), gw(NT)
! exponential collision model
TYPE :: EXPparType
  INTEGER :: sp1, sp2
  REAL(8) :: alpha, A
  REAL(8) :: totxsecpar(3)
  ! sigma_ref (A^2), Et_ref(eV), alpha
  REAL(8) :: collrate(3)
  ! sigma_Tref (A^2), Tref, beta
END TYPE EXPparType
TYPE(EXPparType),target :: EXPpar(13)
INTEGER, ALLOCATABLE :: EXPCOLPair(:,:)

INTEGER :: temp_int(5)
real(8) :: temp_real(5)

INTERFACE EXPCOL_POT
  module procedure EXPCOL_POT_IPAIR
  module procedure EXPCOL_POT_LSMS
END INTERFACE EXPCOL_POT

contains
  subroutine EXPCOL_INIT()
    USE GAS, only: GASCODE, MSP, INONVHS
    implicit none
    integer :: i

    call cgqf(2*NT, 2, 1.0d0, 1.0d0, -1.0d0, 1.0d0, gauss_x, gauss_w)
    do i = NT+1,2*NT
      ! only take all points with x > 0
      gx(i-NT) = gauss_x(i)
      gw(i-NT) = gauss_w(i)
    end do

    IF (GASCODE .ne. 8) THEN
      STOP "Only gascode = 8 works with EXPCOL model"
    END IF
    ALLOCATE(EXPCOLPair(MSP,MSP))
    EXPCOLPair = 0
    ! The following data is generated by python script

    ! Xsec =totxsecpar(1)*exp(1-(Et/totxsecpar(2))^totxsecpar(3)) ! A^2
    ! Rate = sqrt(8*k*T/pi/mr)*collrate(1)*1e-20*(T/collrate(2))^collrate(3) !si
    !N2-N2  SSE: 1.674089e-03
    EXPCOLPair(1,1) = 1
    INONVHS(1,1) = 2
    EXPpar(1)%sp1 = 1; EXPpar(1)%sp2 = 1
    EXPpar(1)%alpha = 3.160000d+00
    EXPpar(1)%A = 2.290000d+03
    EXPpar(1)%totxsecpar = (/2.507358016232497d+02, 5.063345225872244d-06, 6.106254902258598d-02/)
    EXPpar(1)%collrate = (/1.271457411427805d+02, 2.730000d+02, -1.276419988028019d-01/)

    !N2-N  SSE: 1.809947e-03
    EXPCOLPair(1,2) = 2; EXPCOLPair(2,1) = 2
    INONVHS(1,2) = 2; INONVHS(2,1) = 2
    EXPpar(2)%sp1 = 1; EXPpar(2)%sp2 = 2
    EXPpar(2)%alpha = 3.310000d+00
    EXPpar(2)%A = 6.200000d+02
    EXPpar(2)%totxsecpar = (/2.041874969443462d+02, 6.584258603075655d-06, 6.483460208133425d-02/)
    EXPpar(2)%collrate = (/1.014173134533129d+02, 2.730000d+02, -1.396703798073307d-01/)

    !N2-O2  SSE: 1.884107e-03
    EXPCOLPair(1,3) = 3; EXPCOLPair(3,1) = 3
    INONVHS(1,3) = 2; INONVHS(3,1) = 2
    EXPpar(3)%sp1 = 1; EXPpar(3)%sp2 = 3
    EXPpar(3)%alpha = 3.020000d+00
    EXPpar(3)%A = 1.430000d+03
    EXPpar(3)%totxsecpar = (/2.569702083027821d+02, 8.025301230901390d-06, 6.324566206685799d-02/)
    EXPpar(3)%collrate = (/1.328010841984329d+02, 2.730000d+02, -1.317266718307649d-01/)

    !N2-O  SSE: 8.126512e-04
    EXPCOLPair(1,4) = 4; EXPCOLPair(4,1) = 4
    INONVHS(1,4) = 2; INONVHS(4,1) = 2
    EXPpar(4)%sp1 = 1; EXPpar(4)%sp2 = 4
    EXPpar(4)%alpha = 5.120000d+00
    EXPpar(4)%A = 1.135000d+04
    EXPpar(4)%totxsecpar = (/1.189844871183966d+02, 8.979606905376216d-07, 5.436973638063227d-02/)
    EXPpar(4)%collrate = (/5.642696893546214d+01, 2.730000d+02, -1.154919116834621d-01/)

    !N2-NO  SSE: 1.305621e-03
    EXPCOLPair(1,5) = 5; EXPCOLPair(5,1) = 5
    INONVHS(1,5) = 2; INONVHS(5,1) = 2
    EXPpar(5)%sp1 = 1; EXPpar(5)%sp2 = 5
    EXPpar(5)%alpha = 3.610000d+00
    EXPpar(5)%A = 5.780000d+03
    EXPpar(5)%totxsecpar = (/2.134520265014123d+02, 2.713413428398328d-06, 5.775123483594338d-02/)
    EXPpar(5)%collrate = (/1.065667419343105d+02, 2.730000d+02, -1.203162289118382d-01/)

    !N-O2  SSE: 1.141850e-03
    EXPCOLPair(2,3) = 6; EXPCOLPair(3,2) = 6
    INONVHS(2,3) = 2; INONVHS(3,2) = 2
    EXPpar(6)%sp1 = 2; EXPpar(6)%sp2 = 3
    EXPpar(6)%alpha = 4.130000d+00
    EXPpar(6)%A = 3.870000d+03
    EXPpar(6)%totxsecpar = (/1.630045411489098d+02, 1.830399541957495d-06, 5.776633336730437d-02/)
    EXPpar(6)%collrate = (/7.835324391532143d+01, 2.730000d+02, -1.233829009260243d-01/)

    !N-O  SSE: 1.849318e-03
    EXPCOLPair(2,4) = 7; EXPCOLPair(4,2) = 7
    INONVHS(2,4) = 2; INONVHS(4,2) = 2
    EXPpar(7)%sp1 = 2; EXPpar(7)%sp2 = 4
    EXPpar(7)%alpha = 3.410000d+00
    EXPpar(7)%A = 3.480000d+02
    EXPpar(7)%totxsecpar = (/1.837102196761945d+02, 6.855802188600021d-06, 6.645367929011418d-02/)
    EXPpar(7)%collrate = (/8.984493529798932d+01, 2.730000d+02, -1.457544969235835d-01/)

    !N-NO  SSE: 1.073775e-03
    EXPCOLPair(2,5) = 8; EXPCOLPair(5,2) = 8
    INONVHS(2,5) = 2; INONVHS(5,2) = 2
    EXPpar(8)%sp1 = 2; EXPpar(8)%sp2 = 5
    EXPpar(8)%alpha = 4.210000d+00
    EXPpar(8)%A = 5.333000d+03
    EXPpar(8)%totxsecpar = (/1.630636708160226d+02, 1.396932444288509d-06, 5.659715076228259d-02/)
    EXPpar(8)%collrate = (/7.775878439755058d+01, 2.730000d+02, -1.209190897726498d-01/)

    !O2-O2  SSE: 2.182076e-03
    EXPCOLPair(3,3) = 9
    INONVHS(3,3) = 2
    EXPpar(9)%sp1 = 3; EXPpar(9)%sp2 = 3
    EXPpar(9)%alpha = 2.850000d+00
    EXPpar(9)%A = 8.200000d+02
    EXPpar(9)%totxsecpar = (/2.663120843169848d+02, 1.368006136046089d-05, 6.601000293969347d-02/)
    EXPpar(9)%collrate = (/1.408549017303121d+02, 2.730000d+02, -1.369066052595328d-01/)

    !O2-O  SSE: 8.730286e-04
    EXPCOLPair(3,4) = 10; EXPCOLPair(4,3) = 10
    INONVHS(3,4) = 2; INONVHS(4,3) = 2
    EXPpar(10)%sp1 = 3; EXPpar(10)%sp2 = 4
    EXPpar(10)%alpha = 4.850000d+00
    EXPpar(10)%A = 1.025000d+04
    EXPpar(10)%totxsecpar = (/1.304826660229918d+02, 1.048298471567289d-06, 5.483228200498218d-02/)
    EXPpar(10)%collrate = (/6.229565363151878d+01, 2.730000d+02, -1.161953941888959d-01/)

    !O2-NO  SSE: 1.214724e-03
    EXPCOLPair(3,5) = 11; EXPCOLPair(5,3) = 11
    INONVHS(3,5) = 2; INONVHS(5,3) = 2
    EXPpar(11)%sp1 = 3; EXPpar(11)%sp2 = 5
    EXPpar(11)%alpha = 3.780000d+00
    EXPpar(11)%A = 7.620000d+03
    EXPpar(11)%totxsecpar = (/1.994013057735732d+02, 2.485868550623418d-06, 5.702576039816101d-02/)
    EXPpar(11)%collrate = (/9.976306434042252d+01, 2.730000d+02, -1.182918288496699d-01/)

    !O-NO  SSE: 1.236652e-03
    EXPCOLPair(4,5) = 12; EXPCOLPair(5,4) = 12
    INONVHS(4,5) = 2; INONVHS(5,4) = 2
    EXPpar(12)%sp1 = 4; EXPpar(12)%sp2 = 5
    EXPpar(12)%alpha = 3.950000d+00
    EXPpar(12)%A = 3.140000d+03
    EXPpar(12)%totxsecpar = (/1.725577515613599d+02, 2.405159933956255d-06, 5.875724771821587d-02/)
    EXPpar(12)%collrate = (/8.393646414987728d+01, 2.730000d+02, -1.250446186579113d-01/)

    !NO-NO  SSE: 1.650648e-03
    EXPCOLPair(5,5) = 13
    INONVHS(5,5) = 2
    EXPpar(13)%sp1 = 5; EXPpar(13)%sp2 = 5
    EXPpar(13)%alpha = 3.260000d+00
    EXPpar(13)%A = 2.160000d+03
    EXPpar(13)%totxsecpar = (/2.319680090497029d+02, 5.958044583750658d-06, 6.156689631938237d-02/)
    EXPpar(13)%collrate = (/1.187751445857812d+02, 2.730000d+02, -1.281349445215118d-01/)


  end subroutine EXPCOL_INIT

  subroutine EXPCOL_END()
    implicit none
    deallocate(EXPCOLPair)
  end subroutine EXPCOL_END

  subroutine EXPCOL_TOTXSEC(LS, MS, ET, XSEC, IERROR)
    implicit none
    integer, intent(in) :: LS, MS
    real(8), intent(in) :: ET  !should be in eV
    real(8) :: COLXSEC, XSEC
    real(8),pointer :: p(:)
    integer :: IERROR, iexpair
    iexpair = EXPCOLPair(LS, MS)
    if (iexpair .eq. 0) then
      ! no exponential type data was found
      IERROR = 1
      XSEC = 0.0d0
    else
      IERROR = 0
      p => EXPpar(iexpair)%totxsecpar
      XSEC = p(1)*dexp(1.0d0 - (ET/p(2))**p(3))
    end if
  end subroutine EXPCOL_TOTXSEC

  subroutine EXPCOL_RATE(LS, MS, T, rate)
    use calc,only : PI, BOLTZ, AVOG
    use gas,only : SPM
    implicit none
    integer, intent(in) :: LS, MS
    real(8), intent(in) :: T
    real(8), intent(out) :: rate
    ! calculate collision rates in SI unit: m^3/s
    real(8), pointer :: p(:)
    integer :: iexpair
    iexpair = EXPCOLPair(LS, MS)
    p => EXPpar(iexpair)%collrate
    rate = dsqrt(8*BOLTZ*T/SPM(1,LS,MS)) * p(1)*1d-20*(T/p(2))**p(3)
  end subroutine EXPCOL_RATE

  function EXPCOL_POT_IPAIR(IPAIR, R)
    implicit none
    ! calculate the potential energy
    ! R should be in Angstrom, return in eV
    integer, intent(in) :: IPAIR
    real(8) :: R, EXPCOL_POT_IPAIR
    EXPCOL_POT_IPAIR = EXPpar(IPAIR)%A * dexp(-EXPpar(IPAIR)%alpha * R)
  end function EXPCOL_POT_IPAIR

  function EXPCOL_POT_LSMS(LS, MS, R)
    implicit none
    integer :: LS, MS, ipair
    real(8) :: R, EXPCOL_POT_LSMS
    ipair = EXPCOLPair(LS, MS)
    if (ipair == 0) then
      write(*,*) "No exponential parameters for ",LS,MS
    else
      EXPCOL_POT_LSMS = EXPpar(IPAIR)%A * dexp(-EXPpar(IPAIR)%alpha * R)
    end if
  end function EXPCOL_POT_LSMS

  function EXPCOL_DENOM(ipair, r, b, Et)
    implicit none
    integer :: ipair
    real(8) :: b, Et, r, EXPCOL_DENOM
    EXPCOL_DENOM = 1.0d0 - b**2/r**2 - EXPCOL_POT(ipair, r)/Et
  end function EXPCOL_DENOM

  function EXPCOL_DENOMS(r)
    implicit none
    real(8) :: r, EXPCOL_DENOMS
    EXPCOL_DENOMS = 1.0d0 - temp_real(2)**2/r**2 - EXPCOL_POT(temp_int(1), r)/temp_real(1)
  end function EXPCOL_DENOMS

  function EXPCOL_SOLVERM(LS, MS, b, Et) result(RM)
    implicit none
    integer :: LS, MS
    real(8) :: b, Et, RM, lx, rx, tol=1.d-8, fa0, fb0
    real(8),external :: brent
    temp_int(1) = EXPCOLPair(LS, MS)
    temp_real(1) = Et
    temp_real(2) = b

    lx = 1.0d0
    fa0 = EXPCOL_DENOMS(lx)
    do while (fa0 > 0.0d0 .and. lx .ge. 1.0d-3)
      lx = lx / 1.2d0
      fa0 = EXPCOL_DENOMS(lx)
    end do
    if (lx < 1.0d-3) then
      write(*, '("Fail to solve for min ",I2,1X,I2," Et= ",G11.4," b= ",F7.4)') LS, MS, Et, b
      stop
    end if

    rx = 5.0d0
    fb0 = EXPCOL_DENOMS(rx)
    do while (fb0 < 0.0d0 .and. rx .le. 1.0d2)
      rx = rx * 1.2d0
      fb0 = EXPCOL_DENOMS(rx)
    end do
    if (rx > 1.0d2) then
      write(*, '("Fail to solve for max ",I2,1X,I2," Et= ",G11.4," b= ",F7.4)') LS, MS, Et, b
      stop
    end if

    RM = brent(lx, rx, EXPCOL_DENOMS, tol, fa0, fb0)
  end function EXPCOL_SOLVERM

  function EXPCOL_Chi(LS, MS, b, Et)
    ! calculate chi based on b and Et
    use calc, only : PI
    implicit none
    real(8) :: RM, R(NT), f(NT), EXPCOL_Chi, b, Et
    integer :: i, ipair,  LS, MS
    real(8) :: a

    ipair = EXPCOLPair(LS, MS)
    RM = EXPCOL_SOLVERM(LS, MS, b, Et)
    a = 0.0d0
    do i = 1,NT
      R(i) = RM/gx(i)
      f(i) = dsqrt(1.0d0-gx(i)**2)/dsqrt(EXPCOL_DENOM(ipair,R(i),b,Et))
      a = a + f(i)*gw(i)
    end do

    EXPCOL_Chi = PI-2.0d0*b/RM*a
  end function EXPCOL_Chi

  subroutine EXPCOL_Scatter(LS,MS,BMAX,ET0,VR,VRC,VRCP,IDT)
    ! VR: adjusted speed
    ! VRC: adjusted velocity
    ! VRCP: scattered velocity
    USE CALC, only : EVOLT, PI
    USE GAS, only : SPM
    implicit none
    integer :: LS, MS, IDT, IERROR
    real(8) :: VRC(3), VR, VRCP(3), BMAX,RANF,ET0
    real(8) :: b, CHI, CCHI, SCHI, EPSI, CEPSI, SEPSI, D

    call ZGF(RANF, IDT)
    b = dsqrt(RANF)*BMAX

    CHI = EXPCOL_Chi(LS,MS,b,ET0)
    CCHI = DCOS(CHI); SCHI = DSIN(CHI)
    CALL ZGF(RANF,IDT)
    EPSI = RANF*2.0D0*PI
    CEPSI = DCOS(EPSI)
    SEPSI = DSIN(EPSI)

    D=DSQRT(VRC(2)**2+VRC(3)**2)
    VRCP(1) = CCHI*VRC(1) + SCHI*SEPSI*D
    VRCP(2) = CCHI*VRC(2) + SCHI*(VR*VRC(3)*CEPSI-VRC(1)*VRC(2)*SEPSI)/D
    VRCP(3) = CCHI*VRC(3) - SCHI*(VR*VRC(2)*CEPSI+VRC(1)*VRC(3)*SEPSI)/D
  end subroutine EXPCOL_Scatter

  function EXPCOL_VT(LS, MS, EV, EC)
    use CALC, ONLY: EVOLT
    implicit none
    real(8) :: EV, EC, EXPCOL_VT, Etref, alpha
    integer,intent(in) :: LS, MS
    integer :: ipair
    ipair = EXPCOLPair(LS, MS)
    Etref = EXPpar(ipair)%totxsecpar(2)*EVOLT
    alpha = EXPpar(ipair)%totxsecpar(3)

    EXPCOL_VT = (1.0d0 - EV/EC)
    EXPCOL_VT = EXPCOL_VT * dexp(-((EC - EV)/Etref)**alpha + (EC/Etref)**alpha)
  end function EXPCOL_VT

  function EXPCOL_PMAXEQ(x)
    ! x = Et/Ec
    implicit none
    real(8) :: EXPCOL_PMAXEQ,x
    ! temp_real(1) : zeta
    ! temp_real(2) : alpha
    ! temp_real(3) : Ec
    EXPCOL_PMAXEQ = 1.d0 + (1.0d0 - temp_real(1)/2)*x/(1.0d0-x) - &
      & temp_real(2)*(temp_real(3)*x/temp_real(4))**temp_real(2)
  end function EXPCOL_PMAXEQ

  function EXPCOL_RT(LS,MS,EC,RDOF,IDT)
    ! get a probability E_rot/Ec
    ! EC should be in SI
    ! RDOF dof of rotational energy
    use CALC, only : EVOLT
    implicit none
    integer,intent(in) :: LS,MS,RDOF
    real(8),intent(in) :: EC
    integer :: IDT, samp_t, ipair, i
    real(8) :: EXPCOL_RT,lx,rx,RANF,xs,PROB,xout,fa0,fb0,a,b
    real(8),parameter :: tol = 1.0d-4
    real(8),external :: brent
    ipair = EXPCOLPair(LS,MS)
    temp_real(1) = dble(RDOF)
    temp_real(2) = EXPpar(ipair)%totxsecpar(3)
    temp_real(3) = EC
    temp_real(4) = EXPpar(ipair)%totxsecpar(2)*EVOLT

    if (RDOF == 2) then
      samp_t = 1
    else
      fb0 = EXPCOL_PMAXEQ(0.99d0)
      if (fb0 > 0) then
        samp_t = 1
      else
        samp_t = 2
        fa0 = EXPCOL_PMAXEQ(0.01d0)
      end if
    end if

    if (samp_t == 1) then
      xs = 1.0d0
      a = xs*dexp(-(EC*xs/temp_real(4))**temp_real(2))
      b = 0.d0
    else
      xs = brent(0.01d0,0.99d0,EXPCOL_RT,tol,fa0,fb0)
      b = dble(RDOF)*0.5d0-1.0d0
      a = (1.0d0-xs)**b*xs*dexp(-(EC*xs/temp_real(4))**temp_real(2))
    end if

    i = 1
    do while (i > 0 .and. i < 101)
      call ZGF(xout,IDT)
      PROB = (1.0d0 - xout)**b*xout*dexp(-(EC*xout/temp_real(4))**temp_real(2))/a
      call ZGF(RANF,IDT)
      if (PROB > RANF) then
        i = -1
      else
        i = i+1
      end if
    end do

    if (i == -1) then
      EXPCOL_RT = 1.0d0 - xout
    else
      write(*, "(A)") "EXPCOL_RT fails for EC = ",EC, " RDOF=", RDOF
      stop
    end if
  end function EXPCOL_RT
END MODULE EXPCOL
!********************************************************************
!
MODULE MFDSMC
!
!--declares of the variables used in MF-DSMC model
INTEGER :: IMF,IMFS, NMFpair=0,IMFdia
INTEGER, ALLOCATABLE :: NMFANG(:),IMFpair(:,:)
REAL(8),ALLOCATABLE,DIMENSION(:,:)   :: MFRMASS(:,:),NMFET0,NMFET,NMFETR
REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: NMFER0, NMFEV0, NMFER, NMFEV
REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: NMFERR, NMFEVR
REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:) :: NMFVT0, NMFVT,NMFVTR

! variable to store AHO model parameters
type :: aho_dt
  integer :: iaho
  ! iaho = 0 : SHO phase
  !      = 1 : AHO phase calculated from QCT
  !      = 2 : AHO phase calculated using Morse parameters
  integer :: vmax, mnphase
  real(8),pointer :: vphase(:,:,:)
  real(8) :: morse(5)
end type aho_dt
type(aho_dt),pointer :: aho_dat(:)

! variable to store Morse parameters


contains
  subroutine MF_SET_AHO()
    USE GAS, only : IVMODEL, MSP, GASCODE, ISPV
    USE CALC, only: EVOLT
    implicit none

    integer :: i
    CHARACTER(len=20),ALLOCATABLE :: FILENAME(:)

    IF (GASCODE .ne. 8) THEN
      STOP "Only gascode = 8 works with MF model"
    END IF
    IF (IMF < 2) THEN
      STOP "IMF < 2 but program gets into MF_SET_AHO"
    END IF

    ALLOCATE(aho_dat(MSP), FILENAME(MSP))
    do i=1,MSP
      FILENAME(i) = 'na'
      aho_dat(i)%iaho = 0
    end do

    IF (IMF == 2) THEN
      ! hard coded
      FILENAME(1) = 'n2.bin'
      FILENAME(3) = 'oo.bin'
      FILENAME(5) = 'no.bin'
      OPEN(3,FILE="DS1VD.txt", ACCESS="APPEND")
      write(3,'(a)') "Start setting AHO model for MFDSMC"
      do i=1,MSP
        if (ISPV(i) == 1 .and. FILENAME(i) .ne. 'na' .and. IVMODEL(i,1) ==1 ) then
          aho_dat(i)%iaho = 1
          open(10,file=FILENAME(i),form='unformatted', action='read')
          read(10) aho_dat(i)%vmax
          if (aho_dat(i)%vmax < IVMODEL(i,2)) then
            write(*,'("Specie: ", I2, " don''t have enought vlevel ")') i
            stop
          end if
          read(10) aho_dat(i)%mnphase
          allocate(aho_dat(i)%vphase(aho_dat(i)%mnphase,2,aho_dat(i)%vmax+1))
          read(10) aho_dat(i)%vphase(:,1,:)
          read(10) aho_dat(i)%vphase(:,2,:)
          ! write(*,*) aho_dat(i)%vphase(aho_dat(i)%mnphase,2,aho_dat(i)%vmax+1)
          close(10)
          write(3,'(" SP:", I3, " Vmax:", I3, " MNPHASE: ", I8)') i,aho_dat(i)%vmax, aho_dat(i)%mnphase
        end if
      end do
      close(3)
    ELSE IF (IMF == 3) THEN
      ! hard coded setting's for Morse parameter
      ! The first three parameters are De, Rm and beta
      !     V(r) = De*(1-exp(-beta(x-Rm)))^2 - De
      ! The last one is the reduced mass of the molecule

      ! N2
      aho_dat(1)%iaho = 2
      aho_dat(1)%morse(1) = 9.904795d0*EVOLT  !Joule, Morse zero point energy
      aho_dat(1)%morse(2) = 1.0977d0          !A
      aho_dat(1)%morse(3) = 2.81971d10        !m^{-1}
      aho_dat(1)%morse(4) = 2.32587d-26       !kg
      aho_dat(1)%morse(5) = 9.82163d0*EVOLT   !Joule, QCT zero point energy
      ! morse(5) check around line 6182

      ! NO
      aho_dat(5)%iaho = 2
      aho_dat(5)%morse(1) = 6.623d0*EVOLT     !Joule
      aho_dat(5)%morse(2) = 1.159d0           !A
      aho_dat(5)%morse(3) = 2.83d10           !m^{-1}
      aho_dat(5)%morse(4) = 2.480328d-26      !kg
      aho_dat(5)%morse(5) = 6.55879d0*EVOLT   !Joule, check around line 6293

      ! O2
      aho_dat(3)%iaho = 2
      aho_dat(3)%morse(1) = 5.211d0*EVOLT     !Joule
      aho_dat(3)%morse(2) = 1.207d0           !A
      aho_dat(3)%morse(3) = 2.78d10           !m^{-1}
      aho_dat(3)%morse(4) = 2.65676d-26       !kg
      aho_dat(3)%morse(5) = 5.21275d0*EVOLT   !Joule, check around line 6236
    ENDIF
    DEALLOCATE(FILENAME)
  end subroutine MF_SET_AHO

  subroutine MF_CLEAN_AHO()
    implicit none
    integer :: i
    if (associated(aho_dat)) then
      do i=1,size(aho_dat)
        if (associated(aho_dat(i)%vphase)) then
          deallocate(aho_dat(i)%vphase)
        end if
      end do
      deallocate(aho_dat)
    end if
  end subroutine MF_CLEAN_AHO

  subroutine MF_SAMPLE_PHASE(LS, R, Vlevel, Evib)
    use CALC, only : PI
    implicit none
    integer,intent(in) :: LS, Vlevel
    real(8),intent(in) :: Evib
    real(8) :: R
    ! LS: specie of the molecule
    ! R: a random number between 0 and 1

    ! variables for QCT based AHO
    integer :: ii, ij, itry, imaxtry = 5
    logical :: ierror

    ! variables for Morse based AHO
    real(8) :: rho, theta0, R2

    if (IMF == 1) then
      R = R*PI
      ! it shouldn't make difference if just use R*2.0*PI
    else
      if (aho_dat(LS)%iaho == 0) then
        R = R*PI
      else if (aho_dat(LS)%iaho == 1) then
        ierror = .false.
        do itry = 1, imaxtry
          call binary_search(ii,ij,R,aho_dat(LS)%vphase(:,1,Vlevel+1),&
            1, aho_dat(LS)%mnphase, ierror)
          if (.not. ierror) then
            exit
          end if
        end do
        if ( .not. ierror) then
          if (ii == ij) then
            R = aho_dat(LS)%vphase(ii,2,Vlevel+1)
          else
            R = R - aho_dat(LS)%vphase(ii,1,Vlevel+1)
            R = R/(aho_dat(LS)%vphase(ii,1,Vlevel+1) - aho_dat(LS)%vphase(ij,1,Vlevel+1))
            R = R*(aho_dat(LS)%vphase(ii,2,Vlevel+1) - aho_dat(LS)%vphase(ij,2,Vlevel+1))
            R = R + aho_dat(LS)%vphase(ii,2,Vlevel+1)
          end if
        else
          write(*,'(A, G15.6)') "Failed to sampling for ", R
        end if
      else if (aho_dat(LS)%iaho == 2) then
        ! omega0 = aho_dat(LS)%morse(3)*dsqrt(aho_dat(LS)%morse(1)/aho_dat(LS)%morse(4))/1d15
        ! rad/fs
        rho = 1.0d0 + (Evib - aho_dat(LS)%morse(5))/aho_dat(LS)%morse(1)
        if (rho < 0.0d0) then
          rho = 0.0d0
        else
          rho = dsqrt(rho)
        endif
        theta0 = dasin(rho)
        R2 = dcos(theta0)*dcos(2*PI*R - theta0)/(1 + rho*dsin(2*PI*R-theta0))
        IF (R2 < -1.0d0) R2 = -1.0d0
        IF (R2 > 1.0d0)  R2 = 1.0d0
        if (R < theta0/2.0d0/PI+0.5d0) then
          R = dacos(R2)
        else
          R = 2.0d0*PI - dacos(R2)
        end if
      end if
    end if
  end subroutine MF_SAMPLE_PHASE

  subroutine MF_CALC_COLL(K,NMFCALL,MFcoll,IDT)
    use gas,only : IREA,SP
    implicit none
    real(8) :: MFcoll(2),AA
    integer,intent(in) :: K, NMFCALL
    integer :: II,IDT

    IF (NMFCALL == 1) THEN
      DO II=1,2
        IF (IREA(II,K) == 1 .or. IREA(II,K) == 3) THEN
                  ! homo-nuclear molecule like N2 or O2
                  ! MFcoll is the mass of atom
          MFcoll(II) = SP(5,IREA(II,K))/2.0
        ELSE IF (IREA(II,K) == 2 .or. IREA(II,K) == 4) THEN
          MFcoll(II) = SP(5,IREA(II,K))
        ELSE IF (IREA(II,K) == 5) THEN  !randomly select one
                  ! NO dissocation reaction
          CALL ZGF(AA,IDT)
          IF (AA <0.5d0) THEN
            MFcoll(II) = 2.325d-26 ! N atom
          ELSE
            MFcoll(II) = 2.656d-26
          END IF
        ELSE
          STOP "CHECK_RXSECTION: MFcoll is not found"
        END IF
      END DO
    ELSE
      ! at this point MFcoll should already be calculated
      AA = MFcoll(1)
      MFcoll(1) = MFcoll(2)
      MFcoll(2) = AA
    END IF
  end subroutine MF_CALC_COLL

  subroutine MF_SAMPLE_ANGLE(K,NMFCALL,MFANG,MFCtheta,viblevel,MFV1,MFV2,IDT)
    use calc, only:PI
    use gas, only : IREA
    implicit none
    integer :: II, IDT
    integer,intent(in) :: K,viblevel(2), NMFCALL
    real(8),intent(in) :: MFV1, MFV2
    real(8) :: MFANG(8),MFCtheta,AA
    ! K --- index of reaction
    ! nmfcall -- the time that mf model is called
    ! mfang -- angles sampled for MF model
    ! mfctheta -- cos(theta)
    ! viblevel -- the vibrational quantum number of the colliding particles

    ! MFANG:
    !  1. gamma_1
    !  2. gamma_2
    !  3. theta
    !  4. phi_0
    !  5. phi_1
    !  6. beta_1
    !  7. beta_2
    !  8. delta

    IF (NMFCALL == 1) THEN
      !-- first time call MF model, generate angles
      DO II = 1,3
        CALL ZGF(MFANG(II),IDT)
        MFANG(II) = MFANG(II)*PI
      END DO
      CALL ZGF(MFANG(4),IDT)
      CALL MF_SAMPLE_PHASE(IREA(1,K),MFANG(4),viblevel(1),MFV1)
      IF (NMFANG(K) .ge. 8) THEN !collider is diatom
        CALL ZGF(MFANG(5),IDT)
        CALL MF_SAMPLE_PHASE(IREA(2,K),MFANG(5),viblevel(2),MFV2)
        DO II=6,8
          CALL ZGF(MFANG(II),IDT)
          MFANG(II) = MFANG(II)*PI
        END DO
        MFANG(3) = MFANG(3)*2.0d0    ! For atom-diatom, no symmetry
        MFANG(6) = MFANG(6)-PI*0.5d0 ! beta1 [-pi/2, pi/2]
        MFANG(7) = MFANG(7)*2.0d0    ! beta2 [0,2*pi]
        MFANG(8) = MFANG(8)*2.0d0    ! delta [0,2*pi]
      END IF
      MFCtheta = DCOS(MFANG(3))
    ELSE IF(NMFCALL == 2) THEN
      !--this is the second call
      AA=MFANG(4); MFANG(4)=MFANG(5); MFANG(5)=AA !switch phi
      ! we don't change other angles
      MFCtheta = DCOS(MFANG(6))*DCOS(MFANG(7))
    ELSE
      STOP
    END IF
  end subroutine MF_SAMPLE_ANGLE

  subroutine MF_EVAL_F(K,NMFCALL,MFV1,MFV2,MFR1,MFR2,MFcoll, MFANG, MFCtheta,MFF,IDT)
    ! calculate the threshold function
    use calc, only:PI
    USE GAS, ONLY : IREA, REA,ISPV
    implicit none
    integer,intent(in) :: K, NMFCALL
    real(8),intent(in) :: MFV1,MFV2,MFR1,MFR2,MFcoll(2),MFANG(8),MFCtheta
    real(8),intent(out) :: MFF
    integer :: IDT
    real(8) :: MFDSTAR, AA, BB, CC, DD, MFbeta
    ! K:
    !     number of the reaction
    ! MFV1, MFR1
    !     vibrational/rotational energy of the molecule to be dissociated AB
    ! MFV2, MFR2
    !     vibrational/rotational energy, vibrational level of the colliding partner CD
    ! MFcoll(1), MFcoll(2)
    !     mass of atom B, mass of atom C
    ! MFANG, MFCtheta
    !     MF angles sampled for AB dissociation
    ! NMFCALL
    !     =1 : the MFANG is sampled for AB + CD -> A + B + CD
    !     =2 : the MFANG is sampled for AB + CD -> AB + C + D
    MFDSTAR = REA(2,K) - MFR1 + 2.0d0*MFR1**1.5d0/(3.0d0 * dsqrt(3.0d0*2.0d0*REA(2,K)))
    AA = MFDSTAR - MFV1*DSIN(MFANG(4))**2
    IF (AA .LE. 0.0D0 )THEN
      MFF = 0.0D0
    ELSE
      AA = (DSQRT(AA)+DSQRT(MFV1)*DCOS(MFANG(4)))/MFCtheta
      AA = AA*(MFcoll(1)/MFcoll(2)+1.0d0)

      BB =-2.0d0*MFRMASS(2,K)/MFcoll(1)*DSQRT(MFV1)*DCOS(MFANG(4))
      ! MFRMASS(2,K) the mass of the molecule to be dissociated
      BB = BB*MFCtheta

      CC = 0.0d0
      IF (ISPV(IREA(2,K)) == 1) THEN
        IF (NMFCALL == 1) THEN
          CC = DCOS(MFANG(8))*DCOS(MFANG(7))*DSIN(MFANG(6))+DSIN(MFANG(8))*DSIN(MFANG(7))
          ! cos(delta) * cos(beta_2) * cos(beta_1) + sin(delta)*sin(beta_2)
          DD = DCOS(MFANG(5))*DCOS(MFANG(6))*DCOS(MFANG(7))
          ! cos(phi_1) * cos(beta_1) * cos(beta_2)
          CC = DSQRT(MFR2)*CC
          DD = DSQRT(MFV2)*DD
          CC = CC-DD
        ELSE
          CALL ZGF(MFbeta,IDT)
          MFbeta = PI*(MFbeta-0.5d0) ! sample beta angle
          CC = DSQRT(MFR2)*DCOS(MFbeta)*DSIN(MFANG(3))
          DD = DSQRT(MFV2)*DCOS(MFANG(5))*DCOS(MFANG(3))
          CC = -CC-DD
        END IF
        CC = -2.0d0*CC*DSQRT(MFRMASS(2,K)*MFRMASS(3,K))/MFcoll(2)
      END IF

      MFF = (AA+BB+CC)**2
      MFF = MFF*MFRMASS(1,K)/MFRMASS(2,K)/(DCOS(MFANG(1))*DCOS(MFANG(2)))**2
      MFF = MFF / 4.d0
    END IF
  end subroutine MF_EVAL_F

!
!--IMF only applied when QCTMODEL = 1, 0 for TCE, 1 for MFDSMC, 2 for MFDSMC-AHO
!--IMFS = 0, not sample anything related to MF model =1 sample
!--nonVHS flag to control nonVHS model for N2+O
!--NMFANG = number of angle to sample for MF model, <8 for atom-diatom, >= 8 for diatom-diatom
!--NMFpair number of molecular pairs to use MF dissociation model
!--IMFpair store the index of different molecular pair
!--MFRMASS store reduced mass

END MODULE MFDSMC
!
!
!********************************************************************
!
MODULE OUTPUT
!
!--declares the variables associated with the sampling and output
!
IMPLICIT NONE
!
INTEGER :: NSAMP,NMISAMP,NDISSOC,NRECOMB,NBINS,NSNAP,NSPDF,NSCELLS,IENERS  !--isebasti: included NBINS,NSNAP,NSPDF
INTEGER :: NOUT
INTEGER, DIMENSION(0:100) :: NDISSL
INTEGER(8), ALLOCATABLE :: NDROT(:,:),NDVIB(:,:,:,:),NSVEC(:)
REAL(KIND=8):: TISAMP,XVELS,YVELS,AVDTM,DBINV(3),DBINC(3),DBINE(3),ENERS,ENERS0 !--isebasti: DBINV,C,ENERS,ENERS0 included
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: COLLS,WCOLLS,CLSEP,REAC,SREAC
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: VAR,VARS,CSSS,CST,PDF,BIN !--isebasti:CST,PDF,BIN included
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CS,VARSP,PDFS,BINS !--isebasti:PDFS,BINS included
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: CSS
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: EVREM
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
!--CS(8+K,N,L) sampled sum of vibrational energy of
!              species L in cell N for mode K
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
!--WCOLLS(N) cumulative weighted number
!--CLSEP(N) the total collision partner separation distance in sampling cell N
!
!--The following variables apply in the sampling of distribution functions
!--(some are especially for the dissociation of oxygen
!
!--IPDF 0 no pdf sampling, 1 with pdf sampling
!--NSPDF cell considered in pdf sampling (see SAMPLE_PDF subroutine)
!--NDISSOC the number of dissociations
!--NRECOMB the number of recombinations
!--NDISSL(L) the number of dissociations from level
!--NDROT(L,N) number of species L in rotational energy range N
!--NDVIB(NS,K,L,N) number of species L in vibrational level N and mode K within cell NS
!--ENERS total energy sampling
!--ENERS0 previous sampled value
!--IPRS post-reaction vibrational levels are populated using
!       0 LB, 1 pre-reaction sampled pdf, 2 for Bondar's approach
!       Hydrogen-Oxygen Detonation Study by DSMC method (2011)
!
END MODULE OUTPUT
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
INTEGER :: IRUN,NSEED,N,I,J,IRETREM,ISET, NTHREADS, ITHREAD
REAL(KIND=8) :: A,WCLOCK(5)
INTEGER(KIND=8) :: COUNT0, COUNT1, COUNT_RATE
REAL(8) :: CALC_TIME

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
IREAC=2      !0 reaction on and rate sampling off (standard case)
!             1 reaction on and rate sampling on (samples only NSPDF cell, requires ISF>0)
!             2 reaction off and rate sampling on (called by run_ireac_ds1v.sh, samples only NSPDF cell, requires ISF=0 & NSEED>= 0)
QCTMODEL=1   !0 for LB+SHO vibrational levels, 1 for TCE+LB+AHO, 2 MEQCT+AHO (for O2+O and N2+O)
IMF = 1      ! 0: do not use MF model, 1: use MF+SHO 2: use MF+AHO(QCT ladder) 3: use MF+AHO(Morse Potential)
IMFdia = 1   ! 0: for collision of same molecules, dissociate the one with higher vibrational energy
             ! 1: for collision of same molecules, dissociate the one gives lower threshold energy
IMFS = 1     ! 0 for not sampling 1 for sampling
nonVHS = 1
!-- nonVHS 0 for VHS/VSS model
!          1 for VHS/VSS model+QCT N2O model
!          2 exponential model, fallback to 0 if parameters not found
!
!--variables for vibratioal sampling
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
  ISF=1                    !to restart a case that is already ISF=0 as an unsteady-case
  IF (NSEED == 9999) ISF=0 !special code for shockwave sampling
  CALL OXYGEN_NITROGEN     !special code to update reaction rates
  CALL INIT_REAC           !special code to update reaction rates
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
  IF (IREAC == 0) THEN !extracting MeQct Zv(T) values
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
  IF (IREAC == 2) THEN
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
  !CALL WRITE_RESTART
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
OPEN(119, FILE="RunningTime.dat", STATUS="REPLACE")
WRITE(119,"(10(A,1X,I2,2X))") "# IRELAX:",IRELAX,"IREAC:",IREAC,"nonVHS:",nonVHS,"GASCODE:",GASCODE
WRITE(119,"(10(A,1X,I2,2X))") "# QCTMODEL:",QCTMODEL,"IMF:",IMF,"IMFdia:",IMFdia,"IMFS:",IMFS
WRITE(119,"(A,G14.6)") '# SAMPRAT = ', SAMPRAT
WRITE(119,"(A,G14.6)") '# TPOUT = ', TPOUT
WRITE(119,"(A,I5)") '# Nthreads = ', NTHREADS
WRITE(119,"(10(A,1X))") 'VARIABLES = ','"NOUT",', '"CALC_TIME (s)",', '"FTIME (s)",', &
                  &   '"DTM (s)",', '"NMOL,"', '"NSAMP",', '"TISAMP"'
CALL SYSTEM_CLOCK(COUNT=COUNT0, COUNT_RATE=COUNT_RATE)
CLOSE(119)
DO WHILE (FTIME < TLIM)
  OPEN (9,FILE='DIAG.TXT',ACCESS='APPEND')
!
!$ WCLOCK(2)=omp_get_wtime()
!
  DO N=1,INT(TPOUT)
!
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
      IF (ISF == 0) CALL SAMPLE_PDF
      IF ((ISF == 1).AND.(N/TPOUT >= (1.d0-FRACSAM))) CALL SAMPLE_PDF
      IF ((ISF == 2).AND.(N >= INT(TPOUT-10))) CALL SAMPLE_PDF
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
!
  END DO
!
  !$ WCLOCK(1)=omp_get_wtime()
  CALL OUTPUT_RESULTS
  CALL SYSTEM_CLOCK(COUNT=COUNT1)
  CALC_TIME = dble(COUNT1-COUNT0) / dble(COUNT_RATE)
  OPEN(119, FILE="RunningTime.dat", POSITION="APPEND")
  WRITE(119,118) NOUT, CALC_TIME, FTIME, DTM, NM, NSAMP, TISAMP
  CLOSE(119)
118  FORMAT(I10,2X, 10(G14.6,2X))
  IF (IREAC .ne. 2) THEN
    ! do not write restart for ireac = 2, file is too huge
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
!
END DO
CALL MF_CLEAN_AHO()
!
STOP
END PROGRAM DS1
!
!***************************************************************************
!
SUBROUTINE READ_DATA
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE MFDSMC, only : IMF, IMFS, IMFdia
!
IMPLICIT NONE
!
INTEGER :: NVERD,MVERD,N,K
!
WRITE (*,*) 'Reading the data file DS1VD.in'
OPEN (4,FILE='DS1VD.in')
OPEN (3,FILE='DS1VD.txt')
!
WRITE (3,*) 'Data summary for program DS1'
!
READ (4,*) NVERD
WRITE (3,*) 'The n in version number n.m is',NVERD
READ (4,*) MVERD
WRITE (3,*) 'The m in version number n.m is',MVERD
!
READ (4,*) AMEG
WRITE (3,*) 'The approximate number of megabytes for the calculation is',AMEG
!
READ (4,*) IFX
IF (IFX == 0) WRITE (3,*) 'Plane flow'
IF (IFX == 1) WRITE (3,*) 'Cylindrical flow'
IF (IFX == 2) WRITE (3,*) 'Spherical flow'
JFX=IFX+1
!
READ (4,*) XB(1)
WRITE (3,*) 'The minimum x coordinate is',XB(1)
!
READ (4,*) ITYPE(1)
IF (ITYPE(1) == 0) WRITE (3,*) 'The minimum x coordinate is a stream boundary'
IF (ITYPE(1) == 1) WRITE (3,*) 'The minimum x coordinate is a plane of symmetry'
IF (ITYPE(1) == 2) WRITE (3,*) 'The minimum x coordinate is a solid surface'
IF (ITYPE(1) == 3) WRITE (3,*) 'The minimum x coordinate is a vacuum'
!
READ (4,*) XB(2)
WRITE (3,*) 'The maximum x coordinate is',XB(2)
!
READ (4,*) ITYPE(2)
IF (ITYPE(2) == 0) WRITE (3,*) 'The maximum x coordinate is a stream boundary'
IF (ITYPE(2) == 1) WRITE (3,*) 'The maximum x coordinate is a plane of symmetry'
IF (ITYPE(2) == 2) WRITE (3,*) 'The maximum x coordinate is a solid surface'
IF (ITYPE(2) == 3) WRITE (3,*) 'The maximum x coordinate is a vacuum'
!
IF (IFX > 0) THEN
  READ (4,*) IWF
  IF (IWF == 0) WRITE (3,*) 'There are no radial weighting factors'
  IF (IWF == 1) WRITE (3,*) 'There are radial weighting factors'
  IF (IWF == 1) THEN
    READ (4,*) WFM
    WRITE (3,*) 'The maximum value of the weighting factor is ',WFM+1.D00
    WFM=WFM/XB(2)
  END IF
END IF
!
READ (4,*) JCD              !--isebasti: reading JCD before GASCODE
IF (JCD == 0) WRITE (3,*) ' Any chemical reactions are based on the continuum rate equations'
IF (JCD == 1) WRITE (3,*) ' Any chemical reactions are based on quantum kinetic theory'
!
READ (4,*) GASCODE
WRITE (3,*) GASCODE
IF (GASCODE == 1) THEN
  WRITE (3,*) 'Hard sphere gas'
  CALL HARD_SPHERE
END IF
IF (GASCODE == 2) THEN
  WRITE (3,*) 'Argon'
  CALL ARGON
END IF
IF (GASCODE == 3) THEN
  WRITE (3,*) 'Ideal gas'
  READ (4,*) IRM
  IF (IRM == 1) WRITE (3,*) ' N2 (Trans+Rot; Zr const)'
  IF (IRM == 2) WRITE (3,*) ' N2 (Trans+Rot+Vib; Zv const)'
  IF (IRM == 3) WRITE (3,*) ' N2 (Trans+Rot+Vib; Zv(T))'
  IF (IRM == 4) WRITE (3,*) ' H2O (Trans+Rot+Vib)'
  IF (IRM == 5) WRITE (3,*) ' H2O+H2O (Trans+Rot+Vib)'
  IF (IRM == 6) WRITE (3,*) ' H2O+H (Trans+Rot+Vib)'
  CALL IDEAL_GAS
END IF
IF (GASCODE == 4) THEN
  WRITE (3,*) ' Real oxygen'
  CALL REAL_OXYGEN
END IF
IF (GASCODE == 5) THEN
  WRITE (3,*) 'Ideal air'
  CALL IDEAL_AIR
END IF
IF (GASCODE == 6) THEN
  WRITE (3,*) 'Real air @ 7.5 km'
  CALL REAL_AIR
END IF
IF (GASCODE == 7) THEN
  WRITE (3,*) 'Helium-xenon mixture with mol. wt. of argon'
  CALL HELIUM_XENON
END IF
IF (GASCODE == 8) THEN
  WRITE (3,*) 'Oxygen-hydrogen'
  READ (4,*) IRM
  IF (IRM == 0)    WRITE (3,*) ' Reactions are not included'
  IF (IRM == 1)    WRITE (3,*) ' H2/O2 full mechanism by Shatalov (2009)'
  IF (IRM == 2)    WRITE (3,*) ' H2/O2 full mechanism by Davidenko (2006)'
  IF (IRM >= 100)  WRITE (3,*) ' A single forward reaction'
  IF (IRM >= 200)  WRITE (3,*) ' Both forward and reverse reactions'
  CALL OXYGEN_NITROGEN
END IF
!
WRITE (3,*) 'Approach for post-reaction vibrational energy redistribution:',IPRS
WRITE (3,*) 'Approach for reaction sampling and thermal relaxation:',IREAC,IRELAX
WRITE (3,*) 'Model for vibrational energies and VT relaxation:',QCTMODEL
WRITE (3,*) ' by species:',IVMODEL(:,1)
!
READ (4,*) IGS
IF (IGS==0) WRITE (3,*) 'The flowfield is initially a vacuum'
IF (IGS==1) WRITE (3,*) 'The flowfield is initially the stream(s) or reference gas'
IF (IGS==2) WRITE (3,*) 'The flowfield is initially the stream(s) or reference gas with Temp-Nden distributions' !--isebasti: included
IF (IGS==3) WRITE (3,*) 'The flowfield is initially the stream(s) or reference gas with a specific vibrational population' !--isebasti: included
IF (IGS==4) WRITE (3,*) 'The flowfield is initially the stream(s) or reference gas with a specific bimodal T-R energy dist' !--isebasti: included
READ (4,*) ISECS
IF (ISECS == 0) WRITE (3,*) 'There is no secondary stream initially at XS > 0'
IF (ISECS == 1) WRITE (3,*) 'There is a secondary stream applied initially at XS > 0 (XB(2) must be > 0)'
IF (ISECS == 1) THEN  !--isebasti: this and line below were modified
  IF ((IWF == 1).AND.(IFX > 0)) THEN
    WRITE (3,*) 'There cannot be a secondary stream when weighting factors are present'
    STOP
  END IF
END IF
READ (4,*) XS
WRITE (3,*) 'The secondary stream boundary/auxiliar interface is at x=',XS
!
DO K=1,2
  IF ((K == 1).OR.(ITYPE(2) == 0)) THEN  !--isebasti: read XB(2) stream properties here even when ISECS=0
    IF (K == 1) WRITE (3,*) 'The XB(1) stream and/or reference gas properties are:-'
    IF (K == 2) WRITE (3,*) 'The XB(2) stream and/or secondary stream reference gas properties are:-'
!
    READ (4,*) FND(K)
    WRITE (3,*) '    The stream number density is',FND(K)
!
    READ (4,*) FTMP(K)
    WRITE (3,*) '    The stream temperature is',FTMP(K)
!
    READ (4,*) FRTMP(K)
    WRITE (3,*) '    The stream rotational temperature is',FRTMP(K)
!
    IF (MMVM > 0) THEN
      READ (4,*) FVTMP(K)
      WRITE (3,*) '    The stream vibrational temperature is',FVTMP(K)
    END IF
!
    READ (4,*) VFX(K)
    WRITE (3,*) '    The stream velocity in the x direction is',VFX(K)
    IF ((NVERD > 1).OR.((NVERD == 1).AND.(MVERD > 6))) THEN
      READ (4,*) VFY(K)
      WRITE (3,*) '    The stream velocity in the y direction is',VFY(K)
    END IF
!
    DO N=1,MSP
      READ (4,*) FSP(N,K)
      WRITE (3,*) '    The fraction of species',N,' is',FSP(N,K)
    END DO
  END IF
!
  IF (ITYPE(K) == 2) THEN
    IF (K == 1) WRITE (3,*) 'The minimum x boundary is a surface with the following properties'
    IF (K == 2) WRITE (3,*) 'The maximum x boundary is a surface with the following properties'
!
    READ (4,*) TSURF(K)
    WRITE (3,*) '     The temperature of the surface is',TSURF(K)
!
    READ (4,*) FSPEC(K)
    WRITE (3,*) '     The fraction of specular reflection at this surface is',FSPEC(K)
!
    READ (4,*) VSURF(K)
    WRITE (3,*) '     The velocity in the y direction of this surface is',VSURF(K)
  END IF
END DO
!
IF (ITYPE(1) == 0) THEN
  READ (4,*) IREM
  IF (IREM == 0) THEN
    WRITE (3,*) ' There is no molecule removal'
    XREM=XB(1)-1.D00
    FREM=0.D00
  ELSE IF (IREM == 1) THEN
    READ (4,*) XREM
    WRITE (3,*) ' Full removal of molecules between',XREM,' and',XB(2)
    FREM=1.D00
  ELSE IF (IREM == 2) THEN
    WRITE (3,*) ' Molecule removal is specified whenever the program is restarted'
    XREM=XB(1)-1.D00
    FREM=0.D00
  END IF
ELSE
  XREM=XB(1)-1.D00
  FREM=0.D00
END IF
!
!--set the speed of the outer boundary
!
IVB=0
VELOB=0.D00
IF (ITYPE(2) == 1) THEN
  READ (4,*) IVB
  IF (IVB == 0) WRITE (3,*) ' The outer boundary is stationary'
  IF (IVB == 1) THEN
    WRITE (3,*) ' The outer boundary moves with a constant speed'
    READ (4,*) VELOB
    WRITE (3,*) ' The speed of the outer boundary is',VELOB
  END IF
END IF
!
WRITE (3,*) 'Computational Parameters'
READ (4,*) MOLSC
WRITE (3,*) ' The target number of molecules per sampling cell is',MOLSC
READ (4,*) NCIS
WRITE (3,*) ' The number of collision cells per sampling cell is ',NCIS
READ (4,*) CPDTM
WRITE (3,*) ' Maximum collisions in a time step is',CPDTM
READ (4,*) TPDTM
WRITE (3,*) ' Maximum sampling cell transits in a time step is',TPDTM
READ (4,*) NNC
IF (NNC == 0) WRITE (3,*) ' Collision partners are selected randomly from the collision cell  '
IF (NNC == 1) WRITE (3,*) ' Nearest neighbor collisions are employed '
READ (4,*) IMTS
IF (IMTS == 0) WRITE (3,*) ' The move time step is uniform over the cells  '
IF (IMTS == 1) WRITE (3,*) ' The move time step can vary over the cells '
READ (4,*) SAMPRAT
WRITE (3,*) ' The number of time steps in a sampling interval is ',SAMPRAT
READ (4,*) OUTRAT
WRITE (3,*) ' The number of sampling intervals in an output interval is ',OUTRAT
READ (4,*) ISF
IF (ISF == 0) WRITE (3,*) ' The sampling is for an eventual steady flow '
IF (ISF == 1) WRITE (3,*) ' The sampling is for a continuing unsteady flow (sampling size is depends on FRACSAM*OUTRAT)'
IF (ISF == 2) WRITE (3,*) ' The sampling is for a continuing unsteady flow (sampling size is at most 10)'
READ (4,*) FRACSAM
WRITE (3,*) ' Any unsteady sample/update is over the final ',FRACSAM,' of the output interval'
!
!READ (4,*) JCD !--isebasti: reading JCD before GASCODE (see above)
!IF (JCD == 0) WRITE (3,*) ' Any chemical reactions are based on the continuum rate equations'
!IF (JCD == 1) WRITE (3,*) ' Any chemical reactions are based on quantum kinetic theory'
!
!--isebasti: lines bellow are commented to force TCE model and avoid allocation bugs
!IF ((GASCODE == 4).OR.(GASCODE == 6).OR.(GASCODE == 8)) THEN
!--deallocate EARLIER ALLOCATION
!DEALLOCATE (NEX,NSPEX,SPEX,ISPEX,TREACG,PSF,TREACL,TNEX,STAT=ERROR)
!IF (ERROR /= 0) THEN
!  WRITE (*,*)'PROGRAM COULD NOT DEALLOCATE INITIAL GAS VARIABLES',ERROR
!END IF
!!--set the reaction data
!  IF (GASCODE == 4) CALL REAL_OXYGEN
!  IF (GASCODE == 6) CALL REAL_AIR
!  IF (GASCODE == 8) CALL OXYGEN_NITROGEN
!END IF
!
READ (4,*) IRECOM
IF (IRECOM == 0) WRITE (3,*) ' Any recombination reactions are not included'
IF (IRECOM == 1) WRITE (3,*) ' Any recombination reactions are included'
!
READ (4,*) IFI !--isebasti: included IFI option
IF (IFI == 0) WRITE (3,*) ' Forced ignition is not allowed',IFI
IF (IFI >  0) WRITE (3,*) ' Forced ignition is allowed',IFI

!=================== ALL the following is commented out

READ(4,*) nonVHS
WRITE(3, *) ' noVHS = ',nonVHS


!=== set model in DS1 head region
READ (4,*) IMF !--han: flag to control Macheret-Fridman model
WRITE(3, *) ' IMF = ',IMF
IF (IMF == 0) THEN
  WRITE(3,*) '   Macheret-Fridman model will not be applied anywhere'
ELSE IF (IMF == 1) THEN
  WRITE(3,*) '   Macheret-Fridman-Monte-Carlo will be used for dissociation'
ELSE IF (IMF == 2) THEN
  WRITE(3,*) '   Macheret-Fridman-Monte-Carlo model with AHO vibrational phase calculated from QCT'
ELSE IF (IMF == 3) THEN
  WRITE(3,*) '   MAcheret-Fridman-Monte-Carlo model with AHO vibrational phase calculated from Morse parameter'
ELSE
  WRITE(3,*) '   ERROR: Wrong input for Macheret-Fridman model'
  stop
ENDIF

IF (IMF .NE. 0) THEN
  READ(4, *) IMFdia
  IF (IMFdia == 0) THEN
    WRITE(3,*) ' IMFdia=0, dissociate the diatomic molecule with higher vibrational energy'
  ELSE IF (IMFdia == 1) THEN
    WRITE(3,*) ' IMFdia=1, dissociate the diatomic molecule with lower threshold energy'
  ELSE
    WRITE(3,*) ' WRONG INPUT FOR IMFdia'
    STOP
  ENDIF
ENDIF

IF ( IMF .EQ. 0 ) IMFS = 0
IF (QCTMODEL == 2 .AND. IMF .NE. 0)THEN
  WRITE(3,*) ' WARNING: MF MODEL BE NOT BE USED!'
ENDIF
!==============================================================
!
CLOSE (3)
CLOSE (4)
!
RETURN
END SUBROUTINE READ_DATA
!
!***************************************************************************
!
SUBROUTINE SET_INITIAL_STATE
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE MFDSMC,only : IMF, IMFS
!
IMPLICIT NONE
!
INTEGER :: I,J,L,K,KK,KN,NMI,II,III,INC,NSET,NSC,NC,NDES,IDES(10000),IV,IDT=0 !--isebasti: included NC,IDT
INTEGER(KIND=8) :: N,M
REAL(KIND=8) :: AA,A,B,BB,SN,XMIN,XMAX,WFMIN,DENG,TDF,DNF,VMP2,RANF,EVIB,TVAR(10000),BVMP,BTEMP !--isebasti: included TDF,DNF,VMP2,RANF,EVIB
REAL(KIND=8), DIMENSION(3) :: DMOM
REAL(KIND=8), DIMENSION(3,2) :: VB
REAL(KIND=8), DIMENSION(2) :: ROTE
REAL(8),EXTERNAL :: ERF,GAM
!
!--NSET the alternative set numbers in the setting of exact initial state
!--DMOM(N) N=1,2,3 for x,y and z momentum sums of initial molecules
!--DENG the energy sum of the initial molecules
!--VB alternative sets of velocity components
!--ROTE alternative sets of rotational energy
!--NMI the initial number of molecules
!--INC counting increment
!
DMOM=0.D00
DENG=0.D00
!--set the number of molecules, divisions etc. based on stream 1
!
NMI=10000*AMEG+0.999999999D00
NDIV=NMI/MOLSC !--MOLSC molecules per division
WRITE (9,*) 'The number of divisions is',NDIV
!
MDIV=NDIV
ILEVEL=0
!
ALLOCATE (JDIV(0:ILEVEL,MDIV),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR JDIV ARRAY',ERROR
ENDIF
!
DDIV=(XB(2)-XB(1))/DFLOAT(NDIV)
NCELLS=NDIV
WRITE (9,*) 'The number of sampling cells is',NCELLS
NCCELLS=NCIS*NDIV
WRITE (9,*) 'The number of collision cells is',NCCELLS
!
!IF (IFX == 0) XS=0.D00  !--isebasti: do not reset XS read from READ_DATA
!
IF (ISECS == 0) THEN
  IF (IFX == 0) FNUM=(XB(2)-XB(1))*FND(1)/DFLOAT(NMI)
  IF (IFX == 1) FNUM=PI*(XB(2)**2-XB(1)**2)*FND(1)/DFLOAT(NMI)
  IF (IFX == 2) FNUM=1.3333333333333333333333D00*PI*(XB(2)**3-XB(1)**3)*FND(1)/DFLOAT(NMI)
ELSE
  IF (IFX == 0) FNUM=((XS-XB(1))*FND(1)+(XB(2)-XS)*FND(2))/DFLOAT(NMI)
  IF (IFX == 1) FNUM=PI*((XS**2-XB(1)**2)*FND(1)+(XB(2)**2-XS**2)*FND(2))/DFLOAT(NMI)
  IF (IFX == 2) FNUM=1.3333333333333333333333D00*PI*((XS**3-XB(1)**3)*FND(1)-(XB(2)**3-XS**3)*FND(2))/DFLOAT(NMI)
END IF
!
FTIME=0.D00
TLIM=1.D20
!
TOTMOV=0.D00
TOTCOL=0.D00
PCOLLS=0.D00
TOTDUP=0.D00
TCOL=0.D00
CSCR=0.D00
TDISS=0.D00
TRECOMB=0.D00
TFOREX=0.D00
TREVEX=0.D00
TREACG=0
TREACL=0
TNEX=0.D00
!
DO N=1,NDIV
  JDIV(0,N)=-N
END DO
!
ALLOCATE (CELL(4,NCELLS),ICELL(NCELLS),CCELL(5,NCCELLS),ICCELL(3,NCCELLS),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR CELL ARRAYS',ERROR
ENDIF
!
ALLOCATE (COLLS(NCELLS),WCOLLS(NCELLS),CLSEP(NCELLS),REAC(MNRE),SREAC(MNSR),VAR(21,NCELLS), &
          VARSP(0:11,NCELLS,MSP),VARS(0:32+MSP,2),CS(0:8+MMVM,NCELLS,MSP),CSS(0:8,2,MSP,2), &  !-isebasti: correcting CS allocation
          CSSS(6,2),CST(0:4,NCELLS),BINS(0:NBINS,5,MSP),BIN(0:NBINS,5),PDFS(0:NBINS,5,MSP),PDF(0:NBINS,5), &
          NDROT(MSP,100),NDVIB(NSCELLS,0:MMVM,MSP,0:100),STAT=ERROR) !--isebasti: CST,PDFS,PDF,BINS,BIN included
ALLOCATE (EVREM(MNRE),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR SAMPLING ARRAYS',ERROR
ENDIF
!
!
CALL INITIALISE_SAMPLES
!
!--Set the initial cells

DO N=1,NCELLS
  CELL(2,N)=XB(1)+DFLOAT(N-1)*DDIV
  CELL(3,N)=CELL(2,N)+DDIV
  CELL(1,N)=CELL(2,N)+0.5D00*DDIV
  IF (IFX == 0) CELL(4,N)=CELL(3,N)-CELL(2,N)    !--calculation assumes unit cross-section
  IF (IFX == 1) CELL(4,N)=PI*(CELL(3,N)**2-CELL(2,N)**2)  !--assumes unit length of full cylinder
  IF (IFX == 2) CELL(4,N)=1.33333333333333333333D00*PI*(CELL(3,N)**3-CELL(2,N)**3)    !--flow is in the full sphere
  ICELL(N)=NCIS*(N-1)
  DO M=1,NCIS
    L=ICELL(N)+M
    XMIN=CELL(2,N)+DFLOAT(M-1)*DDIV/DFLOAT(NCIS)
    XMAX=XMIN+DDIV/DFLOAT(NCIS)
    IF (IFX == 0) CCELL(1,L)=XMAX-XMIN
    IF (IFX == 1) CCELL(1,L)=PI*(XMAX**2-XMIN**2)  !--assumes unit length of full cylinder
    IF (IFX == 2) CCELL(1,L)=1.33333333333333333333D00*PI*(XMAX**3-XMIN**3)    !--flow is in the full sphere
    CCELL(2,L)=0.D00
    ICCELL(3,L)=N
  END DO
  VAR(8,N)=FTMP(1)
  VAR(9,N)=FRTMP(1)
  VAR(10,N)=FVTMP(1)
  VAR(11,N)=FTMP(1)
END DO
!
IF (IWF == 0) AWF=1.D00
IF (IWF == 1) THEN
!--FNUM must be reduced to allow for the weighting factors
  A=0.D00
  B=0.D00
  DO N=1,NCELLS
    A=A+CELL(4,N)
    B=B+CELL(4,N)/(1.+WFM*CELL(1,N)**IFX)
  END DO
  AWF=A/B
  FNUM=FNUM*B/A
END IF
!
WRITE (9,*) 'FNUM is',FNUM
!
!--set the information on the molecular species
!
A=0.D00
B=0.D00
DO L=1,MSP
  A=A+SP(5,L)*FSP(L,1)
  B=B+(3.+ISPR(1,L))*FSP(L,1)
  VMP(L,1)=SQRT(2.D00*BOLTZ*FTMP(1)/SP(5,L))
  IF ((ITYPE(2)== 0).OR.(ISECS == 1)) VMP(L,2)=SQRT(2.D00*BOLTZ*FTMP(2)/SP(5,L))
  VNMAX(L)=3.*VMP(L,1)
  IF (L == 1) THEN
    VMPM=VMP(L,1)
  ELSE
    IF (VMP(L,1) > VMPM) VMPM=VMP(L,1)
  END IF
END DO
WRITE (9,*) 'VMPM =',VMPM
FDEN=A*FND(1)
FPR=FND(1)*BOLTZ*FTMP(1)
FMA=VFX(1)/DSQRT((B/(B+2.D00))*BOLTZ*FTMP(1)/A)
DO L=1,MSP
  DO M=1,MSP
    SPM(4,L,M)=0.5D00*(SP(1,L)+SP(1,M))  !d_ref
    SPM(3,L,M)=0.5D00*(SP(3,L)+SP(3,M))  !omega_ref
    SPM(5,L,M)=0.5D00*(SP(2,L)+SP(2,M))  !t_ref
    SPM(1,L,M)=SP(5,L)*(SP(5,M)/(SP(5,L)+SP(5,M)))
    SPM(2,L,M)=0.25D00*PI*(SP(1,L)+SP(1,M))**2
    AA=2.5D00-SPM(3,L,M)
    A=GAM(AA)
    SPM(6,L,M)=1.D00/A
    SPM(8,L,M)=0.5D00*(SP(4,L)+SP(4,M))
    IF ((ISPR(1,L) > 0).AND.(ISPR(1,M) > 0)) THEN
      SPM(7,L,M)=(SPR(1,L)+SPR(1,M))*0.5D00
    END IF
    IF ((ISPR(1,L) > 0).AND.(ISPR(1,M) == 0)) THEN
      SPM(7,L,M)=SPR(1,L)
    END IF
    IF ((ISPR(1,M) > 0).AND.(ISPR(1,L) == 0)) THEN
      SPM(7,L,M)=SPR(1,M)
    END IF
  END DO
END DO
!
IF (GASCODE==8) THEN !special code with wysong2014 data for O2-O collisions
    ! O-O2
    SPM(3,3,4)=0.75d0      !omega_ref
    SPM(4,3,4)=3.442d-10  !d_ref
    SPM(5,3,4)=273.d0     !t_ref
    SPM(2,3,4)=PI*SPM(4,3,4)**2
    SPM(6,3,4)=1.D00/GAM(2.5D0 - SPM(3,3,4))
    DO L = 2,6
      SPM(L,4,3) = SPM(L,3,4)
    ENDDO

    ! N2-O2
    L = 3; M = 1;
    SPM(3,L,M)=0.7343d0      !omega_ref
    SPM(4,L,M)=4.0512d-10  !d_ref
    SPM(5,L,M)=273.d0     !t_ref
    SPM(2,L,M)=PI*SPM(4,L,M)**2
    SPM(6,L,M)=1.D00/GAM(2.5D0 - SPM(3,L,M))
    DO N = 2,6
      SPM(N,M,L) = SPM(N,L,M)
    ENDDO

    ! NO-N2
    L = 1; M = 5;
    SPM(3,L,M)=0.7292554d0      !omega_ref
    SPM(4,L,M)=4.34104d-10  !d_ref
    SPM(5,L,M)=273.d0     !t_ref
    SPM(2,L,M)=PI*SPM(4,L,M)**2
    SPM(6,L,M)=1.D00/GAM(2.5D0 - SPM(3,L,M))
    DO N = 2,6
      SPM(N,M,L) = SPM(N,L,M)
    ENDDO
END IF

!IF (GASCODE == 8 )THEN
  !SPM(3,3,4) = 0.75d0
  !SPM(4,3,4) = 9.61877d-10
  !SPM(5,3,4) = 273.d0
  !SPM(2,3,4) = PI*SPM(4,3,4)**2
  !SPM(3,4,3) = SPM(3,3,4)
  !SPM(4,4,3) = SPM(4,3,4)
  !SPM(5,4,3) = SPM(5,3,4)
  !SPM(2,4,3) = SPM(2,4,3)

  !SPM(3,3,3) = 1.678d0
  !SPM(4,3,3) = 24.926d-10
  !SPM(5,3,3) = 273.d0
  !SPM(2,3,3) = PI*SPM(4,3,3)**2
!END IF



!
IF (MNRE > 0) CALL INIT_REAC
!
IF (MMVM > 0) THEN  !--isebasti: included
  IDL=0
  DO L=1,MSP
    IF (ISPV(L) > 0) THEN
      DO K=1,ISPV(L)
        IDL(K,L)=SPVM(4,K,L)/SPVM(1,K,L)  !assuming SHO
      END DO
    END IF
  END DO
END IF
!
IF (MNRE > 0) THEN
  DO L=1,MSP
    IF (ISPV(L) == 0) DPER(L)=0   !number of different vibrational energy levels
    IF (ISPV(L) == 1) DPER(L)=IDL(K,L)+1.d0
  END DO
!
  CPER(:,:)=1.d0  !CPER
  DO KK=1,MNRE
    IF (REA(5,KK) < 0.d0) THEN !consider only endothermic (e.g., dissocation) reactions
      DO KN=1,2
        L=IREA(KN,KK)
!
        IF (ISPV(L) == 1) THEN
          A=0.d0
          B=0.d0
          DO I=0,IDL(1,L)
            B=B+1.d0
            EVIB=DFLOAT(I)*BOLTZ*SPVM(1,1,L)
            IF (EVIB < DABS(REA(5,KK))) A=A+1.d0
          END DO
          CPER(KK,KN)= DPER(L)/A
        END IF
!
        IF (ISPV(L) == 3) THEN
          II=0
          DO I=0,IDL(1,L)
            DO J=0,IDL(2,L)
             DO K=0,IDL(3,L)
                II=II+1.d0    !number of possible states
                !stores the energy level for vibrational state (I,J,K)
                TVAR(II)=DFLOAT(I)*BOLTZ*SPVM(1,1,L)+DFLOAT(J)*BOLTZ*SPVM(1,2,L)+DFLOAT(K)*BOLTZ*SPVM(1,3,L)
              END DO
            END DO
          END DO
!
          NDES=0
          DO I=1,II
            DO J=1,II
              A=0.999999d0*TVAR(J)
              B=1.000001d0*TVAR(J)
              IF ((TVAR(I) > A ).AND.(TVAR(I) < B )) THEN
                NDES=NDES+1    !number of different energy states
                IDES(NDES)=I   !a different energy state
              END IF
            END DO
          END DO
          DPER(L)=NDES         !number of different vibrational energy levels
!
          A=0.d0
          DO I=1,NDES
            EVIB=TVAR(IDES(I))
            IF (EVIB < DABS(REA(5,KK))) A=A+1.d0  !number of different energy states bellow REA(5,K)
          END DO
          !this is the C constant in Bondar's paper required for the normalization condition (Sum[CPER/DPER]=A*CPER/DPER=1)
          !note that only states with Evib < REA(5,K) contribute to the summation; others states are prohibited
          CPER(KK,KN)=DPER(L)/A
        END IF
        !write(*,*) KK,L,IDL(:,L),II,DPER(L),A,CPER(KK,KN)
      END DO
    END IF
  END DO
!
END IF
!pause
!
IF (MSP == 1) THEN
  RMAS=SPM(1,1,1)
  CXSS=SPM(2,1,1)
  RGFS=SPM(6,1,1)
END IF
!
DO L=1,MSP
  CR(L)=0.
  DO M=1,MSP
    CR(L)=CR(L)+2.D00*SPI*SPM(4,L,M)**2*FND(1)*FSP(M,1)*(FTMP(1)/SPM(5,L,M))** &
        (1.-SPM(3,L,M))*DSQRT(2.*BOLTZ*SPM(5,L,M)/SPM(1,L,M))
  END DO
END DO
A=0.D00
FP=0.D00
DO L=1,MSP
  A=A+FSP(L,1)*CR(L)
  FP=FP+(2./SPI)*FSP(L,1)*VMP(L,1)/CR(L)
END DO
CTM=1.D00/A
WRITE (9,*) 'Approximate collision time in the stream is',CTM
!
WRITE (9,*) 'Approximate mean free path in the stream is',FP
!
!--set the initial time step
DTM=CTM*CPDTM
IF (DABS(VFX(1)) > 1.D-6) THEN
  A=(0.5D00*DDIV/VFX(1))*TPDTM
ELSE
  A=0.5D00*DDIV/VMPM
END IF
IF (IVB == 1) THEN
  B=0.25D00*DDIV/(ABS(VELOB)+VMPM)
  IF (B < A) A=B
END IF
IF (DTM > A) DTM=A
DTSAMP=SAMPRAT*DTM
DTOUT=OUTRAT*DTSAMP
TSAMP=0.D00
TOUT=0.D00
TPOUT=OUTRAT
NOUT=-1
ENTMASS=0.D00
!
WRITE (9,*) 'The initial value of the overall time step is',DTM
!
!--initialise cell quantities associated with collisions
!
DO N=1,NCCELLS
  CCELL(3,N)=DTM/2.D00
  CCELL(4,N)=2.D00*VMPM*SPM(2,1,1)
  CALL ZGF(RANF,IDT)
  CCELL(2,N)=RANF
  CCELL(5,N)=0.D00
END DO

!-- Han: enlarge CCELL(4) in case of error
IF (nonVHS .ne. 0) THEN
  CCELL(4,:) = CCELL(4,:)*1.2D0
END IF
!
!--set the entry quantities
!
DO K=1,2
  IF (ITYPE(K) == 0) THEN
    DO L=1,MSP
      IF (K == 1) SN=VFX(1)/VMP(L,1)
      IF (K == 2) SN=-VFX(2)/VMP(L,2)
      AA=SN
      A=1.D00+ERF(AA)
      BB=DEXP(-SN**2)
      ENTR(3,L,K)=SN
      ENTR(4,L,K)=SN+SQRT(SN**2+2.D00)
      ENTR(5,L,K)=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,K)))
      ENTR(6,L,K)=3.D00*VMP(L,K)
      B=BB+SPI*SN*A
      ENTR(1,L,K)=(FND(K)*FSP(L,K)*VMP(L,K))*B/(FNUM*2.D00*SPI)
      ENTR(2,L,K)=0.D00
    END DO
  END IF
END DO
!
!--initialize variables for updating boundary conditions
!
UVFX(:)=VFX(:)  !--isebasti
UFND(:)=FND(:)
UFTMP(:)=FTMP(:)
UVMP(:,:)=VMP(:,:)
!
!--Set the uniform stream
!
MNM=1.1D00*NMI
!
IF (MMVM > 0) THEN
  ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),     &
       IPVIB(0:MMVM,MNM),STAT=ERROR)
  IPVIB=0
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
NM=0
!
IF ((IGS==1).OR.(IGS==3)) THEN
  WRITE (*,*) 'Setting the initial gas'
  DO L=1,MSP
    ROTE=0.
    DO K=1,ISECS+1
      IF (ISECS == 0) THEN         !--no secondary stream
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)
        XMIN=XB(1)
        XMAX=XB(2)
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (K == 1) THEN
          M=IDINT(DFLOAT(NMI)*((XS**JFX-XB(1)**JFX)*FND(1)/A)*FSP(L,1))
          XMIN=XB(1)
          XMAX=XS
        ELSE
          M=IDINT(DFLOAT(NMI)*((XB(2)**JFX-XS**JFX)*FND(2)/A)*FSP(L,2))
          XMIN=XS
          XMAX=XB(2)
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,VMP(L,K),IDT)
            END DO
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,FRTMP(K),PROT(NM),IDT)
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,FVTMP(K),IPVIB(J,NM),J,IDT)
                  IF (IGS==3) THEN
                    !forcing a specific vibrational distribution
                    CALL ZGF(RANF,IDT)
                    !IF (IPVIB(J,NM)==5  .AND. RANF <= 0.6d-0) IPVIB(J,NM)=10 !Tv at 2930K
                    !IF (IPVIB(J,NM)>=15 .AND. RANF <= 0.2d-0) IPVIB(J,NM)=10 !Tv at 2930K
                    !IF (IPVIB(J,NM)>=5) IPVIB(J,NM)=10 !Tv at 2820K
                    !IF (RANF <= 1d-6) IPVIB(J,NM)=59 !trick
                  END IF
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
    END DO
  END DO
END IF
!
!--isebasti: included IGS=2 option
IF (IGS==2) THEN
  WRITE (*,*) 'Setting the initial gas with given temp TDF(NC) and number density DNF(NC) distributions'
  DO L=1,MSP
    ROTE=0.
    DO NC=1,NCELLS
   !DO K=1,ISECS+1
      XMIN=CELL(2,NC)
      XMAX=CELL(3,NC)
      TDF=FTMP(1)+2.d3*DEXP(-((CELL(1,NC)-XS)/0.2d-3)**2.d0)
      DNF=FND(1)*(FTMP(1)/TDF)     !nden is Gaussian
      TDF=FTMP(1)                  !T is uniform
      IF (ISECS == 0) THEN         !--no secondary stream
        K=1
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)*(CELL(4,NC)/(XB(2)-XB(1)))*(DNF/FND(1))
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (CELL(1,NC) < XS) THEN
          K=1
          M=IDINT(DFLOAT(NMI)*((XMAX**JFX-XMIN**JFX)*DNF/A)*FSP(L,1))
        ELSE
          K=2
          M=IDINT(DFLOAT(NMI)*((XMAX**JFX-XMIN**JFX)*DNF/A)*FSP(L,2))
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            VMP2=DSQRT(2.D00*BOLTZ*TDF/SP(5,L))  !based on TDF
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,VMP2,IDT)
            END DO
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,TDF,PROT(NM),IDT)   !based on TDF
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,TDF,IPVIB(J,NM),J,IDT)    !based on TDF
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
   !END DO
    END DO
  END DO
!
  WRITE (9,*) 'DMOM',DMOM
  WRITE (9,*) 'DENG',DENG
END IF
!
!--isebasti: included IGS=4 option
IF (IGS==4) THEN
  WRITE (*,*) 'Setting the initial gas with bimodal energy distribution'
  DO L=1,MSP
    ROTE=0.
    DO K=1,ISECS+1
      IF (ISECS == 0) THEN         !--no secondary stream
        M=(DFLOAT(NMI)*FSP(L,1)*AWF)
        XMIN=XB(1)
        XMAX=XB(2)
      ELSE
        A=(XS**JFX-XB(1)**JFX)*FND(1)+(XB(2)**JFX-XS**JFX)*FND(2)
        IF (K == 1) THEN
          M=IDINT(DFLOAT(NMI)*((XS**JFX-XB(1)**JFX)*FND(1)/A)*FSP(L,1))
          XMIN=XB(1)
          XMAX=XS
        ELSE
          M=IDINT(DFLOAT(NMI)*((XB(2)**JFX-XS**JFX)*FND(2)/A)*FSP(L,2))
          XMIN=XS
          XMAX=XB(2)
        END IF
      END IF
      IF ((K == 1).OR.(ISECS == 1)) THEN
        III=0
        WFMIN=1.D00+WFM*XB(1)**IFX
        N=1
        INC=1
        DO WHILE (N < M)
          A=(XMIN**JFX+((DFLOAT(N)-0.5D00)/DFLOAT(M))*(XMAX-XMIN)**JFX)**(1.D00/DFLOAT(JFX))
          IF (IWF == 0) THEN
            B=1.D00
          ELSE
            B=WFMIN/(1.D00+WFM*A**IFX)
            IF ((B < 0.1D00).AND.(INC == 1)) INC=10
            IF ((B < 0.01D00).AND.(INC == 10)) INC=100
            IF ((B < 0.001D00).AND.(INC == 100)) INC=1000
            IF ((B < 0.0001D00).AND.(INC == 1000)) INC=10000
          END IF
          CALL ZGF(RANF,IDT)
          IF (B*DFLOAT(INC) > RANF) THEN
            NM=NM+1
            PX(NM)=A
            IPSP(NM)=L
            PTIM(NM)=0.
            IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),KK)
            IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),KK,PTIM(NM))
            !--for bimodal distribution
            CALL ZGF(RANF,IDT)
            IF (RANF > 0.5d0) THEN
              BTEMP=FTMP(1)   !based on translational reference temperature
            ELSE
              BTEMP=FRTMP(1)  !based on rotational reference temperature
            END IF
            BVMP=SQRT(2.D00*BOLTZ*BTEMP/SP(5,L))
            DO KK=1,3
              CALL RVELC(PV(KK,NM),A,BVMP,IDT) !for bimodal distribution
            END DO
            !--end
            PV(1,NM)=PV(1,NM)+VFX(K)
            PV(2,NM)=PV(2,NM)+VFY(K)
            IF (ISPR(1,L) > 0) CALL SROT(L,BTEMP,PROT(NM),IDT) !for bimodal distribution
            IF (MMVM > 0) THEN
              IPVIB(0,NM)=0
              IF (ISPV(L) > 0) THEN
                DO J=1,ISPV(L)
                  CALL SVIB(L,FVTMP(K),IPVIB(J,NM),J,IDT) !for Boltzmann distribution
                END DO
              END IF
            END IF
          END IF
          N=N+INC
        END DO
      END IF
    END DO
  END DO
END IF
!
!--SPECIAL CODING FOR INITIATION OF EXPLOSION IN H2-02 MIXTURE (FORCED IGNITION CASES)
IF ((GASCODE == 8).AND.(IFI > 0)) THEN !--isebasti: IFC>0 allow forced ignition
  IF(IGS == 1) A=0.5D00 !0.5 for p=0.1 atm;
  IF(IGS == 2) A=1.0D00
  !
  !--set the vibrational levels of A% random molecules to 5
  IF (IFI == 1) THEN
    M=0.01D00*A*NM
    DO N=1,M
      CALL ZGF(RANF,IDT)
      K=INT(RANF*DFLOAT(NM))+1
      L=IPSP(K)       !--isebasti: modified to increase all vibrational degrees of freedom
      DO KK=1,ISPV(L) !--isebasti
        IPVIB(KK,K)=5 !--isebasti
      END DO          !--isebasti
    END DO
  END IF
  !
  !--isebasti: set vibrational levels of molecules within XMIN and XMAX to dissociation level
  IF (IFI == 2) THEN
    XMIN=XS-0.005D00*A*(XB(2)-XB(1))
    XMAX=XS+0.005D00*A*(XB(2)-XB(1))
    DO N=1,NM
      IF ((PX(N) >= XMIN).AND.(PX(N) <= XMAX)) THEN
        L=IPSP(N)
        IF (ISPV(L) > 0) IPVIB(1,N)=IDL(1,L) !setting only one mode
      END IF
    END DO
  END IF
  !
  !--isebasti: set molecules within XMIN and XMAX to energy states consistent with temperature TDF
  IF (IFI == 3) THEN
    XMIN=XS-0.005D00*A*(XB(2)-XB(1))
    XMAX=XS+0.005D00*A*(XB(2)-XB(1))
    DO N=1,NM
      IF ((PX(N) >= XMIN).AND.(PX(N) <= XMAX)) THEN
        L=IPSP(N)
        TDF=5.D03
        B=SQRT(2.D00*BOLTZ*TDF/SP(5,L))  !B=VMP at temperature TDF
        DO K=1,3
          CALL RVELC(PV(K,N),BB,B,IDT)
        END DO
        IF (PX(N) < XS) PV(1,N)=PV(1,N)+VFX(1)
        IF (PX(N) > XS) PV(1,N)=PV(1,N)+VFX(2)
        IF (ISPR(1,L) > 0) CALL SROT(L,TDF,PROT(N),IDT)
        IF (MMVM > 0) THEN
          DO K=1,ISPV(L)
            CALL SVIB(L,TDF,IPVIB(K,N),K,IDT)
          END DO
        END IF
      END IF
    END DO
END IF
!
END IF

CLOSE(9)
!
!
RETURN
END SUBROUTINE SET_INITIAL_STATE
!
!***************************************************************************
!******************************GAS DATABASE*********************************
!***************************************************************************
!
SUBROUTINE HARD_SPHERE
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=1
MMRM=0
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.0D-10    !reference diameter
SP(2,1)=273.       !reference temperature
SP(3,1)=0.5        !viscosity-temperature index
SP(4,1)=1.         !reciprocal of VSS scattering parameter (1 for VHS)
SP(5,1)=5.D-26     !mass
ISPR(1,1)=0        !number of rotational degrees of freedom
!
RETURN
END SUBROUTINE HARD_SPHERE
!
!***************************************************************************
!
SUBROUTINE ARGON
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=1
MMRM=0
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
SP(1,1)=4.11D-10   !--isebasti: use VSS; VHS value is 4.17D-10
SP(2,1)=273.15
SP(3,1)=0.81
SP(4,1)=1.D0/1.4D0 !--isebasti: use VSS
SP(5,1)=6.63D-26
ISPR(1,1)=0
ISPR(2,1)=0
!
RETURN
END SUBROUTINE ARGON
!
!***************************************************************************
!
SUBROUTINE IDEAL_GAS
!
USE GAS
USE CALC
!
IMPLICIT NONE
INTEGER :: J !--isebasti: J included
!
IF(IRM == 1) THEN
  MSP=1; MMRM=2; MMVM=0
END IF
!
IF(IRM >= 2) THEN
  MSP=1; MMRM=2; MMVM=1
END IF
!
IF(IRM == 4) THEN
  MSP=1; MMRM=3; MMVM=3
END IF
!
IF(IRM == 5) THEN
  MSP=2; MMRM=3; MMVM=3
END IF
!
IF(IRM == 6) THEN
  MSP=2; MMRM=3; MMVM=3
END IF
!
MNRE=0
MTBP=0
MEX=0
MMEX=0
MNSR=0
!
CALL ALLOCATE_GAS
!
IF(IRM == 1) THEN !N2 (Trans+Rot; Zr const)
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
END IF
!
IF(IRM == 2) THEN !N2 (Trans+Rot+Vib; Zv const)
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=5.D00
SPVM(3,1,1)=-1.D00
SPVM(4,1,1)=113500.D00
END IF
!
IF(IRM == 3) THEN !N2 (Trans+Rot+Vib; Zv(T))
SP(1,1)=4.17D-10
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=52600.D00
SPVM(3,1,1)=3371.D00
SPVM(4,1,1)=113500.D00
END IF
!
IF(IRM == 4) THEN !H2O (Trans+Rot+Vib; Zv(T))
J=1
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=0.75D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
END IF
!
IF(IRM == 5) THEN !H2O+H2O (Trans+Rot+Vib; Zv const)
J=1
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.0D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
!
J=2
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.0D00        !-estimate from SMILE
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=5.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=5.D00 !20000.D00   !--estimate
SPVM(3,1,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=5.D00 !20000.D00   !--estimate
SPVM(3,2,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=5.D00 !20000.D00   !--estimate
SPVM(3,3,J)=-1.D00 !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
END IF
!
IF(IRM == 6) THEN !H2O+H (Trans+Rot+Vib)
J=1  !--species 1 is water vapor H2O
SP(1,J)=4.5D-10       !--estimate
SP(2,J)=273.D00
SP(3,J)=1.00D00      !-estimate
SP(4,J)=1.0D00
SP(5,J)=2.99D-26
ISPR(1,J)=3
ISPR(2,J)=0
SPR(1,J)=1.D00
ISPV(J)=3
SPVM(1,1,J)=5261.D00  !--symmetric stretch mode
SPVM(2,1,J)=10. !20000.D00   !--estimate
SPVM(3,1,J)=-1. !2500.D00    !--estimate
SPVM(4,1,J)=60043.83D00
SPVM(1,2,J)=2294.D00  !--bend mode
SPVM(2,2,J)=10. !20000.D00   !--estimate
SPVM(3,2,J)=-1. !2500.D00    !--estimate
SPVM(4,2,J)=60043.83D00
SPVM(1,3,J)=5432.D00  !--asymmetric stretch mode
SPVM(2,3,J)=10. !20000.D00   !--estimate
SPVM(3,3,J)=-1. !2500.D00    !--estimate
SPVM(4,3,J)=60043.83D00
!
J=2  !--species 2 is atomic hydrogen H
SP(1,J)=2.5D-10      !--estimate
SP(2,J)=273.D00
SP(3,J)=0.8D00
SP(4,J)=1.D00
SP(5,J)=1.67D-27
ISPR(1,J)=0
ISPV(J)=0
END IF
!
RETURN
END SUBROUTINE IDEAL_GAS
!
!***************************************************************************
!
SUBROUTINE REAL_OXYGEN
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=1
MMVM=1
MNRE=0
MTBP=0
MEX=0
MMEX=0
MNSR=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.07D-10
SP(2,1)=273.D00
SP(3,1)=0.77D00
SP(4,1)=1.D00
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0             ! 0,1 for constant,polynomial rotational relaxation collision number
SPR(1,1)=5.             ! the collision number or the coefficient of temperature in the polynomial (if a polynomial, the coeff. of T^2 is in spr_db(3  )
ISPV(1)=1               ! the number of vibrational modes
SPVM(1,1,1)=2256.D00          ! the characteristic vibrational temperature
SPVM(2,1,1)=90000.D00        ! a constant Zv, or the reference Zv
IF (IZV == 1) SPVM(2,1,1)=5.D00*SPVM(2,1,1)   !--to allow for the reduction in relaxation time when based on quantized collision temperature
SPVM(3,1,1)=2256.D00        ! -1 for a constant Zv, or the reference temperature
SPVM(4,1,1)=59500.D00       ! characteristic dissociation temperature
ISPVM(1,1,1)=2
ISPVM(2,1,1)=2
!
!--species 2 is atomic oxygen
SP(1,2)=3.D-10
SP(2,2)=273.D00
SP(3,2)=0.8D00
SP(4,2)=1.D00
SP(5,2)=2.656D-26
ISPR(1,2)=0
ISPV(2)=0     !--must be set!
!
!--set data needed for recombination
ISPRC=0
ISPRC(2,2)=1              !--O+O -> O2  recombined species code for an O+O recombination
SPRC(1,2,2,1)=0.04D00     !--UNADJUSTED
IF (ITCV == 1) SPRC(1,2,2,1)=SPRC(1,2,2,1)*0.23
SPRC(2,2,2,1)=-1.3D00
SPRC(1,2,2,2)=0.07D00     !--UNADJUSTED
IF (ITCV == 1) SPRC(1,2,2,2)=SPRC(1,2,2,2)*0.28
SPRC(2,2,2,2)=-1.2D00
!
NSPEX=0
SPEX=0.D00
ISPEX=0
!
RETURN
END SUBROUTINE REAL_OXYGEN
!
!***************************************************************************
!
SUBROUTINE IDEAL_AIR
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=1
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=4.07D-10
SP(2,1)=273.
SP(3,1)=0.77
SP(4,1)=1.
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.
SP(1,2)=4.17D-10
SP(2,2)=273.
SP(3,2)=0.74
SP(4,2)=1.
SP(5,2)=4.65D-26
ISPR(1,2)=2
ISPR(2,2)=0
SPR(1,2)=5.
RETURN
END SUBROUTINE IDEAL_AIR
!
!***************************************************************************
!
SUBROUTINE REAL_AIR
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=5
MMRM=1
MMVM=1

!
IF (JCD == 0) THEN
  MNRE=23
  MTBP=4
  MEX=0
  MMEX=0
END IF
IF (JCD == 1) THEN
  MNRE=0
  MTBP=0
  MEX=4
  MMEX=1
END IF
!
MNSR=0
CALL ALLOCATE_GAS
!--species 1 is oxygen
SP(1,1)=4.07D-10
SP(2,1)=273.D00
SP(3,1)=0.77D00
SP(4,1)=1.d00
SP(5,1)=5.312D-26
ISPR(1,1)=2
ISPR(2,1)=0
SPR(1,1)=5.D00
ISPV(1)=1               ! the number of vibrational modes
SPVM(1,1,1)=2256.D00          ! the characteristic vibrational temperature
SPVM(2,1,1)=18000.D00  !90000.D00        ! a constant Zv, or the reference Zv
SPVM(3,1,1)=2256.D00        ! -1 for a constant Zv, or the reference temperature
SPVM(4,1,1)=59500.D00
ISPVM(1,1,1)=3
ISPVM(2,1,1)=3
!--species 2 is nitrogen
SP(1,2)=4.17D-10
SP(2,2)=273.D00
SP(3,2)=0.74D00
SP(4,2)=1.D00
SP(5,2)=4.65D-26
ISPR(1,2)=2
ISPR(2,2)=0
SPR(1,2)=5.D00
ISPV(2)=1
SPVM(1,1,2)=3371.D00
SPVM(2,1,2)=52000.D00     !260000.D00
SPVM(3,1,2)=3371.D00
SPVM(4,1,2)=113500.D00
ISPVM(1,1,2)=4
ISPVM(2,1,2)=4
!--species 3 is atomic oxygen
SP(1,3)=3.D-10
SP(2,3)=273.D00
SP(3,3)=0.8D00
SP(4,3)=1.D00
SP(5,3)=2.656D-26
ISPR(1,3)=0
ISPV(3)=0
!--species 4 is atomic nitrogen
SP(1,4)=3.D-10
SP(2,4)=273.D00
SP(3,4)=0.8D00
SP(4,4)=1.0D00
SP(5,4)=2.325D-26
ISPR(1,4)=0
ISPV(4)=0
!--species 5 is NO
SP(1,5)=4.2D-10
SP(2,5)=273.D00
SP(3,5)=0.79D00
SP(4,5)=1.0D00
SP(5,5)=4.98D-26
ISPR(1,5)=2
ISPR(2,5)=0
SPR(1,5)=5.D00
ISPV(5)=1
SPVM(1,1,5)=2719.D00
SPVM(2,1,5)=14000.D00   !70000.D00
SPVM(3,1,5)=2719.D00
SPVM(4,1,5)=75500.D00
ISPVM(1,1,5)=3
ISPVM(2,1,5)=4
!--following data is required if JCD=1 (new reaction model)
IF (JCD == 1) THEN
!--set the recombination data for the molecule pairs
  ISPRC=0    !--data os zero unless explicitly set
  SPRC=0.D00
  ISPRC(3,3)=1    !--O+O -> O2  recombined species code for an O+O recombination
  SPRC(1,3,3,1)=0.04D00
  SPRC(2,3,3,1)=-1.3D00
  SPRC(1,3,3,2)=0.07D00
  SPRC(2,3,3,2)=-1.2D00
  SPRC(1,3,3,3)=0.08D00
  SPRC(2,3,3,3)=-1.2D00
  SPRC(1,3,3,4)=0.09D00
  SPRC(2,3,3,4)=-1.2D00
  SPRC(1,3,3,5)=0.065D00
  SPRC(2,3,3,5)=-1.2D00
  ISPRC(4,4)=2
  SPRC(1,4,4,1)=0.15D00
  SPRC(2,4,4,1)=-2.05D00
  SPRC(1,4,4,2)=0.09D00
  SPRC(2,4,4,2)=-2.1D00
  SPRC(1,4,4,3)=0.16D00
  SPRC(2,4,4,3)=-2.0D00
  SPRC(1,4,4,4)=0.17D00
  SPRC(2,4,4,4)=-2.0D00
  SPRC(1,4,4,5)=0.17D00
  SPRC(2,4,4,5)=-2.1D00
  ISPRC(3,4)=5
  SPRC(1,3,4,1)=0.3D00
  SPRC(2,3,4,1)=-1.9D00
  SPRC(1,3,4,2)=0.4D00
  SPRC(2,3,4,2)=-2.0D00
  SPRC(1,3,4,3)=0.3D00
  SPRC(2,3,4,3)=-1.75D00
  SPRC(1,3,4,4)=0.3D00
  SPRC(2,3,4,4)=-1.75D00
  SPRC(1,3,4,5)=0.15D00
  SPRC(2,3,4,5)=-1.9D00
  ISPRC(4,3)=5
  SPRC(1,4,3,1)=0.3D00
  SPRC(2,4,3,1)=-1.9D00
  SPRC(1,4,3,2)=0.4D00
  SPRC(2,4,3,2)=-2.0D00
  SPRC(1,4,3,3)=0.3D00
  SPRC(2,4,3,3)=-1.75D00
  SPRC(1,4,3,4)=0.3D00
  SPRC(2,4,3,4)=-1.75D00
  SPRC(1,4,3,5)=0.15D00
  SPRC(2,4,3,5)=-1.9D00
!--set the exchange reaction data
  SPEX=0.D00
  ISPEX=0
  NSPEX=0
  NSPEX(2,3)=1
  NSPEX(3,2)=1
  NSPEX(5,4)=1
  NSPEX(4,5)=1
  NSPEX(5,3)=1
  NSPEX(3,5)=1
  NSPEX(1,4)=1
  NSPEX(4,1)=1
!--N2+O->NO+N
  ISPEX(1,0,2,3)=2
  ISPEX(1,1,2,3)=5
  ISPEX(1,2,2,3)=4
  ISPEX(1,3,2,3)=1
  ISPEX(1,4,2,3)=1
  SPEX(1,1,2,3)=0.15D00
  SPEX(2,1,2,3)=0.D00
  SPEX(3,1,2,3)=-5.175D-19
  NEX(1,2,3)=1
  ISPEX(1,0,3,2)=2
  ISPEX(1,1,3,2)=5
  ISPEX(1,2,3,2)=4
  ISPEX(1,3,3,2)=1
  ISPEX(1,4,3,2)=1
  SPEX(1,1,3,2)=0.15D00
  SPEX(2,1,3,2)=0.D00
  SPEX(3,1,3,2)=-5.175D-19
  NEX(1,3,2)=1
!--NO+N->N2+0
  ISPEX(1,0,5,4)=5
  ISPEX(1,1,5,4)=2
  ISPEX(1,2,5,4)=3
  ISPEX(1,3,5,4)=1
  ISPEX(1,4,5,4)=1
  SPEX(1,1,5,4)=0.033D00
  SPEX(2,1,5,4)=0.8D00
  SPEX(3,1,5,4)=5.175D-19
  NEX(1,5,4)=2
  ISPEX(1,0,4,5)=5
  ISPEX(1,1,4,5)=2
  ISPEX(1,2,4,5)=3
  ISPEX(1,3,4,5)=1
  ISPEX(1,4,4,5)=1
  SPEX(1,1,4,5)=0.033D00
  SPEX(2,1,4,5)=0.8D00
  SPEX(3,1,4,5)=5.175D-19
  NEX(1,4,5)=2
!--NO+0->O2+N
  ISPEX(1,0,5,3)=5
  ISPEX(1,1,5,3)=1
  ISPEX(1,2,5,3)=4
  ISPEX(1,3,5,3)=1
  ISPEX(1,4,5,3)=1
  SPEX(1,1,5,3)=0.05D00
  SPEX(2,1,5,3)=0.7D00
  SPEX(3,1,5,3)=-2.719D-19
  NEX(1,5,3)=3
  ISPEX(1,0,3,5)=5
  ISPEX(1,1,3,5)=1
  ISPEX(1,2,3,5)=4
  ISPEX(1,3,3,5)=1
  ISPEX(1,4,3,5)=1
  SPEX(1,1,3,5)=0.05D00
  SPEX(2,1,3,5)=0.7D00
  SPEX(3,1,3,5)=-2.719D-19
  NEX(1,3,5)=3
!--O2+N->NO+O
  ISPEX(1,0,1,4)=1
  ISPEX(1,1,1,4)=5
  ISPEX(1,2,1,4)=3
  ISPEX(1,3,1,4)=1
  ISPEX(1,4,1,4)=1
  SPEX(1,1,1,4)=0.D00
  SPEX(2,1,1,4)=0.D00
  SPEX(3,1,1,4)=2.719D-19
  NEX(1,1,4)=4
  ISPEX(1,0,4,1)=1
  ISPEX(1,1,4,1)=5
  ISPEX(1,2,4,1)=3
  ISPEX(1,3,4,1)=1
  ISPEX(1,4,4,1)=1
  SPEX(1,1,4,1)=0.D00
  SPEX(2,1,4,1)=0.D00
  SPEX(3,1,4,1)=2.719D-19
  NEX(1,4,1)=4
!
END IF
IF (JCD == 0) THEN
!--the following data is required if JCD=0 (old reaction model)
! REACTION 1 IS O2+N->2O+N
  LE(1)=1
  ME(1)=4
  KP(1)=3
  LP(1)=3
  MP(1)=4
  CI(1)=1.
  AE(1)=8.197D-19
  AC(1)=5.993D-12
  BC(1)=-1.
  ER(1)=-8.197D-19
!--REACTION 2 IS O2+NO>2O+NO
  LE(2)=1
  ME(2)=5
  KP(2)=3
  LP(2)=3
  MP(2)=5
  CI(2)=1.
  AE(2)=8.197D-19
  AC(2)=5.993D-12
  BC(2)=-1.
  ER(2)=-8.197D-19
!--REACTION 3 IS O2+N2>2O+N2
LE(3)=1
ME(3)=2
KP(3)=3
LP(3)=3
MP(3)=2
CI(3)=1.5
AE(3)=8.197D-19
AC(3)=1.198D-11
BC(3)=-1.
ER(3)=-8.197D-19
!--REACTION 4 IS 2O2>2O+O2
LE(4)=1
ME(4)=1
KP(4)=3
LP(4)=3
MP(4)=1
CI(4)=1.5
AE(4)=8.197D-19
AC(4)=5.393D-11
BC(4)=-1.
ER(4)=-8.197D-19
!--REACTION 5 IS O2+O>3O
LE(5)=1
ME(5)=3
KP(5)=3
LP(5)=3
MP(5)=3
CI(5)=1.
AE(5)=8.197D-19
AC(5)=1.498D-10
BC(5)=-1.
ER(5)=-8.197D-19
!--REACTION 6 IS N2+O>2N+O
LE(6)=2
ME(6)=3
KP(6)=4
LP(6)=4
MP(6)=3
CI(6)=0.5
AE(6)=1.561D-18
AC(6)=3.187D-13
BC(6)=-0.5
ER(6)=-1.561D-18
!--REACTION 7 IS N2+O2>2N+O2
LE(7)=2
ME(7)=1
KP(7)=4
LP(7)=4
MP(7)=1
CI(7)=0.5
AE(7)=1.561D-18
AC(7)=3.187D-13
BC(7)=-0.5
ER(7)=-1.561D-18
!--REACTION 8 IS N2+NO>2N+NO
LE(8)=2
ME(8)=5
KP(8)=4
LP(8)=4
MP(8)=5
CI(8)=0.5
AE(8)=1.561D-18
AC(8)=3.187D-13
BC(8)=-0.5
ER(8)=-1.561D-18
!--REACTION 9 IS 2N2>2N+N2
LE(9)=2
ME(9)=2
KP(9)=4
LP(9)=4
MP(9)=2
CI(9)=1.
AE(9)=1.561D-18
AC(9)=7.968D-13
BC(9)=-0.5
ER(9)=-1.561D-18
!--REACTION 10 IS N2+N>3N
LE(10)=2
ME(10)=4
KP(10)=4
LP(10)=4
MP(10)=4
CI(10)=1.
AE(10)=1.561D-18
AC(10)=6.9E12
BC(10)=-1.5
ER(10)=-1.561D-18
!--REACTION 11 IS NO+N2>N+O+N2
LE(11)=5
ME(11)=2
KP(11)=4
LP(11)=3
MP(11)=2
CI(11)=1.
AE(11)=1.043D-18
AC(11)=6.59D-10
BC(11)=-1.5
ER(11)=-1.043D-18
!--REACTION 12 IS NO+O2>N+O+O2
LE(12)=5
ME(12)=1
KP(12)=4
LP(12)=3
MP(12)=1
CI(12)=1.
AE(12)=1.043D-18
AC(12)=6.59D-10
BC(12)=-1.5
ER(12)=-1.043D-18
!--REACTION 13 IS NO+NO>N+O+NO
LE(13)=5
ME(13)=5
KP(13)=4
LP(13)=3
MP(13)=5
CI(13)=1.
AE(13)=1.043D-18
AC(13)=1.318D-8
BC(13)=-1.5
ER(13)=-1.043D-18
!--REACTION 14 IS NO+O>N+O+O
LE(14)=5
ME(14)=3
KP(14)=4
LP(14)=3
MP(14)=3
CI(14)=1.
AE(14)=1.043D-18
AC(14)=1.318D-8
BC(14)=-1.5
ER(14)=-1.043D-18
!--REACTION 15 IS NO+N>2N+O
LE(15)=5
ME(15)=4
KP(15)=4
LP(15)=3
MP(15)=4
CI(15)=1.
AE(15)=1.043D-18
AC(15)=1.318D-8
BC(15)=-1.5
ER(15)=-1.043D-18
!--REACTION 16 IS NO+O>O2+N
LE(16)=5
ME(16)=3
KP(16)=0
LP(16)=1
MP(16)=4
CI(16)=0.
AE(16)=2.719D-19
AC(16)=5.279D-21
BC(16)=1.
ER(16)=-2.719D-19
!--REACTION 17 IS N2+O>NO+N
LE(17)=2
ME(17)=3
KP(17)=0
LP(17)=5
MP(17)=4
CI(17)=0.
AE(17)=5.175D-19
AC(17)=1.120D-16
BC(17)=0.
ER(17)=-5.175D-19
!--REACTION 18 IS O2+N>NO+O   !Exothermic Exchange
LE(18)=1
ME(18)=4
KP(18)=0
LP(18)=5
MP(18)=3
CI(18)=0.
AE(18)=4.968D-20
AC(18)=1.598D-18
BC(18)=0.5
ER(18)=2.719D-19
!--REACTION 19 IS NO+N>N2+O   !Exothermic Exchange
LE(19)=5
ME(19)=4
KP(19)=0
LP(19)=2
MP(19)=3
CI(19)=0.
AE(19)=0.
AC(19)=2.49D-17
BC(19)=0.
ER(19)=5.175D-19
!--REACTION 20 IS O+O+M1>O2+M1
LE(20)=3
ME(20)=3
KP(20)=-1
LP(20)=1
MP(20)=-1
CI(20)=0.
AE(20)=0.
AC(20)=8.297D-45
BC(20)=-0.5
ER(20)=8.197D-19
!--REACTION 21 IS N+N+M2>N2+M2
LE(21)=4
ME(21)=4
KP(21)=-1
LP(21)=2
MP(21)=-2
CI(21)=0.
AE(21)=0.
AC(21)=3.0051D-44
BC(21)=-0.5
ER(21)=1.561D-18
!--REACTION 22 IS N+N+N>N2+N
LE(22)=4
ME(22)=4
KP(22)=-1
LP(22)=2
MP(22)=-3
CI(22)=0.
AE(22)=0.
AC(22)=6.3962D-40
BC(22)=-1.5
ER(22)=1.5637D-18
!--REACTION 23 IS N+O+M3>NO+M3
LE(23)=4
ME(23)=3
KP(23)=-1
LP(23)=5
MP(23)=-4
CI(23)=0.
AE(23)=0.
AC(23)=2.7846D-40
BC(23)=-1.5
ER(23)=1.043D-18
!
  THBP(1,1)=9.
  THBP(1,2)=2.
  THBP(1,3)=25.
  THBP(1,4)=1.
  THBP(1,5)=1.
  THBP(2,1)=1.
  THBP(2,2)=2.5
  THBP(2,3)=1.
  THBP(2,4)=0.
  THBP(2,5)=1.
  THBP(3,1)=0.
  THBP(3,2)=0.
  THBP(3,3)=0.
  THBP(3,4)=1.
  THBP(3,5)=0.
  THBP(4,1)=1.
  THBP(4,2)=1.
  THBP(4,3)=20.
  THBP(4,4)=20.
  THBP(4,5)=20.
END IF
RETURN
END SUBROUTINE REAL_AIR
!
!******************************************************************************************************************************************
!
SUBROUTINE HELIUM_XENON
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
MSP=2
MMRM=0
MMVM=0
MNRE=0
MTBP=0
MNSR=0
MEX=0
MMEX=0
!
CALL ALLOCATE_GAS
!
SP(1,1)=2.33D-10
SP(2,1)=273.
SP(3,1)=0.66
SP(4,1)=1.
SP(5,1)=6.65D-27
ISPR(1,1)=0
ISPR(2,1)=0
SP(1,2)=5.74D-10
SP(2,2)=273.
SP(3,2)=0.85
SP(4,2)=1.
SP(5,2)=21.8D-26
ISPR(1,2)=0
ISPR(2,2)=0
RETURN
END SUBROUTINE HELIUM_XENON
!
!***************************************************************************
!
SUBROUTINE OXYGEN_NITROGEN
!
USE GAS
USE CALC
USE MOLECS
!
IMPLICIT NONE
REAL(KIND=8) :: VDC=-1.d0, ZV=1.d12 !1.d12
INTEGER :: I,J,ZCHECK,IRATE=0,JRATE=0
CHARACTER (LEN=16) :: FILENAME
!--VDC>=0 for VDC model; <0 for tce model
!--IRATE=0 uncorrected rates; 1 corrected rates for recombination and exchange reactions
!--JRATE= same as IRATE but for dissociation reactions
!
MSP=8
MMRM=3
MMVM=3
!
IF (JCD == 1) THEN
  MNRE=0
  MTBP=0
  MEX=16
  MMEX=3
END IF
IF (JCD == 0) THEN
  IF(IRM == 0) MNRE=0     !Reactions are not included
  IF(IRM == 1) MNRE=37    !Shatalov
  IF(IRM == 2) MNRE=34    !Davidenko
  IF(IRM >= 100) MNRE=1   !testing only one reaction
  IF(IRM >= 200) MNRE=2   !testing only two reactions
  IF(IRM >= 400) MNRE=4   !testing only four reactions
  MTBP=0
  MEX=0
  MMEX=0
END IF
!
MNSR=0
!
CALL ALLOCATE_GAS
!
!--species 1 is hydrogen H2
!SP(1,1)=2.92D-10        !reference diameter  !vss/vhs  2.88/2.92
!SP(2,1)=273.D00         !reference temperature
!SP(3,1)=0.67D00         !viscosity index
!SP(4,1)=1.d00           !reciprocal of VSS scattering parameter; vss 1/1.35 (with BUG!)
!SP(5,1)=3.34D-27        !molecular mass
!ISPR(1,1)=2             !number of rotational degrees of freedom of species L
!ISPR(2,1)=0             !0,1 for constant, polynomial rotational relaxations collision number (Zrot)
!SPR(1,1)=5.d0           !=100 from Haas1995 !constant Zrot or constant in 2nd order polynomial in temperature
!ISPV(1)=1               !the number of vibrational modes
!SPVM(1,1,1)=6159.D00    !the characteristic vibrational temperature
!SPVM(2,1,1)=ZV  !18000.D00   !constant Zv, or reference Zv for mode K --estimate
!SPVM(3,1,1)=-1. !6159.D00    !-1 for constant Zv, or reference temperature
!SPVM(4,1,1)=52438.76D00 !charactheristic dissociation temperature; based on heats of formation
!ISPVM(1,1,1)=2          !the species code of the 1st dissociation product
!ISPVM(2,1,1)=2          !the species code of the 2nd dissociation product
!--species 2 is atomic hydrogen H
!SP(1,2)=2.33D-10         !SMILE/estimate 2.33/2.5
!SP(2,2)=273.D00
!SP(3,2)=0.75D00          !SMILE/estimate 0.75/0.80
!SP(4,2)=1.D00
!SP(5,2)=1.67D-27
!ISPR(1,2)=0
!ISPV(2)=0
!--species 1 is nitrogen N2
SP(1,1)=4.17d-10        !same as in DS1V
SP(2,1)=273.D00
SP(3,1)=0.74D00
SP(4,1)=1.D00
SP(5,1)=4.65D-26
ISPR(1,1)=2
ISPR(2,1)=0             !1 for polinomial Zr, 0 for constant Zr
SPR(1,1)=1.d0 !5.D00
ISPV(1)=1
SPVM(1,1,1)=3371.D00
SPVM(2,1,1)=ZV  !52000.D00 !260000.D00
SPVM(3,1,1)=-1.  !3371.D00
SPVM(4,1,1)=113500.D00
ISPVM(1,1,1)=2
ISPVM(2,1,1)=2
!--species 2 is atomic nitrogen N
SP(1,2)=3.d-10          !same as inDS1V
SP(2,2)=273.D00
SP(3,2)=0.8D00
SP(4,2)=1.0D00
SP(5,2)=2.325D-26
ISPR(1,2)=0
ISPV(2)=0
!--species 3 is oxygen O2
!SP(1,3)=3.985d-10       !wysong2014/vss/vhs 3.985/4.01/4.07
!SP(2,3)=273.D00         !t_ref
!SP(3,3)=0.71D00         !wysong2014

SP(1,3)=4.1515d-10     !Han correct by fitting esposito's data
SP(2,3)=273.D0
SP(3,3)=0.7318d0

SP(4,3)=1.d0            !vss 1./1.4d0 !there is some bug with vss model
SP(5,3)=5.312D-26       !vss correction in collision subroutine must be tested
ISPR(1,3)=2
ISPR(2,3)=0             !1 for polinomial Zr, 0 for constant Zr
SPR(1,3)=1.0 !5.d0           !constant Zr
ISPV(3)=1               !the number of vibrational modes
SPVM(1,1,3)=2273.54     !obtained from nist 2256.D00    !the characteristic vibrational temperature
SPVM(2,1,3)=ZV  !18000.D00   !constant Zv, or the reference Zv
SPVM(3,1,3)=-1. !2256.d0     !-1 for a constant Zv, or the reference temperature
SPVM(4,1,3)=59971.4D00  !the characteristic dissociation temperature
ISPVM(1,1,3)=4
ISPVM(2,1,3)=4
!--species 4 is atomic oxygen O
SP(1,4)=3.458d-10        !wysong2014
SP(2,4)=273.D00          !t_ref
SP(3,4)=0.76D00          !wysong2014
SP(4,4)=1.D00
SP(5,4)=2.656D-26
ISPR(1,4)=0
SPR(1,4)=5.d0
ISPV(4)=0
!--species 5 is hydroxy OH
!SP(1,5)=3.50d-10          !--SMILE/estimate 3.50/4.10
!SP(2,5)=273.D00
!SP(3,5)=0.75D00         !-estimate
!SP(4,5)=1.0D00
!SP(5,5)=2.823D-26
!ISPR(1,5)=2
!ISPR(2,5)=0
!SPR(1,5)=5.D00
!ISPV(5)=1
!SPVM(1,1,5)=5360.D00
!SPVM(2,1,5)=ZV  !18000.D00   !--estimate
!SPVM(3,1,5)=-1. !5360.D00    !--estimate
!SPVM(4,1,5)=51497.18D00
!ISPVM(1,1,5)=2
!ISPVM(2,1,5)=4
!--species 5 is nitric oxide NO
SP(1,5)=4.2D-10        !same as in DS1V
SP(2,5)=273.D00
SP(3,5)=0.79D00
SP(4,5)=1.0D00
SP(5,5)=4.981D-26
ISPR(1,5)=2
ISPR(2,5)=0
SPR(1,5)=5.D00
ISPV(5)=1
SPVM(1,1,5)=2719.D00
SPVM(2,1,5)=ZV  !14000.D00   !70000.D00
SPVM(3,1,5)=-1. !2719.D00
SPVM(4,1,5)=75500.D00
ISPVM(1,1,5)=2
ISPVM(2,1,5)=4
!--species 6 is water vapor H2O
SP(1,6)=4.5d-10         !--estimate
SP(2,6)=273.D00
SP(3,6)=1.0D00         !SMILE/estimate 1.0/0.75
SP(4,6)=1.0D00
SP(5,6)=2.99D-26
ISPR(1,6)=3
ISPR(2,6)=0
SPR(1,6)=10.d0        !--estimate Alexeenko (2003)
ISPV(6)=3
SPVM(1,1,6)=5261.D00           !--symmetric stretch mode
SPVM(2,1,6)=10.d0 !ZV !10.d0   !--estimate Alexeenko (2003); for same species collisions
SPVM(3,1,6)=-2.                !-2 indicates that Zv is different for same species collisions
SPVM(4,1,6)=60043.83D00
SPVM(1,2,6)=2294.D00  !--bend mode
SPVM(2,2,6)=250.d0 !ZV !250.d0    !--estimate Alexeenko (2003); for different species collisions
SPVM(3,2,6)=SPVM(3,1,6)
SPVM(4,2,6)=60043.83D00
SPVM(1,3,6)=5432.D00   !--asymmetric stretch mode
SPVM(2,3,6)=SPVM(2,2,6) !ZV !250.d0
SPVM(3,3,6)=SPVM(3,2,6)
SPVM(4,3,6)=60043.83D00
ISPVM(1,1,6)=2
ISPVM(2,1,6)=5
ISPVM(1,2,6)=2
ISPVM(2,2,6)=5
ISPVM(1,3,6)=2
ISPVM(2,3,6)=5
!--species 7 is hydroperoxy HO2
SP(1,7)=5.5d-10       !--estimate
SP(2,7)=273.D00
SP(3,7)=0.75D00      !-estimate
SP(4,7)=1.0D00
SP(5,7)=5.479D-26
ISPR(1,7)=2    !--assumes that HO2 is linear
ISPR(2,7)=0
SPR(1,7)=10.D00   !estimate
ISPV(7)=3
SPVM(1,1,7)=4950.D00
SPVM(2,1,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,1,7)=-1. ! 2500.D00    !--estimate
SPVM(4,1,7)=24988.08D00
SPVM(1,2,7)=2000.D00
SPVM(2,2,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,2,7)=-1. !2500.D00    !--estimate
SPVM(4,2,7)=24988.08D00
SPVM(1,3,7)=1580.D00
SPVM(2,3,7)=ZV !250 !20000.D00   !--estimate
SPVM(3,3,7)=-1. !2500.D00    !--estimate
SPVM(4,3,7)=24988.08D00
ISPVM(1,1,7)=2
ISPVM(2,1,7)=3
ISPVM(1,2,7)=2
ISPVM(2,2,7)=3
ISPVM(1,3,7)=2
ISPVM(2,3,7)=3
!--species 8 is hydroperoxy CO2
SP(1,8)=5.687766d-10       !SMILE
SP(2,8)=273.d0
SP(3,8)=0.93d0
SP(4,8)=1.0d0
SP(5,8)=7.309999d-26
ISPR(1,8)=2    !--assumes that CO2 is linear
ISPR(2,8)=0
SPR(1,8)=10.d0 !estimate
ISPV(8)=3
SPVM(1,1,8)=945.d0
SPVM(2,1,8)=ZV !250.d0    !--estimate; for different species collisions
SPVM(3,1,8)=-1.       !-2 indicates that Zv is different for same species collisions
SPVM(4,1,8)=64015.d0
SPVM(1,2,8)=1903.d0
SPVM(2,2,8)=ZV !10.d0     !--for same species collisions
SPVM(3,2,8)=-1.
SPVM(4,2,8)=64015.d0
SPVM(1,3,8)=3329.d0
SPVM(2,3,8)=ZV !10.d0
SPVM(3,3,8)=-1.
SPVM(4,3,8)=64015.d0
ISPVM(1,1,8)=2
ISPVM(2,1,8)=3
ISPVM(1,2,8)=2
ISPVM(2,2,8)=3
ISPVM(1,3,8)=2
ISPVM(2,3,8)=3
!--Species 8 is argon
!SP(1,8)=4.17D-10  !vss/vhs 4.11/4.17
!SP(2,8)=273.15
!SP(3,8)=0.81
!SP(4,8)=1.        !vss 1./1.4
!SP(5,8)=6.63D-26
!ISPR(1,8)=0
!ISPV(8)=0
!
!ISPV(:)=0 !test mechanism without vibrational modes
!
IF (JCD ==1) THEN    !--Q-K model
!all the lines were deletd --isebasti
END IF
!
!--the following data is required if JCD=0 (TCE reaction model)
!
IF (JCD == 0 .AND. IRM == 2) THEN
!--rate data is based on Davidenko2006_SystematicStudySupersonicCombustion
!--------------------------------------------------------------
!--REACTION 1 IS H+H+H2>2H2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
  IREV(J)=7
!--REACTION 2 IS H+H+H>H2+H
  J=2
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
  IREV(J)=8
!--REACTION 3 IS H+H+O2>H2+O2
  J=3
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
  IREV(J)=9
!--REACTION 4 IS H+H+O>H2+O
  J=4
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
  IREV(J)=10
!--REACTION 5 IS H+H+OH>H2+OH
  J=5
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
  IREV(J)=11
!--REACTION 6 IS H+H+H2O>H2+H2O
  J=6
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
  IREV(J)=12
!--------------------------------------------------------------
! REACTION 7 IS H2+H2->2H+H2
  J=7
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
!--REACTION 8 IS H2+H>2H+H
  J=8
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=2
!--REACTION 9 IS H2+O2>2H+O2
  J=9
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=3
!--REACTION 10 IS H2+O>2H+O
  J=10
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=4
!--REACTION 11 IS H2+OH>2H+OH
  J=11
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=5
!--REACTION 12 IS H2+H2O>2H+H2O
  J=12
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=6
!--------------------------------------------------------------
!--REACTION 13 IS OH+H+H2>H2O+H2
  J=13
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
  IREV(J)=19
!--REACTION 14 IS OH+H+H>H2O+H
  J=14
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
  IREV(J)=20
!--REACTION 15 IS OH+H+O2>H2O+O2
  J=15
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
  IREV(J)=21
!--REACTION 16 IS OH+H+O>H2O+O
  J=16
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
  IREV(J)=22
!--REACTION 17 IS OH+H+OH>H2O+OH
  J=17
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
  IREV(J)=23
!--REACTION 18 IS OH+H+H2O>H2O+H2O
  J=18
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
  IREV(J)=24
!--------------------------------------------------------------
!--REACTION 19 IS H2O+H2>OH+H+H2
  J=19
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=13
!--REACTION 20 IS H2O+H>OH+H+H
  J=20
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=14
!--REACTION 21 IS H2O+O2>OH+H+O2
  J=21
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=15
!--REACTION 22 IS H2O+O>OH+H+O
  J=22
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=16
!--REACTION 23 IS H2O+OH>OH+H+OH
  J=23
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=17
!--REACTION 24 IS H2O+H2O>OH+H+H2O
  J=24
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=18
!--------------------------------------------------------------
!--REACTION 25 IS OH+OH>H2O+O
  J=25
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
  IREV(J)=26
!--REACTION 26 IS H2O+O>OH+OH
  J=26
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
  IREV(J)=25
!--------------------------------------------------------------
!--REACTION 27 IS H2+O>OH+H
  J=27
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
  IREV(J)=28
!--REACTION 28 IS OH+H>H2+O
  J=28
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
  IREV(J)=27
!--------------------------------------------------------------
!--REACTION 29 IS H2+OH>H20+H
  J=29
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
  IREV(J)=30
!--REACTION 30 IS H2O+H>H2+OH
  J=30
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ  !Dav 9300.d0  !Shat 9030.d0  !Conaire 9252.d0
  AC(J)=7.964d8/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  BC(J)=1.528d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  IF (IRATE == 1) BC(J)=1.939508d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
  IREV(J)=29
!--------------------------------------------------------------
!--REACTION 31 IS O2+H>OH+O
  J=31
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
  IREV(J)=32
!--REACTION 32 IS OH+O>O2+H
  J=32
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
  IREV(J)=31
!--------------------------------------------------------------
!--REACTION 33 IS H2+O2>OH+OH
  J=33
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
  IREV(J)=34
!--REACTION 34 IS OH+OH>H2+O2
  J=34
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
  IREV(J)=33
END IF
!
!--------------------------------------------------------------
!Testing individual reactions
J=1
!
IF (JCD == 0 .AND. IRM == 101) THEN
!--REACTION D1 IS H+H+H2>2H2  (Davidenko)
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 102) THEN
!--REACTION D2 IS H+H+H>H2+H
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 103) THEN
!--REACTION D3 IS H+H+O2>H2+O2
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 104) THEN
!--REACTION D4 IS H+H+O>H2+O
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 105) THEN
!--REACTION D5 IS H+H+OH>H2+OH
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 106) THEN
!--REACTION D6 IS H+H+H2O>H2+H2O
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 107) THEN
! REACTION D7 IS H2+H2->2H+H2
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 108) THEN
!--REACTION D8 IS H2+H>2H+H
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 109) THEN
!--REACTION D9 IS H2+O2>2H+O2
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 110) THEN
!--REACTION D10 IS H2+O>2H+O
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 111) THEN
!--REACTION D11 IS H2+OH>2H+OH
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 112) THEN
!--REACTION D12 IS H2+H2O>2H+H2O
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
END IF
!
IF (JCD == 0 .AND. IRM == 113) THEN
!--REACTION D13 IS OH+H+H2>H2O+H2
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 114) THEN
!--REACTION D14 IS OH+H+H>H2O+H
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 115) THEN
!--REACTION D15 IS OH+H+O2>H2O+O2
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 116) THEN
!--REACTION D16 IS OH+H+O>H2O+O
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 117) THEN
!--REACTION D17 IS OH+H+OH>H2O+OH
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 118) THEN
!--REACTION D18 IS OH+H+H2O>H2O+H2O
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 119) THEN
!--REACTION D19 IS H2O+H2>OH+H+H2
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 120) THEN
!--REACTION D20 IS H2O+H>OH+H+H
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 121) THEN
!--REACTION D21 IS H2O+O2>OH+H+O2
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 122) THEN
!--REACTION D22 IS H2O+O>OH+H+O
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 123) THEN
!--REACTION D23 IS H2O+OH>OH+H+OH
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 124) THEN
!--REACTION D24 IS H2O+H2O>OH+H+H2O
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
END IF
!
IF (JCD == 0 .AND. IRM == 125) THEN
!--REACTION D25 IS OH+OH>H20+O
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
END IF
!
IF (JCD == 0 .AND. IRM == 126) THEN
!--REACTION D26 IS H2O+O>OH+OH
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
END IF
!
IF (JCD == 0 .AND. IRM == 127) THEN
!--REACTION D27 IS H2+O>OH+H
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
END IF
!
IF (JCD == 0 .AND. IRM == 128) THEN
!--REACTION D28 IS OH+H>H2+O
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
END IF
!
IF (JCD == 0 .AND. IRM == 129) THEN
!--REACTION D29 IS H2+OH>H20+H
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
END IF
!
IF (JCD == 0 .AND. IRM == 130) THEN
!--REACTION D30 IS H2O+H>H2+OH
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ
  AC(J)=7.964d8/(AVOG*1000.d0)
  BC(J)=1.528d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.939508d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
END IF
!
IF (JCD == 0 .AND. IRM == 131) THEN
!--REACTION D31 IS O2+H>OH+O
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
END IF
!
IF (JCD == 0 .AND. IRM == 132) THEN
!--REACTION D32 IS OH+O>O2+H
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
END IF
!
IF (JCD == 0 .AND. IRM == 133) THEN
!--REACTION D33 IS H2+O2>OH+OH
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
END IF
!
IF (JCD == 0 .AND. IRM == 134) THEN
!--REACTION D34 IS OH+OH>H2+O2
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
END IF
!
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 140) THEN
!--REACTION 140 IS N2+N>N+N+N (from QCT data)
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=113200.d0*BOLTZ         !from AHO model
  AC(J)=6.6831d-6*1.0d-6        ! fit from Jaffe data
  BC(J)=-0.6996d0
  ER(J)=-AE(J)
END IF
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 141) THEN
!--REACTION 140 IS N2+N2>N2+N+N (from Bender data)
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=117000.d0*BOLTZ         !from AHO model
  AC(J)=4.5d-6*1.0d-6        ! fit from Jaffe data
  BC(J)=-0.675
  ER(J)=-AE(J)
END IF
!--------------------------------------------------------------
!
IF (JCD == 0 .AND. IRM == 150) THEN
!--REACTION 150 IS O2+O>O+O+O (from QCT data)
  LE(J)=3
  ME(J)=4
  KP(J)=4
  LP(J)=4
  MP(J)=4
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT          !from AHO model
  AC(J)=2.500d18/(AVOG*1000.d0)  !AC and BC are fit from QCT data
  BC(J)=-0.565d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN
      AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-0.185509d0
    ELSE
      AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-0.185509d0
    END IF
  END IF
  ER(J)=-5.21275d0*EVOLT !based on heat of formation !-5.21275d0*EVOLT !from AHO model
END IF
!
IF (JCD == 0 .AND. IRM == 160) THEN
!--REACTION 160 IS O2+O2>O+O+O2 (fit from park1990 but using the DSMC value of activation energy)
  LE(J)=3
  ME(J)=3
  KP(J)=4
  LP(J)=4
  MP(J)=3
  CI(J)=VDC
  !AE(J)=5.21275d0*EVOLT           !from AHO model
  !AC(J)=1.9337d22/(AVOG*1000.d0)  !AC and BC from park1990
  !BC(J)=-1.7334d0
  AE(J) = 59380.d0*BOLTZ+(5.21275d0-5.1153d0)*EVOLT
   !Han use Byron's rate but shift to add zeropoint energy
  AC(J) = 1.9392865d21/(AVOG*1000.d0)
  BC(J) = -1.5d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN
      AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-1.1919d0
    ELSE
      AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
      BC(J)=-1.1919d0
    END IF
  END IF
  !ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 161) THEN
!--REACTION 160 IS O2+N2>O+O+N2 (fit from park1990 but using the DSMC value of activation energy)
  LE(J)=3
  ME(J)=1
  KP(J)=4
  LP(J)=4
  MP(J)=1
  CI(J)=VDC
  AE(J) = 5.21275d0*EVOLT
  AC(J) = 1.586878535d15/(AVOG*1000.d0)
  BC(J) = 4.255d-3
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 162) THEN
!--REACTION 160 IS N2+O2>N+N+O2
  LE(J) = 1
  ME(J) = 3
  KP(J) = 2
  LP(J) = 2
  MP(J) = 3
  CI(J) = VDC
  !  Wray: doi:10.2514/5.2018-0240 fit good with Andrienko's QCT
  ! here we refit to match dissociation energy
  AC(J) = 1.1447510d18/(AVOG*1000.d0)
  BC(J) = -0.6823d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 163) THEN
!--REACTION 160 IS NO+N2>N+O+N2
  LE(J) = 5
  ME(J) = 1
  KP(J) = 2
  LP(J) = 4
  MP(J) = 1
  CI(J) = VDC
  ! Park's rate, match well with Andrienko's QCT
  ! 5.0e15, 0 , 75500K
  ! here we refit to match dissociation energy
  AC(J) = 1.08705682d16/(AVOG*1000.d0)
  BC(J) = -0.0766715d0
  AE(J) = 6.55879*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 164) THEN
!--REACTION 160 IS N2+NO>N+N+NO
  LE(J) = 1
  ME(J) = 5
  KP(J) = 2
  LP(J) = 2
  MP(J) = 5
  CI(J) = VDC
  ! Park uses N2+N2 as N2+NO
  ! Wray suggest N2+Ar as N2+NO
  ! here we refit Wray's data to match dissociation energy
  AC(J) = 9.0558360d17/(AVOG*1000.d0)
  BC(J) = -0.662299d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 170) THEN
!--REACTION 170 IS N2+O>2N+O (from QCT data)
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
 ! AE(J)=9.8216d0*EVOLT            !from AHO model
 ! AC(J)=8.9330d16/(AVOG*1000.d0)  !AC and BC are fit from QCT data
 ! BC(J)=-0.3842d0
  !IF (JRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)  !to be corrected
  !IF (JRATE == 1) BC(J)=0.7079989d0                 !to be corrected
 ! ER(J)=-9.8216d0*EVOLT !from AHO model
  AE(J) = 113950.d0*BOLTZ
  AC(J) = 8.934d-13
  BC(J) = -0.4807
  ER(J) = -AE(J)
END IF
!
IF (JCD == 0 .AND. IRM == 180) THEN
!--REACTION 180 IS N2+O>NO+N (from QCT data)
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=VDC
  AE(J)=3.2628d0*EVOLT           !from AHO model
  AC(J)=2.4444d11/(AVOG*1000.d0) !AC and BC are fit from QCT data
  BC(J)=0.7071d0
  !IF (JRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00) !to be corrected
  !IF (JRATE == 1) BC(J)=0.7079989d0                !to be corrected
  ER(J)=-3.2628d0*EVOLT !from AHO model
END IF
!--------------------------------------------------------------
!Testing reaction pairs
!
IF (JCD == 0 .AND. IRM == 201) THEN
!--REACTION 1 IS H+H+H2>2H2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=1.180305d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5999194d0
  ER(J)=7.239220d-19
  IREV(J)=2
! REACTION 7 IS H2+H2->2H+H2
  J=2
  LE(J)=1
  ME(J)=1
  KP(J)=2
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=2.5d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.258131d+12/(AVOG*1000.d0)
      BC(J)=8.798179d-01
    ELSE
      AC(J)=1.124184d19/(AVOG*1000.d0)
      BC(J)=-1.097974d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 202) THEN
!--REACTION 2 IS H+H+H>H2+H
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.673070d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5987639d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 8 IS H2+H>2H+H
  J=2
  LE(J)=1
  ME(J)=2
  KP(J)=2
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=4.897980d+13/(AVOG*1000.d0)
      BC(J)=4.199751d-01
    ELSE
      AC(J)=1.799959d20/(AVOG*1000.d0)
      BC(J)=-1.571824d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 203) THEN
!--REACTION 3 IS H+H+O2>H2+O2
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.763155d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6015061d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 9 IS H2+O2>2H+O2
  J=2
  LE(J)=1
  ME(J)=3
  KP(J)=2
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.199534d+11/(AVOG*1000.d0)
      BC(J)=8.911783d-01
    ELSE
      AC(J)=1.115394d20/(AVOG*1000.d0)
      BC(J)=-1.519741d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 204) THEN
!--REACTION 4 IS H+H+O>H2+O
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.636160d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5978748d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 10 IS H2+O>2H+O
  J=2
  LE(J)=1
  ME(J)=4
  KP(J)=2
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.282898d+14/(AVOG*1000.d0)
      BC(J)=2.374446d-01
    ELSE
      AC(J)=2.110156d20/(AVOG*1000.d0)
      BC(J)=-1.565169d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 205) THEN
!--REACTION 5 IS H+H+OH>H2+OH
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=4.720037d16/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.6002319d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 11 IS H2+OH>2H+OH
  J=2
  LE(J)=1
  ME(J)=5
  KP(J)=2
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.260180d+11/(AVOG*1000.d0)
      BC(J)=1.045312d+00
    ELSE
      AC(J)=6.346795d19/(AVOG*1000.d0)
      BC(J)=-1.443632d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 206) THEN
!--REACTION 6 IS H+H+H2O>H2+H2O
  J=1
  LE(J)=2
  ME(J)=2
  KP(J)=-1
  LP(J)=1
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*9.791d16/(AVOG*1000.d0)**2
  BC(J)=-0.6d0
  IF (IRATE == 1) AC(J)=5.648955d17/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-0.5993142d0
  ER(J)=7.239220d-19
  IREV(J)=2
!--REACTION 12 IS H2+H2O>2H+H2O
  J=2
  LE(J)=1
  ME(J)=6
  KP(J)=2
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=52105.d0*BOLTZ
  AC(J)=12.d0*5.086d16/(AVOG*1000.d0)
  BC(J)=-0.362d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.378291d+11/(AVOG*1000.d0)
      BC(J)=1.322471d+00
    ELSE
      AC(J)=5.874078d14/(AVOG*1000.d0)
      BC(J)=0.1337761d0
    END IF
  END IF
  ER(J)=-7.239220d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 213) THEN
!--REACTION 13 IS OH+H+H2>H2O+H2
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=1
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.5d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.670047d22/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.956666d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 19 IS H2O+H2>OH+H+H2
  J=2
  LE(J)=6
  ME(J)=1
  KP(J)=5
  LP(J)=2
  MP(J)=1
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=2.5d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=6.027961d+15/(AVOG*1000.d0)
      BC(J)=1.632129d-01
    ELSE
      AC(J)=3.292285d33/(AVOG*1000.d0)
      BC(J)=-4.719732d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 214) THEN
!--REACTION 14 IS OH+H+H>H2O+H
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.354999d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932729d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 20 IS H2O+H>OH+H+H
  J=2
  LE(J)=6
  ME(J)=2
  KP(J)=5
  LP(J)=2
  MP(J)=2
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.031394d+13/(AVOG*1000.d0)
      BC(J)=6.817159d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !1.900933d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.858611d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 215) THEN
!--REACTION 15 IS OH+H+O2>H2O+O2
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=3
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.314137d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.931931d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 21 IS H2O+O2>OH+H+O2
  J=2
  LE(J)=6
  ME(J)=3
  KP(J)=5
  LP(J)=2
  MP(J)=3
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=9.889330d+15/(AVOG*1000.d0)
      BC(J)=6.684782d-03
    ELSE
      AC(J)=9.336914d32/(AVOG*1000.d0)
      BC(J)=-4.683550d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 216) THEN
!--REACTION 16 IS OH+H+O>H2O+O
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.347222d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.932318d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 22 IS H2O+O>OH+H+O
  J=2
  LE(J)=6
  ME(J)=4
  KP(J)=5
  LP(J)=2
  MP(J)=4
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=2.637473d+18/(AVOG*1000.d0)
      BC(J)=-6.616960d-01
    ELSE
      AC(J)=8.936d22/(AVOG*1000.d0)  !9.055880d33/(AVOG*1000.d0)
      BC(J)=-1.835d0  !-4.942958d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 217) THEN
!--REACTION 17 IS OH+H+OH>H2O+OH
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=5
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=8.422456d21/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-1.933605d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 23 IS H2O+OH>OH+H+OH
  J=2
  LE(J)=6
  ME(J)=5
  KP(J)=5
  LP(J)=2
  MP(J)=5
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=1.430367d+14/(AVOG*1000.d0)
      BC(J)=4.977865d-01
    ELSE
      AC(J)=1.457072d33/(AVOG*1000.d0)
      BC(J)=-4.735936d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 218) THEN
!--REACTION 18 IS OH+H+H2O>H2O+H2O
  J=1
  LE(J)=5
  ME(J)=2
  KP(J)=-1
  LP(J)=6
  MP(J)=6
  CI(J)=0.
  AE(J)=0.d0
  AC(J)=12.d0*2.212d22/(AVOG*1000.d0)**2
  BC(J)=-2.d0
  IF (IRATE == 1) AC(J)=2.412574d23/(AVOG*1000.d0)**2
  IF (IRATE == 1) BC(J)=-2.019927d0
  ER(J)=8.252333d-19
  IREV(J)=2
!--REACTION 24 IS H2O+H2O>OH+H+H2O
  J=2
  LE(J)=6
  ME(J)=6
  KP(J)=5
  LP(J)=2
  MP(J)=6
  CI(J)=VDC
  AE(J)=59743.d0*BOLTZ
  AC(J)=12.d0*8.936d22/(AVOG*1000.d0)
  BC(J)=-1.835d0
  IF (JRATE == 1) THEN
    IF (CI(J) < 0.) THEN !TCE
      AC(J)=3.472299d+15/(AVOG*1000.d0)
      BC(J)=4.256304d-01
    ELSE
      AC(J)=2.735897d22/(AVOG*1000.d0)
      BC(J)=-1.554205d0
    END IF
  END IF
  ER(J)=-8.252333d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 225) THEN
!--REACTION 25 IS OH+OH>H20+O
  J=1
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=4
  CI(J)=0.
  AE(J)=50.d0*BOLTZ
  AC(J)=1.506d9/(AVOG*1000.d0)
  BC(J)=1.14d0
  IF (IRATE == 1) AC(J)=1.499219d9/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.140777d0
  ER(J)=1.113715d-19
  IREV(J)=2
!--REACTION 26 IS H2O+O>OH+OH
  J=2
  LE(J)=6
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=8613.d0*BOLTZ
  AC(J)=2.220d10/(AVOG*1000.D00)
  BC(J)=1.089d0
  IF (IRATE == 1) AC(J)=9.868296d8/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=1.451492d0
  ER(J)=-1.113715d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 227) THEN
!--REACTION 27 IS H2+O>OH+H
  J=1
  LE(J)=1
  ME(J)=4
  KP(J)=0
  LP(J)=5
  MP(J)=2
  CI(J)=0.
  AE(J)=3163.d0*BOLTZ
  AC(J)=5.119d4/(AVOG*1000.d0)
  BC(J)=2.67d0
  IF (IRATE == 1) AC(J)=2.733975d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.760115d0
  ER(J)=-1.006026d-20   !Bird -0.13d-19
  IREV(J)=2
!--REACTION 28 IS OH+H>H2+O
  J=2
  LE(J)=5
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=4
  CI(J)=0.
  AE(J)=2240.d0*BOLTZ
  AC(J)=2.701d4/(AVOG*1000.d0)
  BC(J)=2.65d0
  IF (IRATE == 1) AC(J)=2.465206d4/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=2.659041d0
  ER(J)=1.006026d-20   !Bird 0.13d-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 229) THEN
!--REACTION 29 IS H2+OH>H20+H
  J=1
  LE(J)=1
  ME(J)=5
  KP(J)=0
  LP(J)=6
  MP(J)=2
  CI(J)=0.
  AE(J)=1660.d0*BOLTZ
  AC(J)=1.024d8/(AVOG*1000.d0)
  BC(J)=1.60d0
  IF (IRATE == 1) AC(J)=1.050693d8/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=1.594793d0
  ER(J)=1.013113d-19  !Bird 1.05D-19
  IREV(J)=2
!--REACTION 30 IS H2O+H>H2+OH
  J=2
  LE(J)=6
  ME(J)=2
  KP(J)=0
  LP(J)=1
  MP(J)=5
  CI(J)=0.
  AE(J)=9300.d0*BOLTZ  !Dav 9300.d0  !Shat 9030.d0  !Conaire 9252.d0
  AC(J)=7.964d8/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  BC(J)=1.528d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  IF (IRATE == 1) AC(J)=2.265672d7/(AVOG*1000.d0)  !Dav 7.964d8   !Shat 4.52d8   !Conaire 2.3d9
  IF (IRATE == 1) BC(J)=1.939508d0  !Dav 1.528d0  !Shat 1.6d0  !Conaire 1.4d0
  ER(J)=-1.013113d-19  !Bird 1.05D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 231) THEN
!--REACTION 31 IS O2+H>OH+O
  J=1
  LE(J)=3
  ME(J)=2
  KP(J)=0
  LP(J)=5
  MP(J)=4
  CI(J)=0.
  AE(J)=8456.d0*BOLTZ
  AC(J)=1.987d14/(AVOG*1000.D00)
  BC(J)=0.d0
  IF (IRATE == 1) AC(J)=4.243092d13/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.1832357d0
  ER(J)=-1.137477d-19  !Bird -1.17D-19
  IREV(J)=2
!--REACTION 32 IS OH+O>O2+H
  J=2
  LE(J)=5
  ME(J)=4
  KP(J)=0
  LP(J)=3
  MP(J)=2
  CI(J)=0.
  AE(J)=0.  !--isebasti: paper provides -118.d0*BOLTZ
  AC(J)=8.930d11/(AVOG*1000.d0)
  BC(J)=0.338d0
  IF (IRATE == 1) AC(J)=8.889175d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3388361d0
  ER(J)=1.137477d-19  !Bird 1.17D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 233) THEN
!--REACTION 33 IS H2+O2>OH+OH
  J=1
  LE(J)=1
  ME(J)=3
  KP(J)=0
  LP(J)=5
  MP(J)=5
  CI(J)=0.
  AE(J)=24044.d0*BOLTZ
  AC(J)=1.70d13/(AVOG*1000.d0)
  BC(J)=0.
  IF (IRATE == 1) AC(J)=5.357596d11/(AVOG*1000.d0)
  IF (IRATE == 1) BC(J)=0.3949028d0
  ER(J)=-1.238079d-19 !Bird -1.3D-19
  IREV(J)=2
!--REACTION 34 IS OH+OH>H2+O2
  J=2
  LE(J)=5
  ME(J)=5
  KP(J)=0
  LP(J)=1
  MP(J)=3
  CI(J)=0.
  AE(J)=14554.d0*BOLTZ
  AC(J)=4.032d10/(AVOG*1000.D00)
  BC(J)=0.317d0
  IF (IRATE == 1) AC(J)=1.390438d9/(AVOG*1000.D00)
  IF (IRATE == 1) BC(J)=0.7079989d0
  ER(J)=1.238079d-19 !Bird 1.3D-19
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 250) THEN
!--REACTION 150 IS O2+O>O+O+O (from QCT data)
  J=1
  LE(J)=3
  ME(J)=4
  KP(J)=4
  LP(J)=4
  MP(J)=4
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT          !from AHO model
  AC(J)=2.500d18/(AVOG*1000.d0)  !AC and BC are fit from QCT data
  BC(J)=-0.565d0
  IF (IRATE == 1) AC(J)=5.64764d16/(AVOG*1000.d0) !TCE-corrected
  IF (IRATE == 1) BC(J)=-0.185509d0
  ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  IREV(J)=2
!--REACTION 160 IS O2+O2>O+O+O2 (fit from park1990 but using the DSMC value of activation energy)
  J=2
  LE(J)=3
  ME(J)=3
  KP(J)=4
  LP(J)=4
  MP(J)=3
  CI(J)=VDC
  AE(J)=5.21275d0*EVOLT           !from AHO model
  AC(J)=1.9337d22/(AVOG*1000.d0)  !AC and BC from park1990
  BC(J)=-1.7334d0
  IF (IRATE == 1) AC(J)=9.3069d19/(AVOG*1000.d0) !TCE-corrected
  IF (IRATE == 1) BC(J)=-1.1919d0
  ER(J)=-8.276d-19 !based on heat of formation !-5.21275d0*EVOLT !from AHO model
  IREV(J)=1
END IF
!--------------------------------------------------------------
IF (JCD == 0 .AND. IRM == 251) THEN
!--REACTION 1 is O2+N2>O+O+N2
  LE(1) = 3
  ME(1) = 1
  KP(1) = 4
  LP(1) = 4
  MP(1) = 1
  CI(1) = VDC
  ! Park 2E21cc/mol/sec -1.5 59700
  ! Boyd 8.132E-9 cc/sec -0.131 59380
  ! fit1 2.7212E14 cc/mol/sec 0.22006 59380
  ! fit2 1.5869E15 cc/mol/sec 0.004255  60490
  AC(1) = 1.586878535d15/(AVOG*1000.d0)
  BC(1) = 4.255d-3
  AE(1) = 5.21275d0*EVOLT
  ER(1) = -AE(1)
!--REACTION 2 is N2+O2->N+N+O2
  LE(2) = 1
  ME(2) = 3
  KP(2) = 2
  LP(2) = 2
  MP(2) = 3
  CI(2) = VDC
  !  Wray: doi:10.2514/5.2018-0240
  ! here we refit to match dissociation energy
  AC(2) = 1.1447510d18/(AVOG*1000.d0)
  BC(2) = -0.6823
  AE(2) = 9.82163d0*EVOLT
  ER(2) = -AE(2)
END IF
!
IF (JCD == 0 .AND. IRM == 252) THEN
  !--REACTION 160 IS NO+N2>N+O+N2
  J = 1
  LE(J) = 5
  ME(J) = 1
  KP(J) = 2
  LP(J) = 4
  MP(J) = 1
  CI(J) = VDC
  ! Park's rate, match well with Andrienko's QCT
  ! 5.0e15, 0 , 75500K
  ! here we refit to match dissociation energy
  AC(J) = 1.08705682d16/(AVOG*1000.d0)
  BC(J) = -0.0766715d0
  AE(J) = 6.55879*EVOLT
  ER(J) = -AE(J)
  !--REACTION 160 IS N2+NO>N+N+NO
  J = 2
  LE(J) = 1
  ME(J) = 5
  KP(J) = 2
  LP(J) = 2
  MP(J) = 5
  CI(J) = VDC
  ! Park uses N2+N2 as N2+NO
  ! Wray suggest N2+Ar as N2+NO
  ! here we refit Wray's data to match dissociation energy
  AC(J) = 9.0558360d17/(AVOG*1000.d0)
  BC(J) = -0.662299d0
  AE(J) = 9.82163d0*EVOLT
  ER(J) = -AE(J)
END IF
!--------------------------------------------------------------
!
IF (IPRS == 1) THEN !--read pre-collision vibrational pdfs
  DO I=1,MNRE
    DO J=1,ITMAX
      WRITE(FILENAME,700) I,J
      700 FORMAT('DS1VIBF1_',i2.2,'_',i2.2)
      800 CONTINUE
      OPEN (7,FILE=TRIM(FILENAME)//'.BIN',FORM='UNFORMATTED',ERR=800)
      READ (7) FPTEMP(J),FEVIB(I,J,:,:),FPVIB(I,J,:,:,:),ZCHECK
      CLOSE(7)
      IF (ZCHECK /= 1234567) THEN
        WRITE (9,*) 'Failed reading DS1FVIB for IKA-IT-TEMP; Check:',I,J,FPTEMP(J),ZCHECK
        STOP
      ELSE
        WRITE (9,*) 'Success readding DS1FVIB for IKA-IT-TEMP; Check:',I,J,FPTEMP(J),ZCHECK
      END IF
    END DO
  END DO
END IF
!--------------------------------------------------------------
IF (QCTMODEL > 0) THEN
  VIBEN=1.d3            !initializing variable
!
  J=1                   !N2 vibrational levels from: Han's Luo (Purdue AAE) Master Thesis (2016)
  IVMODEL(J,1)=1        !use pre-computed vibrational levels for this species
  IVMODEL(J,2)=59        !max possible vibrational level
  VIBEN(0,J)=-9.758052d0
  VIBEN(1,J)=-9.467844d0
  VIBEN(2,J)=-9.181767d0
  VIBEN(3,J)=-8.899808d0
  VIBEN(4,J)=-8.621953d0
  VIBEN(5,J)=-8.348193d0
  VIBEN(6,J)=-8.078518d0
  VIBEN(7,J)=-7.812921d0
  VIBEN(8,J)=-7.551397d0
  VIBEN(9,J)=-7.293940d0
  VIBEN(10,J)=-7.040548d0
  VIBEN(11,J)=-6.791219d0
  VIBEN(12,J)=-6.545954d0
  VIBEN(13,J)=-6.304754d0
  VIBEN(14,J)=-6.067621d0
  VIBEN(15,J)=-5.834561d0
  VIBEN(16,J)=-5.605578d0
  VIBEN(17,J)=-5.380680d0
  VIBEN(18,J)=-5.159875d0
  VIBEN(19,J)=-4.943173d0
  VIBEN(20,J)=-4.730586d0
  VIBEN(21,J)=-4.522125d0
  VIBEN(22,J)=-4.317806d0
  VIBEN(23,J)=-4.117644d0
  VIBEN(24,J)=-3.921656d0
  VIBEN(25,J)=-3.729860d0
  VIBEN(26,J)=-3.542278d0
  VIBEN(27,J)=-3.358931d0
  VIBEN(28,J)=-3.179843d0
  VIBEN(29,J)=-3.005039d0
  VIBEN(30,J)=-2.834546d0
  VIBEN(31,J)=-2.668394d0
  VIBEN(32,J)=-2.506615d0
  VIBEN(33,J)=-2.349241d0
  VIBEN(34,J)=-2.196308d0
  VIBEN(35,J)=-2.047855d0
  VIBEN(36,J)=-1.903922d0
  VIBEN(37,J)=-1.764553d0
  VIBEN(38,J)=-1.629794d0
  VIBEN(39,J)=-1.499695d0
  VIBEN(40,J)=-1.374309d0
  VIBEN(41,J)=-1.253693d0
  VIBEN(42,J)=-1.137909d0
  VIBEN(43,J)=-1.027022d0
  VIBEN(44,J)=-9.211040d-1
  VIBEN(45,J)=-8.202310d-1
  VIBEN(46,J)=-7.244866d-1
  VIBEN(47,J)=-6.339614d-1
  VIBEN(48,J)=-5.487544d-1
  VIBEN(49,J)=-4.689744d-1
  VIBEN(50,J)=-3.947413d-1
  VIBEN(51,J)=-3.261886d-1
  VIBEN(52,J)=-2.634655d-1
  VIBEN(53,J)=-2.067411d-1
  VIBEN(54,J)=-1.562088d-1
  VIBEN(55,J)=-1.120932d-1
  VIBEN(56,J)=-7.466090d-2
  VIBEN(57,J)=-4.423639d-2
  VIBEN(58,J)=-2.123195d-2
  VIBEN(59,J)=-6.208831d-3
  VIBEN(:,J)=9.821630d0+VIBEN(:,J)  !add the zero point energy
!
  J=3                   !O2 vibrational levels from: Esposito et al, Chemical Physics 351 (2008)
  IVMODEL(J,1)=1        !use pre-computed vibrational levels for this species
  IVMODEL(J,2)=46       !max possible vibrational level
  VIBEN(0,J)=-5.1153d0  !in eV
  VIBEN(1,J)=-4.9221d0
  VIBEN(2,J)=-4.7315d0
  VIBEN(3,J)=-4.5433d0
  VIBEN(4,J)=-4.3578d0
  VIBEN(5,J)=-4.1749d0
  VIBEN(6,J)=-3.9948d0
  VIBEN(7,J)=-3.8175d0
  VIBEN(8,J)=-3.6431d0
  VIBEN(9,J)=-3.4717d0
  VIBEN(10,J)=-3.3034d0
  VIBEN(11,J)=-3.1381d0
  VIBEN(12,J)=-2.9760d0
  VIBEN(13,J)=-2.8172d0
  VIBEN(14,J)=-2.6617d0
  VIBEN(15,J)=-2.5096d0
  VIBEN(16,J)=-2.3610d0
  VIBEN(17,J)=-2.2160d0
  VIBEN(18,J)=-2.0746d0
  VIBEN(19,J)=-1.9369d0
  VIBEN(20,J)=-1.8030d0
  VIBEN(21,J)=-1.6729d0
  VIBEN(22,J)=-1.5469d0
  VIBEN(23,J)=-1.4249d0
  VIBEN(24,J)=-1.3071d0
  VIBEN(25,J)=-1.1935d0
  VIBEN(26,J)=-1.0842d0
  VIBEN(27,J)=-0.97939d0
  VIBEN(28,J)=-0.87913d0
  VIBEN(29,J)=-0.78354d0
  VIBEN(30,J)=-0.69274d0
  VIBEN(31,J)=-0.60685d0
  VIBEN(32,J)=-0.52601d0
  VIBEN(33,J)=-0.45036d0
  VIBEN(34,J)=-0.38004d0
  VIBEN(35,J)=-0.31522d0
  VIBEN(36,J)=-0.25607d0
  VIBEN(37,J)=-0.20277d0
  VIBEN(38,J)=-0.15550d0
  VIBEN(39,J)=-0.11445d0
  VIBEN(40,J)=-0.079818d0
  VIBEN(41,J)=-0.051751d0
  VIBEN(42,J)=-0.030323d0
  VIBEN(43,J)=-0.015397d0
  VIBEN(44,J)=-0.0063748d0
  VIBEN(45,J)=-0.0019261d0
  VIBEN(46,J)=-0.00029275d0
  VIBEN(:,J)=5.21275d0+VIBEN(:,J)  !add the zero point energy
!
  J=5                    !NO vibrational levels from: Han's Luo (Purdue AAE) Master Thesis (2016)
  IVMODEL(J,2)=52        !max possible vibrational level
  IVMODEL(J,1)=1         !use pre-computed vibrational levels for this species
  VIBEN(0,J)=-6.495503d0 !in eV
  VIBEN(1,J)=-6.261055d0
  VIBEN(2,J)=-6.030857d0
  VIBEN(3,J)=-5.804906d0
  VIBEN(4,J)=-5.583201d0
  VIBEN(5,J)=-5.365739d0
  VIBEN(6,J)=-5.152521d0
  VIBEN(7,J)=-4.943543d0
  VIBEN(8,J)=-4.738807d0
  VIBEN(9,J)=-4.538309d0
  VIBEN(10,J)=-4.342050d0
  VIBEN(11,J)=-4.150030d0
  VIBEN(12,J)=-3.962248d0
  VIBEN(13,J)=-3.778704d0
  VIBEN(14,J)=-3.599398d0
  VIBEN(15,J)=-3.424332d0
  VIBEN(16,J)=-3.253505d0
  VIBEN(17,J)=-3.086920d0
  VIBEN(18,J)=-2.924577d0
  VIBEN(19,J)=-2.766480d0
  VIBEN(20,J)=-2.612630d0
  VIBEN(21,J)=-2.463031d0
  VIBEN(22,J)=-2.317685d0
  VIBEN(23,J)=-2.176598d0
  VIBEN(24,J)=-2.039773d0
  VIBEN(25,J)=-1.907217d0
  VIBEN(26,J)=-1.778935d0
  VIBEN(27,J)=-1.654934d0
  VIBEN(28,J)=-1.535221d0
  VIBEN(29,J)=-1.419806d0
  VIBEN(30,J)=-1.308699d0
  VIBEN(31,J)=-1.201909d0
  VIBEN(32,J)=-1.099450d0
  VIBEN(33,J)=-1.001335d0
  VIBEN(34,J)=-9.075788d-1
  VIBEN(35,J)=-8.181986d-1
  VIBEN(36,J)=-7.332134d-1
  VIBEN(37,J)=-6.526442d-1
  VIBEN(38,J)=-5.765146d-1
  VIBEN(39,J)=-5.048512d-1
  VIBEN(40,J)=-4.376836d-1
  VIBEN(41,J)=-3.750456d-1
  VIBEN(42,J)=-3.169752d-1
  VIBEN(43,J)=-2.635155d-1
  VIBEN(44,J)=-2.147164d-1
  VIBEN(45,J)=-1.706349d-1
  VIBEN(46,J)=-1.313377d-1
  VIBEN(47,J)=-9.690349d-2
  VIBEN(48,J)=-6.742670d-2
  VIBEN(49,J)=-4.302306d-2
  VIBEN(50,J)=-2.383927d-2
  VIBEN(51,J)=-1.007100d-2
  VIBEN(52,J)=-2.005209d-3
  VIBEN(:,J)=6.55879d0+VIBEN(:,J)   !add the NO dissociation energy
!
  VIBEN=VIBEN*EVOLT  !convert to Joules
END IF
!--------------------------------------------------------------
!
RETURN
END SUBROUTINE OXYGEN_NITROGEN
!
!***************************************************************************
!*************************END OF GAS DATABASE*******************************
!***************************************************************************
!
SUBROUTINE ALLOCATE_GAS
!
USE GAS
USE CALC
USE MOLECS
USE MFDSMC,only : IMF, IMFS, NMFANG, IMFpair, MFRMASS, NMFER0, NMFET0,NMFEV0,NMFER,&
  NMFET,NMFEV,NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
USE EXPCOL,only : EXPCOL_INIT
!
IMPLICIT NONE
INTEGER :: NMFpair0
!
ALLOCATE (FSP(MSP,2),SP(5,MSP),SPR(3,MSP),SPM(8,MSP,MSP),ISPR(2,MSP),ISPV(MSP),ENTR(6,MSP,2),NRSP(MSP,MSP),      &
          VMP(MSP,2),UVMP(MSP,2),VNMAX(MSP),CR(MSP),TCOL(MSP,MSP),CSCR(MSP,MSP),ISPRC(MSP,MSP),SPRC(2,MSP,MSP,MSP),STAT=ERROR)  !--isebasti: UVMP included
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
IF (MMVM > 0) THEN  !--isebasti: IDL,CPER,DPER,IREV included
  ALLOCATE (SPVM(4,MMVM,MSP),ISPVM(2,MMVM,MSP),IDL(MMVM,MSP),DPER(MSP),IVMODEL(MSP,2),VIBEN(0:100,MSP),STAT=ERROR)
  IF (ERROR /= 0) THEN
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE VIBRATION VARIABLES',ERROR
  END IF
  IVMODEL(:,1)=0   !initializing variable
  IVMODEL(:,2)=100 !initializing variable; must be the same number in the vibrational pdf sampling variables
END IF
!
IF (MNRE > 0) THEN
  ALLOCATE (CI(MNRE),AE(MNRE),AC(MNRE),BC(MNRE),ER(MNRE),REA(6,MNRE),LE(MNRE),ME(MNRE),KP(MNRE),LP(MNRE),MP(MNRE),      &
            IREA(2,MNRE),JREA(2,MNRE,2),IRCD(MNRE,MSP,MSP),NREA(2,MNRE),CPER(MNRE,2),IREV(MNRE),STAT=ERROR)
  ALLOCATE (NPVIB(2,MNRE,3,3,0:100),FPVIB(MNRE,ITMAX,3,3,0:100),FPTEMP(ITMAX),STAT=ERROR) !--isebasti: included
  ALLOCATE (NEVIB(MNRE,2,0:100),FEVIB(MNRE,ITMAX,2,0:100),STAT=ERROR) !--isebasti: included
  IF (ERROR /= 0) THEN
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE REACTION VARIABLES A',ERROR
  END IF
ELSE
  ALLOCATE (NPVIB(1,1,1,1,1),FPVIB(1,1,1,1,1),FPTEMP(1),NEVIB(1,1,1),FEVIB(1,1,1,1),STAT=ERROR) !--isebasti: included
  IF (ERROR /= 0) THEN
    WRITE (*,*)'PROGRAM COULD NOT ALLOCATE REACTION VARIABLES B',ERROR
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
!-- Allocate variables related to MF model
IF (MNRE > 0 .AND. IMF .ne. 0) THEN
  ALLOCATE(NMFANG(MNRE),IMFpair(MSP,MSP),MFRMASS(3,MNRE),STAT = ERROR)
  !--count number of pairs
  NMFpair0 = FLOOR(DBLE(MSP*(MSP+1))/2.0d0)
  IF (IMFS == 1 ) THEN
    ALLOCATE(NMFER0(1000,2,NMFpair0), NMFET0(1000,NMFpair0),NMFEV0(0:100,2,NMFpair0), &
      NMFER(1000,2,NMFpair0), NMFET(1000,NMFpair0),NMFEV(0:100,2,NMFpair0), &
      NMFERR(1000,2,MNRE),   NMFETR(1000,MNRE),  NMFEVR(0:100,2,MNRE), &
      NMFVT0(0:100,1000,2,NMFpair0), NMFVT(0:100,1000,2,NMFpair0),NMFVTR(0:100,1000,2,MNRE),STAT=ERROR)
    IF (ERROR /= 0) THEN
      WRITE(*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR MF RELATED VARIBALES', ERROR
      STOP
    ENDIF
  END IF
ENDIF

!-- Allocate variables related to nonVHS model
ALLOCATE(INONVHS(MSP,MSP), STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE(*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR nonVHS RELATED VARIBALES', ERROR
  STOP
ENDIF
INONVHS = 0
IF (GASCODE .eq. 8) THEN
  IF (nonVHS == 2) THEN
    CALL EXPCOL_INIT()
  ELSE IF (nonVHS == 1) THEN
    INONVHS(1,4) = 1
    INONVHS(4,1) = 1
  END IF
END IF

RETURN
END SUBROUTINE ALLOCATE_GAS
!
!*****************************************************************************
!
SUBROUTINE READ_RESTART
!
USE MOLECS
USE GEOM
USE GAS
USE CALC
USE OUTPUT
USE MFDSMC, only:IMF,IMFS,NMFpair,NMFANG, IMFpair, &
  MFRMASS, NMFER0, NMFET0,NMFEV0,NMFER,&
  NMFET,NMFEV,NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
!
IMPLICIT NONE
!
INTEGER :: ZCHECK,ZCHECK2
ZCHECK = 0; ZCHECK2 = 0
!
101 CONTINUE
OPEN (7,FILE='PARAMETERS.DAT',FORM='UNFORMATTED',ERR=101) !--isebasti: replace binary by unformatted
READ (7) NCCELLS,NCELLS,MMRM,MMVM,MNM,MNRE,MNSR,MSP,MTBP,ILEVEL,MDIV,IRECOM,MMEX,MEX,ISF,NBINS,NSNAP,ITMAX,IMF, IMFS, NMFpair,nonVHS !--isebasti: included ISF,NBINS,NSAP
CLOSE(7)
!
IF (MMVM > 0) THEN
  ALLOCATE (PX(MNM),PTIM(MNM),PROT(MNM),IPCELL(MNM),IPSP(MNM),ICREF(MNM),IPCP(MNM),PV(3,MNM),      &
       IPVIB(0:MMVM,MNM),STAT=ERROR)
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
          VARSP(0:11,NCELLS,MSP),VARS(0:32+MSP,2),CS(0:8+MMVM,NCELLS,MSP),CSS(0:8,2,MSP,2), &  !--isebasti: correcting CS allocation
          CSSS(6,2),CST(0:4,NCELLS),BINS(0:NBINS,5,MSP),BIN(0:NBINS,5),&
          PDFS(0:NBINS,5,MSP),PDF(0:NBINS,5),NDROT(MSP,100),NDVIB(NSCELLS,0:MMVM,MSP,0:100),STAT=ERROR) !--isebasti: CST,BINS,BIN,PDFS,PDF included
IF (MNRE > 0) THEN
  ALLOCATE (NPVIB(2,MNRE,3,3,0:100),FPVIB(MNRE,ITMAX,3,3,0:100),FPTEMP(ITMAX),STAT=ERROR)
  ALLOCATE (NEVIB(MNRE,2,0:100),FEVIB(MNRE,ITMAX,2,0:100),STAT=ERROR) !--isebasti: included
ELSE
  ALLOCATE (NPVIB(1,1,1,1,1),FPVIB(1,1,1,1,1),FPTEMP(1),NEVIB(1,1,1),FEVIB(1,1,1,1),STAT=ERROR) !--isebasti: included
  ALLOCATE (EVREM(MNRE),STAT=ERROR)
END IF
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR SAMPLING ARRAYS',ERROR
ENDIF
!
!
CALL ALLOCATE_GAS
!
102 CONTINUE
OPEN (7,FILE='RESTART.DAT',FORM='UNFORMATTED',ERR=102) !--isebasti: replace binary by unformatted
READ (7) AC,AE,AVDTM,BC,BOLTZ,EVOLT,CCELL,CELL,CI,CLSEP,COLLS, &
         CPDTM,CR,CS,CSS,CSSS,CTM,CXSS,CST,BINS,BIN,PDFS,PDF,NSPDF,NDROT,NDVIB,DDIV,DPI,DTM,DTSAMP,DTOUT, &
         ENTMASS,ENTR,ER,ERROR,ERS,FDEN,FMA,FND,FNUM,FRACSAM,FSP,FP,FPR,FREM,FSPEC, &
         FTMP,FTIME,FRTMP,FVTMP,GASCODE,ICCELL,ICELL,ICREF,IFX,IMTS,IPCELL,IPCP, &
         IPSP,IPVIB,IRCD,IREA,IREM,ISECS,ISF,ISPEX,ISPR,ISPRC,ISPV,ISPVM,ISRCD,ITYPE,IVB,IWF,JCD, &
         JDIV,JREA,KP,LE,LIS,LP,LRS,ME,MOLSC,MP,MVER,NCCELLS,NCELLS, &
         NCIS,NDIV,NEX,NM,NMISAMP,NNC,NOUT,NREA,NRSP,NSAMP,NSPEX,NVER,OUTRAT,SAMPRAT,PI,PROT,PTIM,PV,PX, &
         REAC,RGFS,RMAS,REA,SP,SPEX,SPI,SPM,SPR,SPRC,SPV,SPVM,IDL,CPER,DPER,ENERS,ENERS0,IREV,SREAC,NPVIB,FPVIB,FPTEMP,IRM, &
         TCOL,CSCR,TDISS,TFOREX,TRECOMB,TREVEX,THBP,TISAMP,TPOUT,TREF,TLIM,TOTCOL,PCOLLS,TOTDUP,TOTMOV,     &
         TREACG,TREACL,TOUT,TPDTM,TREF,TSAMP,TSURF,VAR,VARS,VARSP,VELOB,VFX,VFY,UVFX,UFND,UFTMP,UVMP,VMP, &
         VMPM,VNMAX,VSURF,WCOLLS,WFM,XB,XREM,XVELS,YVELS,TNEX,NEVIB,FEVIB,QCTMODEL,IVMODEL,VIBEN,IGS,AMEG,EVREM,ZCHECK
!--isebasti: CST,BINS,BIN,PDFS,PDF,UVFX,UFND,UFTMP,UVMP,IDL,CPER,DPER,ENERS,ENERS0,IREV,NPVIB,FPVIB,FPTEMP included
IF (ZCHECK /= 1234567) THEN
  WRITE (9,*) 'Wrong', NM,' Molecules, Check integer =',ZCHECK
  STOP
END IF

IF (MNRE >0 .AND. IMF .ne. 0) THEN
  READ(7) NMFANG, IMFpair,MFRMASS
  IF (IMFS == 1 .AND. NMFpair > 0)THEN
    READ(7) NMFER0,NMFET0,NMFEV0,NMFER,NMFET,NMFEV,&
      NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
  END IF
  READ(7) ZCHECK2
  IF (ZCHECK2 /= 1234567) THEN
    WRITE (9,*) 'WrongIMFpair, Check integer =',ZCHECK2
    STOP
  END IF
END IF
CLOSE(7)
!
RETURN
END SUBROUTINE READ_RESTART
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
USE MFDSMC, only:IMF,IMFS,NMFpair,NMFANG, IMFpair, &
  MFRMASS, NMFER0, NMFET0,NMFEV0,NMFER,&
  NMFET,NMFEV,NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
!
IMPLICIT NONE
!
INTEGER :: ZCHECK,ZCHECK2
!
ZCHECK=1234567
ZCHECK2 = 1234567
!
101 CONTINUE
OPEN (7,FILE='PARAMETERS.DAT',FORM='UNFORMATTED',ERR=101) !--isebasti: replace binary by unformatted
WRITE (7) NCCELLS,NCELLS,MMRM,MMVM,MNM,MNRE,MNSR,MSP,MTBP,ILEVEL,MDIV,IRECOM,MMEX,MEX,ISF,NBINS,NSNAP,ITMAX,IMF,IMFS,NMFpair,nonVHS !--isebasti: included ISF,NBINS,NSAP,ITMAX
CLOSE(7)
!
102 CONTINUE
OPEN (7,FILE='RESTART.DAT',FORM='UNFORMATTED',ERR=102) !--isebasti: replace binary by unformatted
WRITE (7)AC,AE,AVDTM,BC,BOLTZ,EVOLT,CCELL,CELL,CI,CLSEP,COLLS, &
         CPDTM,CR,CS,CSS,CSSS,CTM,CXSS,CST,BINS,BIN,PDFS,PDF,NSPDF,NDROT,NDVIB,DDIV,DPI,DTM,DTSAMP,DTOUT, &
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
IF (MNRE >0 .AND. IMF .ne. 0) THEN
  WRITE(7) NMFANG, IMFpair,MFRMASS
  IF (IMFS == 1 .AND. NMFpair > 0)THEN
    WRITE(7) NMFER0,NMFET0,NMFEV0,NMFER,NMFET,NMFEV,&
      NMFERR,NMFETR,NMFEVR,NMFVT0,NMFVT,NMFVTR
  END IF
  WRITE(7) ZCHECK2
END IF
CLOSE(7)
!
WRITE (9,*) 'Restart files written'
!
RETURN
END SUBROUTINE WRITE_RESTART
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
INTEGER :: N,L,M,NSC,NCC,ND
REAL(KIND=8) :: X,FRAC,DSC,A,B,C,TIM
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
!
!******************************************************************************
!
FUNCTION GAM(X)
!
!--calculates the Gamma function of X.
!
IMPLICIT NONE
REAL(8) ::X,A,Y,GAM
A=1.D0
Y=X
IF (Y < 1.D0) THEN
  A=A/Y
ELSE
  Y=Y-1.D0
  DO WHILE (Y >= 1.D0)
    A=A*Y
    Y=Y-1.D0
  END DO
END IF
GAM=A*(1.D0-0.5748646D0*Y+0.9512363D0*Y**2-0.6998588D0*Y**3+  &
       0.4245549D0*Y**4-0.1010678D0*Y**5)
!
RETURN
!
END FUNCTION GAM
!
!*****************************************************************************
!
FUNCTION ERF(S)
!
!--evaluates the error function of S
!
IMPLICIT NONE
REAL(8) :: S,B,C,T,D,ERF
B=DABS(S)
IF (B < 4.D0) THEN
  C=DEXP(-B*B)
  T=1.D0/(1.D0+0.3275911D0*B)
  D=1.D0-(0.254829592D0*T-0.284496736D0*T*T+1.421413741D0*T*T*T- &
    1.453152027D0*T*T*T*T+1.061405429D0*T*T*T*T*T)*C
ELSE
  D=1.D0
END IF
IF (S < 0.D0) D=-D
ERF=D
RETURN
END FUNCTION ERF
!
!*****************************************************************************
!
SUBROUTINE RVELC(U,V,VMP,IDT)
!
USE CALC
!
IMPLICIT NONE
!
!--generates two random velocity components U and V in an equilibrium
!--gas with most probable speed VMP
INTEGER :: IDT !included IDT
REAL(KIND=8) :: U,V,VMP,A,B,RANF !--isebasti: included RANF
!
CALL ZGF(RANF,IDT)
A=DSQRT(-DLOG(RANF))
CALL ZGF(RANF,IDT)
B=DPI*RANF
U=A*DSIN(B)*VMP
V=A*DCOS(B)*VMP
RETURN
!
END SUBROUTINE RVELC
!
!*****************************************************************************
!
SUBROUTINE SROT(L,TEMP,ROTE,IDT)
!
!--sets a typical rotational energy ROTE of species L
!
USE CALC
USE GAS
!
IMPLICIT NONE
!
INTEGER :: I,L,IDT !included IDT
REAL(KIND=8) :: A,B,ROTE,ERM,TEMP,RANF !--isebasti: included RANF
!
IF (ISPR(1,L).EQ.2) THEN
  CALL ZGF(RANF,IDT)
  ROTE=-DLOG(RANF)*BOLTZ*TEMP
ELSE
  A=0.5D00*ISPR(1,L)-1.D00
  I=0
  DO WHILE (I == 0)
    CALL ZGF(RANF,IDT)
    ERM=RANF*10.D00
!--there is an energy cut-off at 10 kT
    B=((ERM/A)**A)*DEXP(A-ERM)
    CALL ZGF(RANF,IDT)
    IF (B > RANF) I=1
  END DO
  ROTE=ERM*BOLTZ*TEMP
END IF
!
RETURN
!
END SUBROUTINE SROT
!
!*****************************************************************************
!
SUBROUTINE SVIB(L,TEMP,IVIB,K,IDT)
!
!--sets a typical vibrational state at temp. T of mode K of species L
!
USE GAS
USE CALC
USE MOLECS
!
IMPLICIT NONE
!
INTEGER :: II,K,L,IDT,IVIB,MAXLEVEL !included IDT
REAL(KIND=8) :: EVIB,TEMP,RANF,QVIB,PROB  !--isebasti: included RANF
!
IF (TEMP < 5.d0) THEN
  IVIB=0
  RETURN
END IF
!
!old approach (limited for SHO)
!CALL ZGF(RANF,IDT)
!IVIB=-DLOG(RANF)*TEMP/SPVM(1,K,L) !eqn(11.24)
!
!calculate the vibrational partition function
QVIB=0.d0
MAXLEVEL=IVMODEL(L,2)
DO IVIB=0,MAXLEVEL
  CALL VIB_ENERGY(EVIB,IVIB,K,L)
  QVIB=QVIB+DEXP(-EVIB/(BOLTZ*TEMP))
END DO
!
!sample a level from Boltzmann distribution
II=0
DO WHILE (II == 0)
  CALL ZGF(RANF,IDT)
  IVIB=RANF*(MAXLEVEL+0.99999999D00)
  CALL VIB_ENERGY(EVIB,IVIB,K,L)
  PROB=DEXP(-EVIB/(BOLTZ*TEMP))/QVIB
  CALL ZGF(RANF,IDT)
  IF (PROB > RANF) II=1
END DO
!
RETURN
END SUBROUTINE SVIB
!
!*****************************************************************************
!
SUBROUTINE INIT_REAC
!
!--initialise the reaction variables (in the old reaction model)
!
USE MOLECS
!
USE GAS
USE CALC
USE GEOM
USE OUTPUT
USE MFDSMC, only : IMF,IMFS,MFRMASS, NMFANG,NMFpair, IMFpair,MF_SET_AHO, IMFdia
!
IMPLICIT NONE
!
INTEGER :: I,K,KK,L,M,N,LS,MS,NNRE,NTBP,II,JJ
REAL(KIND=8) :: A,B,C,D,EPS,AL,X,AA,BB
REAL(8),EXTERNAL :: GAM
!
!--I,K,L,M,N working integers
!--A,B working variables
!--LS,MS species
!--EPS symmetry factor
!
!
!--convert the reaction data into the required form
!
IF (MNRE == 0) THEN
  NNRE=1
ELSE
  NNRE=MNRE
END IF
IF (MTBP == 0) THEN
  NTBP=1
ELSE
  NTBP=MTBP
END IF
!
IF (MNRE > 0) THEN
  REA=0.; IREA=0; NREA=0; JREA=0; NRSP=0; IRCD=0 ; REAC=0.; EVREM = 0.d0
  IF ( IMF .ne. 0) THEN
     MFRMASS=-100.d0; NMFANG = 0;
  END IF
!
  OPEN(10, FILE="ChemicalReaction.txt")
  WRITE(10, "(A, I3)") "Number of reactions",MNRE
  DO  N=1,MNRE
    WRITE(10, "(A,I3)") "Reaction ",N
    L=LE(N)
    M=ME(N)
!--the pre-reaction species codes are L and M
    IF ((L > 0).AND.(M > 0)) THEN
      IREA(1,N)=L
      IREA(2,N)=M
      IF (M < L) THEN
        K = L
        L = M
        M = K
      END IF
      NRSP(L,M) = NRSP(L,M) + 1
      IF (L /= M) NRSP(M,L) = NRSP(M,L) + 1
      K = NRSP(L,M)
!--this is the Kth reaction involving species L and M
      IRCD(K,M,L)=N
      IRCD(K,L,M)=N
    END IF

!============= set post reaction info
    K=KP(N)
    L=LP(N)
    M=MP(N)
!--three post-collision particles (dissociation)
    IF (K > 0) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Dissociation'
      WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3,' + ',I3)") LE(N),ME(N),K,L,M
      WRITE(10,*)
      NREA(1,N)=2
      NREA(2,N)=1
      JREA(1,N,1)=K
      JREA(1,N,2)=L
      JREA(2,N,1)=M
    END IF
!--two post-collision particles (exchange)
    IF (K == 0) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Exchange'
      WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3)") LE(N),ME(N),L,M
      WRITE(10,*)
      NREA(1,N)=1
      NREA(2,N)=1
      JREA(1,N,1)=L
      JREA(2,N,1)=M
    END IF
!--one post collision particle (recombination)
    IF (K == -1) THEN
      WRITE(10, "(A,A)") 'ReactionType: ', 'Recombination'
      WRITE(10, "(I3,' + ',I3,' + ',I3,' -> ',I3,' + ',I3)") LE(N),ME(N),M,L,M
      WRITE(10,*)
      NREA(1,N)=1
      NREA(2,N)=0
      JREA(1,N,1)=L
      JREA(2,N,1)=M
    END IF
!
    REA(1,N)=CI(N)
    REA(2,N)=AE(N)
    REA(3,N)=AC(N)
    REA(4,N)=BC(N)
    REA(5,N)=ER(N)
  END DO
  CLOSE(10)
!
END IF
!
DO L=1,MNRE
  DO M=1,MSP
    DO N=1,MSP
      IF (IRCD(L,N,M) == 0) IRCD(L,N,M)=IRCD(L,M,N)
    END DO
  END DO
END DO
!
!--calculate the coefficients of the energy terms in steric factors (see SMILE documentation)
!
NMFpair = 0
IF (IMF .ne. 0 .and. IMFS == 1) THEN
  DO I=1,MSP
    DO K = I,MSP
      NMFpair = NMFpair+1
      IMFpair(I,K) = NMFpair
      IMFpair(K,I) = NMFpair
    END DO
  END DO
END IF


IF (IMF == 2 .or. IMF == 3) THEN
  CALL MF_SET_AHO()
END IF

IF (IMF .ne. 0) THEN
  OPEN(10, FILE="ChemicalReaction.txt", ACCESS="APPEND")
  WRITE(10,*)
  WRITE(10,*)
  IF (IMF == 1) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-SHO  Dia:',IMFdia
  ELSE IF (IMF == 2) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-AHO QCT vphase  Dia:',IMFdia
  ELSE IF (IMF == 3) THEN
    WRITE(10,'(A,1X,I2)') 'Macheret-Fridman model: MF-AHO Morse vphase  Dia:',IMFdia
  END IF
  WRITE(10,'(A)') "Macheret-Fridman model parameter:"
END IF

DO K=1,MNRE
  LS=IREA(1,K)
  MS=IREA(2,K)
  EPS=1.D00
  IF (LS == MS) EPS=2.D00
  I=NREA(1,K)+NREA(2,K)
  AL=SPM(3,LS,MS)-0.5d0  !alpha
  X=REA(4,K)-0.5d0+AL
!
  IF (I == 1) THEN
!--recombination reaction
    C=EPS*SPI*REA(3,K)*SPM(1,LS,MS)**0.5d0
    D=(2.d0**1.5d0)*(BOLTZ**REA(4,K))*SPM(2,LS,MS)*((2.d0-AL)*BOLTZ*SPM(5,LS,MS))**AL
    REA(6,K)=(C/D)
  END IF
  IF (I == 2) THEN
!--exchange reaction
    BB=2.d0-AL
    B=GAM(BB)
    C=EPS*SPI*REA(3,K)*(SPM(1,LS,MS)/(8.d0*BOLTZ))**0.5d0
!   D=SPM(2,LS,MS)*(BOLTZ**X)*(2.d0-AL*SPM(5,LS,MS))**AL  !Alexeenko-SMILE report with bug
!   REA(6,K)=(C/D)*(1.d0/B)                               !Alexeenko-SMILE report with bug
    D=SPM(2,LS,MS)*(BOLTZ**X)*(SPM(5,LS,MS))**AL          !Ivanov-SMILE report
    REA(6,K)=(C/D)                                        !Ivanov-report
  END IF
!--dissociation reaction
  IF (I == 3) THEN
    IF (REA(1,K) < 0.) THEN !TCE
      BB=2.d0-AL
      B=GAM(BB)
      C=EPS*SPI*REA(3,K)*(SPM(1,LS,MS)/(8.d0*BOLTZ))**0.5d0
      D=SPM(2,LS,MS)*(BOLTZ**X)*(SPM(5,LS,MS))**AL        !Ivanov-SMILE report
      REA(6,K)=(C/D)                                      !Ivanov-SMILE report
    ELSE  !VDC
      BB=2.d0-AL
      B=GAM(BB)
      C=EPS*REA(3,K)*SPI*0.5d0*(SPM(1,LS,MS)/(2.d0*BOLTZ))**(0.5d0-AL)
      D=SPM(2,LS,MS)*B
      REA(6,K)=(C/D)
    END IF

    IF (IMF .ne. 0) THEN  !pre set for Macheret Fridman model
      IF (GASCODE .ne. 8)THEN
        STOP "Only gascode = 8 works with MF model"
      END IF
      ! LS is always the one to be dissociated
      ! MFRAMSS(1): reduced mass of whole system
      IF ((LS .gt. 5) .or. (MS .gt. 5)) THEN
        WRITE(10, "(A, I3, A, I3,A)") "Reaction: ",LS, ' + ',MS, ' is not supported by MF model'
        WRITE(10, "(I3,' + ',I3,' -> ',I3,' + ',I3,' + ',I3)") LE(K),ME(K),KP(K),LP(K),MP(K)
        WRITE(10,*)
        MFRMASS(:,K) = 0.0d0
        NMFANG(K) = 0
      ELSE
      !----- calculated reduced mass
      ! reduced mass of A+B
        MFRMASS(1,K) = SP(5,LS)*SP(5,MS)/(SP(5,LS)+SP(5,MS))
        DO II = 1,2
          IF (ISPV(IREA(II,K)) == 0) THEN
          ! single atom
            MFRMASS(II+1,K) = SP(5,IREA(II,K))
          ELSE IF (IREA(II,K) == 5) THEN
          ! NO molecule
            MFRMASS(II+1,K) = SP(5,2)*SP(5,4)/(SP(5,2)+SP(5,4))
          ELSE IF (IREA(II,K) == 3 .or. IREA(II,K) == 1) THEN
          ! nitogen, oxygen and other homonuclear molecule
          ! reduced mass should be 1/4 of the mass of molecule
            MFRMASS(II+1,K) = SP(5,IREA(II,K))/4.0d0
          ELSE
            WRITE(*,*) "Not supported MFRMASS", IREA(II,K)
            STOP
          END IF
        END DO
        WRITE(10, "(A, I3, A, I3)") "Reaction: ",LS, ' + ',MS
        WRITE(10, "(A, 3(G10.3, 1X))") "Reduced mass (kg): ", MFRMASS(1,K), MFRMASS(2,K),MFRMASS(3,K)
      !
      !----- determine reaction type
        NMFANG(K) = 4 ! default: homo-atom
        IF (ISPV(MS) == 1) NMFANG(K) = 8  !default: homo-homo
        IF (LS == 5) NMFANG(K) = NMFANG(K)+1 !hetero-homo or hetero-atom: 5,9
        IF (MS == 5) NMFANG(K) = NMFANG(K)+2
      ! ----- homo   hetero   !dissociated molecule
      !atom    4       5
      !homo    8       9
      !hetero  10     11
      !
      END IF
    END IF
  END IF
!
END DO
IF (IMF .ne. 0) THEN
  CLOSE(10)
END IF
!
RETURN
!
END SUBROUTINE INIT_REAC
!
!*****************************************************************************
!
SUBROUTINE MOLECULES_ENTER
!
!--molecules enter boundary at XB(1) and XB(2) and may be removed behind a wave
!
USE MOLECS
USE GAS
USE CALC
USE GEOM
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: K,L,M,N,NENT,II,J,JJ,KK,NTRY,IDT=0 !included IDT
REAL(KIND=8) :: A,B,AA,BB,U,VN,XI,X,DX,DY,DZ,RANF !--isebasti: included RANF
!
!--NENT number to enter in the time step
!
DO J=1,2       !--J is the end
  IF (ITYPE(J) == 0) THEN
    KK=1 !--the entry surface will normally use the reference gas (main stream) properties
    IF ((J == 2).OR.((ISECS == 1).AND.(XB(2) > 0.D00))) KK=2  !--KK is 1 for reference gas 2 for the secondary stream  !--isebasti: AND replaced by OR
    DO L=1,MSP
      A=ENTR(1,L,J)*DTM+ENTR(2,L,J)
      NENT=A
      ENTR(2,L,J)=A-NENT
      IF (NENT > 0) THEN
        DO M=1,NENT
          IF (NM >= MNM) THEN
            WRITE (*,*) 'EXTEND_MNM from MOLECULES_ENTER'
            CALL EXTEND_MNM(1.1d0)
          END IF
          NM=NM+1
          AA=DMAX1(0.D00,ENTR(3,L,J)-3.D00)
          BB=DMAX1(3.D00,ENTR(3,L,J)+3.D00)
          II=0
          DO WHILE (II == 0)
            CALL ZGF(RANF,IDT)
            B=AA+(BB-AA)*RANF
            U=B-ENTR(3,L,J)
            A=(2.D00*B/ENTR(4,L,J))*DEXP(ENTR(5,L,J)-U*U)
            CALL ZGF(RANF,IDT)
            IF (A > RANF) II=1
          END DO
          PV(1,NM)=B*UVMP(L,KK)  !--isebasti: VMP replaced by UVMP
          IF (J == 2) PV(1,NM)=-PV(1,NM)
!
          CALL RVELC(PV(2,NM),PV(3,NM),UVMP(L,KK),IDT)  !--isebasti: VMP replaced by UVMP
          PV(2,NM)=PV(2,NM)+VFY(J)
!
          IF (ISPR(1,L) > 0) CALL SROT(L,UFTMP(KK),PROT(NM),IDT)  !--isebasti: UFTMP replaced by UFTMP
!
          IF (MMVM > 0) THEN
            IPVIB(0,NM)=0
            DO K=1,ISPV(L)
              CALL SVIB(L,UFTMP(KK),IPVIB(K,NM),K,IDT)  !--isebasti: FTMP replaced by UFTMP
            END DO
          END IF
          IPSP(NM)=L
!--advance the molecule into the flow
          CALL ZGF(RANF,IDT)
          XI=XB(J)
          DX=DTM*RANF*PV(1,NM)
          IF (IFX == 0) X=XI+DX
          IF (IFX > 0) DY=DTM*RANF*PV(2,NM)
          DZ=0.D00
          IF (IFX == 2) DZ=DTM*RANF*PV(3,NM)
          IF (IFX > 0) CALL AIFX(XI,DX,DY,DZ,X,PV(1,NM),PV(2,NM),PV(3,NM),IDT)
          IF (IFX == 0) PX(NM)=X
          PTIM(NM)=FTIME
          IF (IVB == 0) CALL FIND_CELL(PX(NM),IPCELL(NM),JJ)
          IF (IVB == 1) CALL FIND_CELL_MB(PX(NM),IPCELL(NM),JJ,PTIM(NM))
          IPCP(NM)=0
          IF (XREM > XB(1)) ENTMASS=ENTMASS+SP(5,L)
        END DO
      END IF
    END DO
  END IF
END DO
!
!--stagnation streamline molecule removal
IF (XREM > XB(1)) THEN
  ENTMASS=FREM*ENTMASS
  NTRY=0
  DO WHILE ((ENTMASS > 0.D00).AND.(NTRY < 10000))
    NTRY=NTRY+1
    IF (NTRY == 10000) THEN
      WRITE (*,*) 'Unable to find molecule for removal'
      ENTMASS=0.D00
      VNMAX=0.D00
    END IF
    CALL ZGF(RANF,IDT)
    N=NM*RANF+0.9999999D00
    IF (PX(N) > XREM) THEN
      CALL ZGF(RANF,IDT)
      !IF (RANF < ((PX(N)-XREM)/(XB(2)-XREM))**2) THEN
      IF (DABS(VFY(1)) < 1.D-3) THEN
        VN=DSQRT(PV(2,N)*PV(2,N)+PV(3,N)*PV(3,N))   !--AXIALLY SYMMETRIC STREAMLINE
      ELSE
        VN=DABS(PV(3,N))   !--TWO-DIMENSIONAL STREAMLINE
      END IF
      L=IPSP(N)
      IF (VN > VNMAX(L)) VNMAX(L)=VN
      CALL ZGF(RANF,IDT)
      IF (RANF < VN/VNMAX(L)) THEN
        CALL REMOVE_MOL(N)
        ENTMASS=ENTMASS-SP(5,L)
        NTRY=0
      END IF
      !END IF
    END IF
  END DO
END IF

!
END SUBROUTINE MOLECULES_ENTER
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
!
!*****************************************************************************
!
SUBROUTINE MOLECULES_MOVE
!
!--molecule moves appropriate to the time step
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: N,L,M,K,NCI,J,II,JJ,IDT=0 !included IDT,iflag
REAL(KIND=8) :: A,B,X,XI,XC,DX,DY,DZ,DTIM,S1,XM,R,TI,DTC,POB,UR,WFI,WFR,WFRI,RANF !--isebasti: included RANF
!
!--N working integer
!--NCI initial collision cell
!--DTIM time interval for the move
!--POB position of the outer boundary
!--TI initial time
!--DTC time interval to collision with surface
!--UR radial velocity component
!--WFI initial weighting factor
!--WFR weighting factor radius
!--WFRI initial weighting factor radius
!
!$omp parallel &
!$omp private(idt,n,nci,dtim,wfi,ii,ti,xi,dx,x,dz,dy,r,j,l,dtc,xc,s1,wfr,wfri,k,m,jj) &
!$omp reduction(+:totmov,entmass)
!$    idt=omp_get_thread_num()  !thread id
!$omp do !schedule (dynamic,NCCELLS/32)
!
DO N=1,NM
!
  NCI=IPCELL(N)
  IF (IMTS == 0) DTIM=DTM
  IF (IMTS == 1) DTIM=2.D00*CCELL(3,NCI)
  IF (FTIME-PTIM(N) > 0.5*DTIM) THEN
    WFI=1.D00
    IF (IWF == 1) WFI=1.D00+WFM*PX(N)**IFX
    II=0 !--becomes 1 if a molecule is removed
    TI=PTIM(N)
    PTIM(N)=TI+DTIM
    TOTMOV=TOTMOV+1
!
    XI=PX(N)
    DX=DTIM*PV(1,N)
    X=XI+DX
!
    IF (IFX > 0) THEN
      DZ=0.D00
      DY=DTIM*PV(2,N)
      IF (IFX == 2) DZ=DTIM*PV(3,N)
      R=DSQRT(X*X+DY*DY+DZ*DZ)
    END IF
!
    IF (IFX == 0) THEN
      DO J=1,2    ! 1 for minimum x boundary, 2 for maximum x boundary
        IF (II == 0) THEN
          IF (((J == 1).AND.(X < XB(1))).OR.((J == 2).AND.(X > (XB(2)+VELOB*PTIM(N))))) THEN  !--molecule crosses a boundary
            IF ((ITYPE(J) == 0).OR.(ITYPE(J) == 3)) THEN
              IF (XREM > XB(1)) THEN
                L=IPSP(N)
                ENTMASS=ENTMASS-SP(5,L)
              END IF
              IPCELL(N)=-IPCELL(N) !molecule is marked for removel !--isebasti: original use CALL REMOVE_MOL(N); !N=N-1
              II=1
            END IF
!
            IF (ITYPE(J) == 1) THEN
              IF ((IVB == 0).OR.(J == 1)) THEN
                X=2.D00*XB(J)-X
                PV(1,N)=-PV(1,N)
              ELSE IF ((J == 2).AND.(IVB == 1)) THEN
                DTC=(XB(2)+TI*VELOB-XI)/(PV(1,N)-VELOB)
                XC=XI+PV(1,N)*DTC
                PV(1,N)=-PV(1,N)+2.*VELOB
                X=XC+PV(1,N)*(DTIM-DTC)
              END IF
            END IF
!
            IF (ITYPE(J) == 2) THEN
              CALL REFLECT(N,J,X,IDT)
                IF((J == 2).AND.(X < XB(1))) THEN !--isebasti: included to fix bug
                  DO WHILE ((X<XB(1)).OR.(X>XB(2)))
                   !WRITE(9,*) 'REFLECTED MOLEC CROSSING OPPOSITE BOUNDARY:',FTIME,N,X,PV(1,N)
                    IF (X < XB(1)) CALL REFLECT(N,1,X,IDT) !molecules reflected at xb2 crossed xb1
                    IF (X > XB(2)) CALL REFLECT(N,2,X,IDT) !molecules reflected at xb1 crossed xb2
                   !WRITE(9,*) 'REFLECTED MOLEC CROSSING OPPOSITE BOUNDARY:',FTIME,N,X,PV(1,N)
                  END DO
                END IF
            END IF
          END IF
        END IF
      END DO
    ELSE         !--cylindrical or spherical flow
!--check boundaries
      IF ((X < XB(1)).AND.(XB(1) > 0.D00)) THEN
        CALL RBC(XI,DX,DY,DZ,XB(1),S1)
        IF (S1 < 1.D00) THEN     !--intersection with inner boundary
          IF (ITYPE(1) == 2) THEN !--solid surface
            DX=S1*DX
            DY=S1*DY
            DZ=S1*DZ
            CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
            CALL REFLECT(N,1,X,IDT)
          ELSE
            IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
            II=1
          END IF
        END IF
      ELSE IF ((IVB == 0).AND.(R > XB(2))) THEN
        CALL RBC(XI,DX,DY,DZ,XB(2),S1)
        IF (S1 < 1.D00) THEN     !--intersection with outer boundary
          IF (ITYPE(2) == 2) THEN !--solid surface
            DX=S1*DX
            DY=S1*DY
            DZ=S1*DZ
            CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
            X=1.001D00*XB(2)
            DO WHILE (X > XB(2))
              CALL REFLECT(N,2,X,IDT)
            END DO
          ELSE
            IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
            II=1
          END IF
        END IF
      ELSE IF ((IVB == 1).AND.(R > (XB(2)+PTIM(N)*VELOB))) THEN
        IF (IFX == 1) UR=DSQRT(PV(1,N)**2+PV(2,N)**2)
        IF (IFX == 2) UR=DSQRT(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        DTC=(XB(2)+TI*VELOB-XI)/(UR-VELOB)
        S1=DTC/DTIM
        DX=S1*DX
        DY=S1*DY
        DZ=S1*DZ
        CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
        PV(1,N)=-PV(1,N)+2.*VELOB
        X=X+PV(1,N)*(DTIM-DTC)
      ELSE
        CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
      END IF
!
!--DIAGNOSTIC
      IF (X > XB(2)+PTIM(N)*VELOB) THEN
        WRITE (*,*) N,FTIME,X,XB(2)+PTIM(N)*VELOB
      END IF
!
!--Take action on weighting factors
      IF ((IWF == 1).AND.(II >= 0)) THEN
        WFR=WFI/(1.D00+WFM*X**IFX)
        L=0
        WFRI=WFR
        IF (WFR >= 1.D00) THEN
          DO WHILE (WFR >= 1.D00)
            L=L+1
            WFR=WFR-1.D00
          END DO
        END IF
        CALL ZGF(RANF,IDT)
        IF (RANF <= WFR) L=L+1
        IF (L == 0) THEN
          IPCELL(N)=-IPCELL(N) !--isebasti: CALL REMOVE_MOL(N); !N=N-1
          II=1
        END IF
        L=L-1
        IF (L > 0) THEN
          DO K=1,L
!$omp critical
            IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
            NM=NM+1  !--isebasti: for spherical and cylindrical geometries, it will cause a problem because NM value changes.
            PX(NM)=X
            DO M=1,3
              PV(M,NM)=PV(M,N)
            END DO
            IF (MMRM > 0) PROT(NM)=PROT(N)
            IPCELL(NM)=ABS(IPCELL(N))  !--isebasti: I think we should remove ABS
            IPSP(NM)=IPSP(N)
            IPCP(NM)=IPCP(N)
            IF (MMVM > 0) THEN
              DO M=1,MMVM
                IPVIB(M,NM)=IPVIB(M,N)
              END DO
            END IF
            PTIM(NM)=PTIM(N)    !+5.D00*DFLOAT(K)*DTM
!--note the possibility of a variable time advance that may take the place of the duplication buffer in earlier programs
!
            IF (PX(NM) > XB(2)+PTIM(NM)*VELOB) THEN
              WRITE (*,*) 'DUP',NM,FTIME,PX(NM),XB(2)+PTIM(NM)*VELOB
            END IF
!$omp end critical
          END DO
        END IF
      END IF
    END IF
!
    IF (II == 0) THEN
      PX(N)=X
        if ((px(n) > xb(1)).and.(px(n) < xb(2))) then
          continue
        else
!$omp critical(move)
          write (*,*) n,'Outside flowfield at',px(n)
!$omp end critical(move)
        end if
      IF (IVB == 0) CALL FIND_CELL(PX(N),IPCELL(N),JJ)
      IF (IVB == 1) CALL FIND_CELL_MB(PX(N),IPCELL(N),JJ,PTIM(N))
    END IF
!
  END IF
!
END DO
!$omp end do
!$omp end parallel
!
!--isebasti: remove marked molecules
N=0
DO WHILE (N < NM)
  N=N+1
  IF (IPCELL(N) < 0) THEN
    CALL REMOVE_MOL(N)
    N=N-1
  END IF
END DO
!
RETURN
!
END SUBROUTINE MOLECULES_MOVE
!
!*****************************************************************************
!
SUBROUTINE REFLECT(N,J,X,IDT)
!
!--reflects molecule N and samples the surface J properties
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: N,J,L,K,M,IDT !included IDT
REAL(KIND=8) :: A,B,VMPS,DTR,X,XI,DX,DY,DZ,WF,EVIB,RANF !--isebasti: included RANF
!
!--VMPS most probable velocity at the surface temperature
!--DTR time remaining after molecule hits a surface
!
!iflag=0; if ((n == 66430).and.(ftime > 1.5255121136)) iflag=1 !--isebasti: lines for debugging
!
L=IPSP(N)
WF=1.D00
IF (IWF == 1) WF=1.D00+WFM*X**IFX
!$omp critical(reflect1)
CSS(0,J,L,1)=CSS(0,J,L,1)+1.D00
CSS(1,J,L,1)=CSS(1,J,L,1)+WF
CSS(2,J,L,1)=CSS(2,J,L,1)+WF*PV(1,N)*SP(5,L)
CSS(3,J,L,1)=CSS(3,J,L,1)+WF*(PV(2,N)-VSURF(J))*SP(5,L)
CSS(4,J,L,1)=CSS(4,J,L,1)+WF*PV(3,N)*SP(5,L)
A=(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
CSS(5,J,L,1)=CSS(5,J,L,1)+WF*0.5D00*SP(5,L)*A
IF (ISPR(1,L) > 0) CSS(6,J,L,1)=CSS(6,J,L,1)+WF*PROT(N)
IF (MMVM > 0) THEN
  IF (ISPV(L) > 0) THEN
    DO K=1,ISPV(L)
      CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,L)
      CSS(7,J,L,1)=CSS(7,J,L,1)+WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
    END DO
  END IF
END IF
B=DABS(PV(1,N))
CSSS(1,J)=CSSS(1,J)+WF/B
CSSS(2,J)=CSSS(2,J)+WF*SP(5,L)/B
CSSS(3,J)=CSSS(3,J)+WF*SP(5,L)*PV(2,N)/B
!--this assumes that any flow normal to the x direction is in the y direction
CSSS(4,J)=CSSS(4,J)+WF*SP(5,L)*A/B
IF (ISPR(1,L) > 0) THEN
  CSSS(5,J)=CSSS(5,J)+WF*PROT(N)/B
  CSSS(6,J)=CSSS(6,J)+WF*ISPR(1,L)/B
END IF
!$omp end critical(reflect1)
!
CALL ZGF(RANF,IDT)
IF (FSPEC(J) > RANF) THEN !--specular reflection
  X=2.D00*XB(J)-X
  PV(1,N)=-PV(1,N)
  DTR=(X-XB(J))/PV(1,N)
ELSE                      !--diffuse reflection
  VMPS=SQRT(2.D00*BOLTZ*TSURF(J)/SP(5,L))
  DTR=(X-XB(J))/PV(1,N)   !--isebasti: original line DTR=(XB(J)-PX(N))/PV(1,N) was corrected
  !write(*,*) n,px(n),x-px(n),x,dtm,dtr
  CALL ZGF(RANF,IDT)
  PV(1,N)=SQRT(-LOG(RANF))*VMPS
  IF (J == 2) PV(1,N)=-PV(1,N)
  CALL RVELC(PV(2,N),PV(3,N),VMPS,IDT)
  PV(2,N)=PV(2,N)+VSURF(J)
  IF (ISPR(1,L) > 0) CALL SROT(L,TSURF(J),PROT(N),IDT)
  IF (MMVM > 0) THEN
    DO K=1,ISPV(L)
      CALL SVIB(L,TSURF(J),IPVIB(K,N),K,IDT)
    END DO
  END IF
END IF
!
!$omp critical(reflect2)
CSS(2,J,L,2)=CSS(2,J,L,2)-WF*PV(1,N)*SP(5,L)
CSS(3,J,L,2)=CSS(3,J,L,2)-WF*(PV(2,N)-VSURF(J))*SP(5,L)
CSS(4,J,L,2)=CSS(4,J,L,2)-WF*PV(3,N)*SP(5,L)
A=(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
CSS(5,J,L,2)=CSS(5,J,L,2)-WF*0.5D00*SP(5,L)*A
IF (ISPR(1,L) > 0) CSS(6,J,L,2)=CSS(6,J,L,2)-WF*PROT(N)
IF (MMVM > 0) THEN
  IF (ISPV(L).GT.0) THEN
    DO K=1,ISPV(L)
      CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,L)
      CSS(7,J,L,2)=CSS(7,J,L,2)-WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L) !vibrational energy accumulation
    END DO
  END IF
END IF
B=DABS(PV(1,N))
CSSS(1,J)=CSSS(1,J)+WF/B
CSSS(2,J)=CSSS(2,J)+WF*SP(5,L)/B
CSSS(3,J)=CSSS(3,J)+WF*SP(5,L)*PV(2,N)/B
!--this assumes that any flow normal to the x direction is in the y direction
CSSS(4,J)=CSSS(4,J)+WF*SP(5,L)*A/B
IF (ISPR(1,L) > 0) THEN
  CSSS(5,J)=WF*CSSS(5,J)+PROT(N)/B
  CSSS(6,J)=CSSS(6,J)+WF*ISPR(1,L)/B
END IF
!$omp end critical(reflect2)
!
XI=XB(J)
DX=DTR*PV(1,N)
DZ=0.D00
IF (IFX > 0) DY=DTR*PV(2,N)
IF (IFX == 2) DZ=DTR*PV(3,N)
IF (IFX == 0) X=XI+DX
IF (IFX > 0) CALL AIFX(XI,DX,DY,DZ,X,PV(1,N),PV(2,N),PV(3,N),IDT)
!write(*,*) n,xi,x-xi,x,dtm,dtr
!pause
!
RETURN
!
END SUBROUTINE REFLECT
!
!*****************************************************************************
!
SUBROUTINE RBC(XI,DX,DY,DZ,R,S)
!
!--calculates the trajectory fraction S from a point at radius XI with
!----displacements DX, DY, and DZ to a possible intersection with a
!----surface of radius R, IFX=1, 2 for cylindrical, spherical geometry
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
!
IMPLICIT NONE
!
REAL(KIND=8) :: A,B,C,XI,DX,DY,DZ,R,S,DD,S1,S2
!
DD=DX*DX+DY*DY
IF (IFX == 2) DD=DD+DZ*DZ
B=XI*DX/DD
C=(XI*XI-R*R)/DD
A=B*B-C
IF (A >= 0.D00) THEN
!--find the least positive solution to the quadratic
  A=DSQRT(A)
  S1=-B+A
  S2=-B-A
  IF (S2 < 0.D00) THEN
    IF (S1 > 0.D00) THEN
      S=S1
    ELSE
      S=2.D00
    END IF
  ELSE IF (S1 < S2) THEN
    S=S1
  ELSE
    S=S2
  END IF
ELSE
  S=2.D00
!--setting S to 2 indicates that there is no intersection
END IF
!
RETURN
!
END SUBROUTINE RBC
!
!*****************************************************************************
!
SUBROUTINE AIFX(XI,DX,DY,DZ,X,U,V,W,IDT)
!
!--calculates the new radius and realigns the velocity components in
!----cylindrical and spherical flows
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: IDT !included IDT
REAL(KIND=8) :: A,B,C,XI,DX,DY,DZ,X,U,V,W,DR,VR,S,RANF !--isebasti: included RANF
!
IF (IFX == 1) THEN
  DR=DY
  VR=V
ELSE IF (IFX == 2) THEN
  DR=DSQRT(DY*DY+DZ*DZ)
  VR=DSQRT(V*V+W*W)
END IF
A=XI+DX
X=DSQRT(A*A+DR*DR)
S=DR/X
C=A/X
B=U
U=B*C+VR*S
V=-B*S+VR*C
IF (IFX == 2) THEN
  VR=V
  CALL ZGF(RANF,IDT)
  A=DPI*RANF
  V=VR*DSIN(A)
  W=VR*DCOS(A)
END IF
!
RETURN
!
END SUBROUTINE AIFX
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
!
!*****************************************************************************
!
SUBROUTINE INELASTIC_VT(LS, MS, EV, EC, DOF, P)
!
!--calculate P/P_max for VT process
! energy should be given in si
! DOF is the total dof of all mode in EC other than vibrational mode, only need for LB
! EXPCOL type collision model only handles VT
!
USE GAS,ONLY : INONVHS
USE EXPCOL,ONLY : EXPCOL_VT
IMPLICIT NONE
INTEGER, INTENT(IN) :: LS, MS
REAL(KIND=8), INTENT(IN) :: EV, EC, DOF
REAL(KIND=8) :: P

IF (INONVHS(LS,MS) == 0 .OR. INONVHS(LS, MS) == 1) THEN
  ! use regular LB model
  P = (1.d0 - EV/EC)**(DOF*0.5d0 - 1.0d0)  ! bird eq 5.61
ELSE
  P = EXPCOL_VT(LS, MS, EV, EC)
END IF
!
END SUBROUTINE INELASTIC_VT
!
!*****************************************************************************
!
SUBROUTINE INELASTIC_RT(KS, JS, EC, P, IDT)
!
!--calculate P=Erot(KS)/Ec for RT procedure
!
USE GAS,ONLY : INONVHS, SPM, ISPR
USE EXPCOL,ONLY : EXPCOL_RT
IMPLICIT NONE
INTEGER, INTENT(IN) :: KS, JS
INTEGER :: IDT
REAL(KIND=8), INTENT(IN) :: EC
REAL(KIND=8) :: P, RANF, PMAX

IF (INONVHS(KS, JS) == 0 .OR. INONVHS(KS, JS) == 1) THEN
  ! use regular LB model
  IF (ISPR(1,KS) == 2) THEN
    CALL ZGF(RANF, IDT)
    P = 1.0d0 - RANF**(1.0d0 / (2.5d0 - SPM(3,KS,JS))) !EQN 5.46
  ELSE
    CALL LBS(DBLE(ISPR(1,KS))*0.5d0 - 1.0d0, 1.5d0 - SPM(3,KS,JS), P, IDT)
  END IF
ELSE
  P =  EXPCOL_RT(KS,JS,EC,ISPR(1,KS),IDT)
END IF

END SUBROUTINE INELASTIC_RT
!
!*****************************************************************************
!
SUBROUTINE LBS(XMA,XMB,ERM,IDT)
!
!--selects a Larsen-Borgnakke energy ratio using eqn (11.9)
!
IMPLICIT NONE
!
REAL(KIND=8) :: PROB,ERM,XMA,XMB,RANF
INTEGER :: I,N,IDT !--isebasti: included IDT
!
!--I is an indicator
!--PROB is a probability
!--ERM ratio of rotational to collision energy
!--XMA degrees of freedom under selection-1
!--XMB remaining degrees of freedom-1
!
I=0
DO WHILE (I == 0)
  CALL ZGF(RANF,IDT)
  ERM=RANF
  IF ((XMA < 1.D-6).OR.(XMB < 1.D-6)) THEN
!    IF (XMA < 1.E-6.AND.XMB < 1.E-6) RETURN
!--above can never occur if one mode is translational
    IF (XMA < 1.D-6) PROB=(1.D00-ERM)**XMB
    IF (XMB < 1.D-6) PROB=(1.D00-ERM)**XMA
  ELSE
    PROB=(((XMA+XMB)*ERM/XMA)**XMA)*(((XMA+XMB)*(1.D00-ERM)/XMB)**XMB)
  END IF
  CALL ZGF(RANF,IDT)
  IF (PROB > RANF) I=1
END DO
!
RETURN
!
END SUBROUTINE LBS
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
USE MFDSMC,only: IMF, IMFS, NMFEV0 , NMFER0 , NMFET0 ,&
  NMFEV , NMFER , NMFET , NMFEVR , NMFERR , NMFETR ,&
  NMFVT0 , NMFVT , NMFVTR
!
IMPLICIT NONE
!
INTEGER :: N
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
 CCELL(4,:)=SQRT(2.D00*BOLTZ*VAR(8,:)/SP(5,3))*SPM(2,3,3) !--isebasti: included

IF (nonVHS .ne. 0)THEN
  CCELL(4,:) = CCELL(4,:)*1.2d0
ENDIF

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
!
IF (MNRE > 0  .and. IMF == 0 .and. IREAC == 2 .and. IMFS == 1) THEN
   NMFEV0 = 0.; NMFER0 = 0.; NMFET0 = 0.
   NMFEV = 0.; NMFER = 0.; NMFET = 0.
   NMFEVR = 0.; NMFERR = 0.; NMFETR = 0.
   NMFVT0 = 0.; NMFVT = 0.; NMFVTR = 0.
   !NCANGLE = 0; NCRANGLE = 0.
 END IF
!
END SUBROUTINE INITIALISE_SAMPLES
!
!*****************************************************************************
SUBROUTINE SAMPLE_FLOW
!
!--sample the flow properties
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: NC,NCC,LS,N,M,K,L,I,KV
REAL(KIND=8) :: A,TE,TT,WF,EVIB
!
!--NC the sampling cell number
!--NCC the collision cell number
!--LS the species code
!--N,M,K working integers
!--TE total translational energy
!
NSAMP=NSAMP+1
TSAMP=FTIME
!
!$omp parallel &
!$omp private(n,ncc,nc,wf,ls,m,k,evib) &
!$omp reduction(+:cs)
!$omp do
DO N=1,NM
  NCC=IPCELL(N)
  NC=ICCELL(3,NCC)
  WF=1.D00
  IF (IWF == 1) WF=1.D00+WFM*PX(N)**IFX
  IF ((NC > 0).AND.(NC <= NCELLS)) THEN
    IF (MSP > 1) THEN
      LS=ABS(IPSP(N))
    ELSE
      LS=1
    END IF
    CS(0,NC,LS)=CS(0,NC,LS)+1.D00
    CS(1,NC,LS)=CS(1,NC,LS)+WF
    DO M=1,3
      CS(M+1,NC,LS)=CS(M+1,NC,LS)+WF*PV(M,N)
      CS(M+4,NC,LS)=CS(M+4,NC,LS)+WF*PV(M,N)**2
    END DO
    IF (MMRM > 0) CS(8,NC,LS)=CS(8,NC,LS)+WF*PROT(N)
    IF (MMVM > 0) THEN
      IF (ISPV(LS) > 0) THEN
        DO K=1,ISPV(LS)
          CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,LS)
          CS(K+8,NC,LS)=CS(K+8,NC,LS)+WF*EVIB !DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,LS) !vibrational energy accumulation
        END DO
      END IF
    END IF
  ELSE
    WRITE (*,*) 'Illegal sampling cell',NC,NCC,' for MOL',N,' at',PX(N)  !;STOP
  END IF
END DO
!$omp end do
!$omp end parallel
!
A=0.d0
IF ((IENERS > 0).AND.(GASCODE == 8)) CALL ENERGY(GASCODE,0,A)
ENERS=ENERS+A
!
RETURN
END SUBROUTINE SAMPLE_FLOW
!
!*****************************************************************************
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
USE MFDSMC
!
IMPLICIT NONE
!
INTEGER :: I,IJ,J,JJ,K,KK,JR,KR,L,M,N,NN,NMCR,CTIME,NMS(MSP),ISNAP(1,NSNAP),ZCHECK,NSS,MCT,MAXLEVEL,IDT=0 !--isebasti: included ZCHECK,NSS,IDT
INTEGER(KIND=8) :: NNN
REAL(KIND=8) :: AA,BB,BB2,CC,AS,AT,AU,AQ  !--isebasti: included RANF,GAM
REAL(KIND=8) :: A,B,C,C2,D,SDTM,SMCR,DOF,AVW,UU,VDOFM,TVIBM,DTMI,TT,EVIBT,RANF,F(MMVM,0:100,MSP),TVPDF(MMVM,MSP)  !--isebasti: included RANF,F,TVPDF
REAL(KIND=8) :: VPF(MMVM,MSP),EQROT(MSP,100),EQVIBVT(MMVM,MSP,0:100),EQVIBOT(MMVM,MSP,0:100)  !--isebsati: included
REAL(KIND=8) :: DSUM(0:12),SUMS(0:8,2),TV(MMVM,MSP),TVIB(MSP),DF(NCELLS,MMVM,MSP),VDOF(MSP),PPA(MSP),&
                THCOL(MSP,MSP),PSNAP(1,NSNAP),VSNAP(3,NSNAP),SDOF(MSP),SRR(MNRE),EVIB,QVIBVT,QVIBOT,COEF(0:9),ROOTF
INTEGER,PARAMETER :: PQKIND = SELECTED_REAL_KIND(8)
REAL(KIND=PQKIND) :: AQUAD
CHARACTER (LEN=32) :: FILENAME,TNAME
CHARACTER (LEN=4) :: E
REAL(8) :: COLRR(MNRE)
REAL(8),EXTERNAL :: ERF,gamain,GAM
!
!--CTIME  computer time (microseconds)
!--SUMS(N,L) sum over species of CSS(N,J,L,M) for surface properties
!
!--For flowfield properties,where <> indicates sampled sum
!--DSUM(0) the molecular number sum over all species
!--DSUM(1) the weighted number sum over all species
!--DSUM(2) the weighted sum of molecular masses
!--DSUM(3),(4),(5) the weighted sum over species of m*<u>,<v>,<w>
!--DSUM(6) the weighted sum over species of m*(<u**2>+<v**2>+<w**2>)
!--DSUM(7) the weighted sum over species of <u**2>+<v**2>+<w**2>
!--DSUM(8) the weighted sum of rotational energy
!--DSUM(9) the weighted sum of rotational degrees of freedom
!--DSUM(10) the weighted sum over species of m*<u**2>
!--DSUM(11) the weighted sum over species of m*<v**2>
!--DSUM(12) sum over species of m*<w**2>
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
!--SDOF(L) total degrees of freedom for species L
!--PPA particles per atom
!--NMS number per species
!--NSNAP number of particles considered in SNAPSHOT files
!--SRR sampled reaction rate
!
!----------------------------------------------------------------------------
!--set variables
!
NOUT=NOUT+1
IF (NOUT > 999999999) NOUT=NOUT-999999999
! CALL NUMCHAR4 (NOUT,E)
WRITE (*,*) 'Generating files for output interval',NOUT
WRITE (*,*) 'ISF,FTIME,Number of samples',ISF,FTIME,NSAMP
!
!----------------------------------------------------------------------------
!--compute surface properties
!
VARSP=0.D00
IF (IFX == 0) A=FNUM/(FTIME-TISAMP)   !--flow X-section area = unity for 1-D flow
DO JJ=1,2
IF (IFX == 1) A=FNUM/(2.D00*PI*XB(JJ)*(FTIME-TISAMP))
IF (IFX == 2) A=FNUM/(4.D00*PI*XB(JJ)*XB(JJ)*(FTIME-TISAMP))
!--JJ=1 for surface at XB(1), JJ=2 for surface at XB(2)
  IF (ITYPE(JJ) == 2) THEN
    SUMS=0.D00
    DO L=1,MSP
      DO J=0,8
        DO IJ=1,2
          SUMS(J,IJ)=SUMS(J,IJ)+CSS(J,JJ,L,IJ)
        END DO
      END DO
    END DO
!
    VARS(0,JJ)=SUMS(0,1)
    VARS(1,JJ)=SUMS(1,1)
    VARS(2,JJ)=SUMS(1,2)
    VARS(3,JJ)=SUMS(1,1)*A
    VARS(4,JJ)=SUMS(1,2)*A
    VARS(5,JJ)=SUMS(2,1)*A
    VARS(6,JJ)=SUMS(2,2)*A
    VARS(7,JJ)=SUMS(3,1)*A
    VARS(8,JJ)=SUMS(3,2)*A
    VARS(9,JJ)=SUMS(4,1)*A
    VARS(10,JJ)=SUMS(4,2)*A
    VARS(11,JJ)=SUMS(5,1)*A
    VARS(12,JJ)=SUMS(5,2)*A
    VARS(13,JJ)=SUMS(6,1)*A
    VARS(14,JJ)=SUMS(6,2)*A
    VARS(15,JJ)=SUMS(7,1)*A
    VARS(16,JJ)=SUMS(7,2)*A
    VARS(17,JJ)=SUMS(8,1)*A
    VARS(18,JJ)=SUMS(8,2)*A
    IF (CSSS(1,JJ) > 1.D-6) THEN
      VARS(19,JJ)=CSSS(3,JJ)/CSSS(2,JJ)
      VARS(20,JJ)=(CSSS(4,JJ)-CSSS(2,JJ)*VARS(19,JJ)*VARS(19,JJ))/(CSSS(1,JJ)*3.D00*BOLTZ)-TSURF(JJ)
      VARS(19,JJ)=VARS(19,JJ)-VSURF(JJ)
      IF (CSSS(6,JJ) > 0.D00) THEN
        VARS(21,JJ)=(2.D000/BOLTZ)*(CSSS(5,JJ)/CSSS(6,JJ))-TSURF(JJ)
      ELSE
        VARS(21,JJ)=0.D00
      END IF
    ELSE
      VARS(19,JJ)=0.D00
      VARS(20,JJ)=0.D00
      VARS(21,JJ)=0.D00
    END IF
    VARS(22,JJ)=(SUMS(2,1)+SUMS(2,2))*A
    VARS(23,JJ)=(SUMS(3,1)+SUMS(3,2))*A
    VARS(24,JJ)=(SUMS(4,1)+SUMS(4,2))*A
    VARS(25,JJ)=(SUMS(5,1)+SUMS(5,2))*A
    VARS(26,JJ)=(SUMS(6,1)+SUMS(6,2))*A
    VARS(27,JJ)=(SUMS(7,1)+SUMS(7,2))*A
    VARS(28,JJ)=(SUMS(8,1)+SUMS(8,2))*A
    VARS(29,JJ)=VARS(11,JJ)+VARS(13,JJ)+VARS(15,JJ)
    VARS(30,JJ)=VARS(12,JJ)+VARS(14,JJ)+VARS(16,JJ)
    VARS(31,JJ)=VARS(29,JJ)+VARS(30,JJ)
    DO L=1,MSP
      IF (SUMS(1,1) > 0) THEN
        VARS(32+L,JJ)=100.*CSS(1,JJ,L,1)/SUMS(1,1)
      ELSE
        VARS(32+L,JJ)=0.
      END IF
    END DO
  END IF
END DO
!
!----------------------------------------------------------------------------
!--compute flowfield properties
VAR=0.D00
VARSP=0.
SMCR=0
NMCR=0
VDOFM=0.
DO N=1,NCELLS
  A=FNUM/(CELL(4,N)*NSAMP)
  IF (IVB == 1) A=A*((XB(2)-XB(1))/(XB(2)+VELOB*0.5D00*(FTIME+TISAMP)-XB(1)))**(IFX+1)
!--check the above for non-zero XB(1)
  DSUM=0.
  NMCR=NMCR+1
  DO L=1,MSP
    DSUM(0)=DSUM(0)+CS(0,N,L)
    DSUM(1)=DSUM(1)+CS(1,N,L)
    DSUM(2)=DSUM(2)+SP(5,L)*CS(1,N,L)
    DO K=1,3
      DSUM(K+2)=DSUM(K+2)+SP(5,L)*CS(K+1,N,L)
      IF (CS(1,N,L) > 0.1D00) THEN
        VARSP(K+1,N,L)=CS(K+4,N,L)/CS(1,N,L)
!--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
        VARSP(K+8,N,L)=CS(K+1,N,L)/CS(1,N,L)
!--VARSP(9,10,11 are temporarily the mean of the velocities
      END IF
    END DO
    DSUM(6)=DSUM(6)+SP(5,L)*(CS(5,N,L)+CS(6,N,L)+CS(7,N,L))
    DSUM(10)=DSUM(10)+SP(5,L)*CS(5,N,L)
    DSUM(11)=DSUM(11)+SP(5,L)*CS(6,N,L)
    DSUM(12)=DSUM(12)+SP(5,L)*CS(7,N,L)
    IF (CS(1,N,L) > 0.5D00) THEN
      DSUM(7)=DSUM(7)+CS(5,N,L)+CS(6,N,L)+CS(7,N,L)
    END IF
    IF (ISPR(1,L) > 0) THEN
      DSUM(8)=DSUM(8)+CS(8,N,L)
      DSUM(9)=DSUM(9)+CS(1,N,L)*ISPR(1,L)
    END IF
  END DO
  AVW=0.
  DO L=1,MSP
    VARSP(0,N,L)=CS(1,N,L)
    VARSP(1,N,L)=0.D00
    VARSP(6,N,L)=0.
    VARSP(7,N,L)=0.
    VARSP(8,N,L)=0.
    IF (DSUM(1) > 0.1) THEN
      VARSP(1,N,L)=CS(1,N,L)/DSUM(1)  !isebasti: deleted 100* factor
      AVW=AVW+SP(3,L)*CS(1,N,L)/DSUM(1)
      IF ((ISPR(1,L) > 0).AND.(CS(1,N,L) > 0.5)) VARSP(6,N,L)=(2.D00/BOLTZ)*CS(8,N,L)/(DFLOAT(ISPR(1,L))*CS(1,N,L))
    END IF
    VARSP(5,N,L)=0.
    DO K=1,3
      VARSP(K+1,N,L)=(SP(5,L)/BOLTZ)*(VARSP(K+1,N,L)-VARSP(K+8,N,L)**2)
      VARSP(5,N,L)=VARSP(5,N,L)+VARSP(K+1,N,L)
    END DO
    VARSP(5,N,L)=VARSP(5,N,L)/3.D00
    VARSP(8,N,L)=(3.D00*VARSP(5,N,L)+DFLOAT(ISPR(1,L))*VARSP(6,N,L))/(3.D00+DFLOAT(ISPR(1,L))) !isebasti: included according to DSMC.f90
  END DO
!
  IF (IVB == 0) VAR(1,N)=CELL(1,N)
  IF (IVB == 1) THEN
    C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV)      !--new DDIV
    VAR(1,N)=XB(1)+(DFLOAT(N-1)+0.5)*C
  END IF
  VAR(2,N)=DSUM(0)
  IF (DSUM(1) > 0.5) THEN
    VAR(3,N)=DSUM(1)*A    !--number density Eqn. (4.28)
    VAR(4,N)=VAR(3,N)*DSUM(2)/DSUM(1)   !--density Eqn. (4.29)
    VAR(5,N)=DSUM(3)/DSUM(2)    !--u velocity component Eqn. (4.30)
    VAR(6,N)=DSUM(4)/DSUM(2)    !--v velocity component
    VAR(7,N)=DSUM(5)/DSUM(2)    !--w velocity component
    UU= VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2
    IF (DSUM(1) > 1) THEN
      VAR(8,N)=(ABS(DSUM(6)-DSUM(2)*UU))/(3.D00*BOLTZ*DSUM(1))  !Eqn. (4.39)
!--translational temperature
      VAR(19,N)=(ABS(DSUM(10)-DSUM(2)*VAR(5,N)**2))/(BOLTZ*DSUM(1))
      VAR(20,N)=(ABS(DSUM(11)-DSUM(2)*VAR(6,N)**2))/(BOLTZ*DSUM(1))
      VAR(21,N)=(ABS(DSUM(12)-DSUM(2)*VAR(7,N)**2))/(BOLTZ*DSUM(1))
    ELSE
      VAR(8,N)=1.
      VAR(19,N)=1.
      VAR(20,N)=1.
      VAR(21,N)=1.
    END IF
!--rotational temperature
    IF (DSUM(9) > 0.01D00) THEN
      VAR(9,N)=(2.D00/BOLTZ)*DSUM(8)/DSUM(9)    !Eqn. (4.36)
    ELSE
      VAR(9,N)=0.
    END IF
    DOF=(3.D00+DSUM(9)/DSUM(1))
!--vibration temperature default
    VAR(10,N)=FTMP(1)
!--overall temperature based on translation and rotation
    VAR(11,N)=(3.*VAR(8,N)+(DSUM(9)/DSUM(1))*VAR(9,N))/DOF
!--scalar pressure (now (from V3) based on the translational temperature)
    VAR(18,N)=VAR(3,N)*BOLTZ*VAR(8,N)
!
!--Tvib calculations according to DSMC.f90
    IF (MMVM > 0) THEN
      DO L=1,MSP
        VDOF(L)=0.
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (CS(K+8,N,L) > 0.d0) THEN
              AQUAD=CS(K+8,N,L)/CS(1,N,L)   !average vibrational energy (J)
              EVIB=AQUAD/EVOLT              !average vibrational energy (eV)
              IF (IVMODEL(L,1) == 0) THEN
                TV(K,L)=SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/AQUAD)  !--Eqn.(4.45) - assuming SHO
              ELSE
!do i=0,100
!aquad=VIBEN(0,L)*.9999+(0.09745000001d0*EVOLT-VIBEN(0,L))*DFLOAT(i)/100.d0
                  IF (AQUAD >= VIBEN(0,L)) THEN
                    !use secant method to calculate Tvib(evib) based on Boltzmann distribution
                    CALL ROOTF_SECANT(1,L,L,L,1.d3,AQUAD,TV(K,L))
                    !DF(N,K,L)=2.d0*AQUAD/(BOLTZ*TV(K,L)) !--Eqn. (11.28) Bird94 - general definition
                    DF(N,K,L)=2.d0*(SPVM(1,K,L)/TV(K,L))/(DEXP(SPVM(1,K,L)/TV(K,L))-1.d0) !SHO model
                  ELSE
                    TV(K,L)=1.d0 !to avoid division by zero
                    DF(N,K,L)=0.d0
                  END IF
!write(*,*) AQUAD/EVOLT,TV(K,L),SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/DFLOAT(AQUAD))
!pause
!end do
!stop
              END IF
            ELSE
              TV(K,L)=0.
              DF(N,K,L)=0.
            END IF
            VDOF(L)=VDOF(L)+DF(N,K,L)  !--Eqn.(4.49)
          END DO
          TVIB(L)=0.
          DO K=1,ISPV(L)
            IF (VDOF(L) > 1.D-6) THEN
              TVIB(L)=TVIB(L)+TV(K,L)*DF(N,K,L)/VDOF(L)  !--Eqn.(4.50)
            ELSE
              TVIB(L)=SUM(TV(:,L))/DFLOAT(ISPV(L))
            END IF
          END DO
        ELSE
          TVIB(L)=0. !TREF  !--isebasti: TREF is not defined
          VDOF(L)=0.
        END IF
        VARSP(7,N,L)=TVIB(L)
      END DO
      VDOFM=0.
      TVIBM=0.
      A=1.D00 !--isebasti: instead of 0
      DO L=1,MSP
        IF (ISPV(L) > 0) A=A+CS(1,N,L)
      END DO
      DO L=1,MSP
        IF (ISPV(L) > 0) THEN
          VDOFM=VDOFM+VDOF(L)*CS(1,N,L)/A  !--Eqn.(4.51)
          TVIBM=TVIBM+TVIB(L)*CS(1,N,L)/A  !--Eqn.(4.52)
        END IF
      END DO
      VAR(10,N)=TVIBM
    END IF
!
!--convert the species velocity components to diffusion velocities
    DO L=1,MSP
      IF (VARSP(0,N,L) > 0.5) THEN
        DO K=1,3
          VARSP(K+8,N,L)=VARSP(K+8,N,L)-VAR(K+4,N)
        END DO
      ELSE
        DO K=1,3
          VARSP(K+8,N,L)=0.D00
        END DO
      END IF
    END DO
!
!--reset the overall temperature and degrees of freedom (now including vibrational modes)
    IF (MMVM > 0) THEN
      DO L=1,MSP
        SDOF(L)=3.D00+ISPR(1,L)+VDOF(L)
        VARSP(8,N,L)=(3.*VARSP(5,N,L)+ISPR(1,L)*VARSP(6,N,L)+VDOF(L)*VARSP(7,N,L))/SDOF(L)  !species overall T
      END DO
      A=0.D00
      B=0.D00
      DO L=1,MSP
        A=A+SDOF(L)*VARSP(8,N,L)*CS(1,N,L)
        B=B+SDOF(L)*CS(1,N,L)
      END DO
      VAR(11,N)=A/B !mixture overall T
      DOF=DOF+VDOFM !--isebasti: included
    END IF
!
!--Mach number
    VAR(17,N)=DSQRT(VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2)
    VAR(12,N)=VAR(17,N)/SQRT((DOF+2.D00)*VAR(11,N)*(DSUM(1)*BOLTZ/DSUM(2))/DOF)
!--average number of molecules in (collision) cell
    VAR(13,N)=DSUM(0)/NSAMP/DFLOAT(NCIS)
    IF (COLLS(N) > 2.) THEN
!--mean collision time (see page 17 of my Bird94 book)
      VAR(14,N)=0.5D00*(FTIME-TISAMP)*(DSUM(1)/NSAMP)/WCOLLS(N)
!--mean free path (based on r.m.s speed with correction factor based on equilibrium)
      VAR(15,N)=0.92132D00*DSQRT(DABS(DSUM(7)/DSUM(1)-UU))*VAR(14,N)
      VAR(16,N)=CLSEP(N)/(COLLS(N)*VAR(15,N))
    ELSE
!--m.f.p set by nominal values
      VAR(14,N)=1.D10
      VAR(15,N)=1.D10/VAR(3,N)
    END IF
  ELSE
    DO L=3,19
      VAR(L,N)=0.
    END DO
  END IF
!
  IF(N == NSPDF) THEN
    TVPDF(:,:)=TV(:,:)                        !store TV for the cell with pdf sampling
    A=2.*PCOLLS*DFLOAT(NSAMP)/DSUM(1)         !mean collision times
    MCT=A                                     !integer value
!
    IF(VARSP(1,N,1) > VARSP(1,N,3)) THEN
      L=1                                     !for N2 species
    ELSE
      L=3                                     !for O2 species
    END IF
!
    B=0.d0
    QVIBVT=0.d0
    DO M=0,IVMODEL(L,2)
      CALL VIB_ENERGY(EVIBT,M,1,L)
      B=B+EVIBT*DEXP(-EVIBT/(BOLTZ*VAR(8,N)))
      QVIBVT=QVIBVT+DEXP(-EVIBT/(BOLTZ*VAR(8,N)))
    END DO
    EVIBT=(B/QVIBVT)/EVOLT                     !equilibrium vibrational temperature based Tt (ev)
!
    B=1.d0/(VAR(3,N)*VARSP(1,N,4)*CSCR(L,4)/TCOL(L,4)) !mct between collisions of one particular L molecule with any O
    C=1.d0/(VAR(3,N)*VARSP(1,N,L)*CSCR(L,L)/TCOL(L,L)) !mct between collisions of one particular L molecule with any L
!
    OPEN (7,FILE='RELAX.DAT',ACCESS='APPEND') !write temperature relaxation
    WRITE (7,993) FTIME,VAR(8:11,N),VAR(18,N),DSUM(9)/DSUM(1),VDOFM,NM,A,&
                  DABS(0.5*(VAR(8,N)+VAR(9,N))-VAR(10,N)),DABS(VAR(8,N)-VAR(9,N)),DABS(VAR(8,N)-VAR(10,N)),&
                  (TV(:,L),L=1,MSP),EVIB,EVIBT,VAR(14,N),B,C
    CLOSE (7)
  END IF
!
END DO
!
!----------------------------------------------------------------------------
!--write general/surface properties
!
OPEN (3,FILE='DS1GEN.DAT')
!
IF (IFX == 0) WRITE (3,*) 'DSMC program DS1 for a one-dimensional plane flow'
IF (IFX == 1) WRITE (3,*) 'DSMC program DS1 for a cylindrical flow'
IF (IFX == 2) WRITE (3,*) 'DSMC program DS1 for a spherical flow'
WRITE (3,*)
!
WRITE (3,*) 'Interval',NOUT,'Time ',FTIME, ' with',NSAMP,' samples from',TISAMP
WRITE (*,*) 'TOTAL MOLECULES = ',NM
!
NMS=0
DO N=1,NM
  M=IPSP(N)
  NMS(M)=NMS(M)+1
END DO
WRITE (3,*) 'Total simulated molecules =',NM
DO N=1,MSP
  WRITE (*,*) 'SPECIES ',N,' TOTAL = ',NMS(N)
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
CTIME=MCLOCK()
WRITE (3,*) 'Computation time',FLOAT(CTIME)/1000.,'seconds'
WRITE (3,*) 'Collision events per (cpu) second',(TOTCOL-TOTCOLI)*1000.D00/DFLOAT(CTIME)
WRITE (3,*) 'Molecule moves per (cpu) second',(TOTMOV-TOTMOVI)*1000.D00/DFLOAT(CTIME)
WRITE (3,*)
!
IF (IPDF > 0) THEN
  WRITE (3,*) 'Distribution function sampling started at',TISAMP
  WRITE (3,*) 'Total dissociations',NDISSOC
  WRITE (3,*) 'Total recombinations',NRECOMB
  WRITE (3,*) 'Gas temperature',VAR(11,1),' K'
  WRITE (3,*) 'Gas translational temperature',VAR(8,1),' K'
  WRITE (3,*) 'Gas rotational temperature',VAR(9,1),' K'
  WRITE (3,*) 'Gas vibrational temperature',VAR(10,1),' K'
  DO L=1,MSP
    WRITE (3,*) 'Species',L,' overall temperature',VARSP(8,1,L),' K'
    WRITE (3,*) 'Species',L,' translational temperature',VARSP(5,1,L),' K'
    WRITE (3,*) 'Species',L,' rotational temperature',VARSP(6,1,L),' K'
    WRITE (3,*) 'Species',L,' vibrational temperature',VARSP(7,1,L),' K'
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
WRITE (3,994) ' Macro or Coll Temps (ITCV,IEAA,IZV):',ITCV,IEAA,IZV
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
IF (MNRE > 0) THEN
  OPEN (3,FILE='DS1REAC.DAT',ACCESS='APPEND')
!
  IF (JCD == 1) THEN
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
  IF (JCD == 0) THEN
    AS=0.
    SRR=0.
    DO M=1,MNRE
      AS=AS+REAC(M)
    END DO
    IF (AS == 0.) AS=1.d0
    IF (IREAC == 0) THEN
      WRITE (3,993) FTIME,AS,REAC(1:MNRE)/AS  !no reaction rate sampling
    ELSE
      DO K=1,MNRE
        A=REAC(K)*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP) !# of reaction/V/t
        C=TCOL(LE(K),ME(K))*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP)
        IF(KP(K) >= 0) THEN
          ! nA * nB
          B=(VARSP(1,NSPDF,LE(K))*VARSP(1,NSPDF,ME(K))*VAR(3,NSPDF)**2.d0)/(AVOG*1.d3)                             !exchange/dissociation
        ELSE
          B=(VARSP(1,NSPDF,LE(K))*VARSP(1,NSPDF,ME(K))*VARSP(1,NSPDF,MP(K))*VAR(3,NSPDF)**3.d0)/(AVOG*1.d3)**2.d0  !recombination
        END IF
        SRR(K)=A/B !sampled reaction rate (cm3/mol/s)
        COLRR(K) = C/B
        ! COLRR is the collision rate, but not the collision frequency
        ! For different specie, collision frequency = COLRR * n1 * n2
        ! For the same specie, collision frequency = COLRR * n1 * n2 *0.5


       ! IF (VAR(10,NSPDF) < 6.d3)  SRR(K)=SRR(K)/1.d3        !trick to sample low T rates
      END DO
      WRITE (3,993) FTIME,VAR(8:11,NSPDF),VAR(18,NSPDF),AS,REAC(1:MNRE)/AS,SRR(1:MNRE),COLRR(1:MNRE)
    END IF
  END IF
!
  CLOSE(3)

7787 FORMAT('#Reaction: ',I3,' Sp: ',I2,'+',I2,'->',I2,'+',I2,'+',I2)
7786 FORMAT("#",A13,3X,4(A14,3X))
  DO K=1,MNRE
    IF (KP(K) > 0) THEN !only for dissociation reaction
      WRITE(FILENAME,'("REAC_COUNT_",I3.3,".DAT")') K
      IF (K == 1) FILENAME='REAC_COUNT.DAT'
      OPEN(3,FILE=trim(FILENAME),ACCESS = 'APPEND')
      IF (NOUT < 1) THEN
        WRITE(3,7787) K, LE(K),ME(K),KP(K),LP(K),MP(K)
        WRITE(3,7786) 'Nreac','Evrem (K)','Evrem/D','Evrem (eV)','Time'
      END IF
      IF (REAC(K) > 1.0d0) THEN
        WRITE(3,'(5(E14.6,3x))') REAC(K), EVREM(K)/REAC(K)/BOLTZ, EVREM(K)/REAC(K)/REA(2,K), &
          EVREM(K)/EVOLT/REAC(K),FTIME
      ELSE
        WRITE(3,'(5(E14.6,3x))') 0.d0, 0.d0, 0.d0, 0.d0, FTIME
      END IF
      CLOSE(3)
    END IF
  END DO

END IF
!
!
!---------------------------------------------------------------------------
!--write debug information for imf method
7788 FORMAT('ZONE I=',I6,', T="REAC',I3.3,'", SOLUTIONTIME=',G14.6)
7789 FORMAT(F6.2,3X,21(G14.6,3X))
7790 FORMAT('ZONE I = ',i6,', J = ',i6,', T= "REAC',I3.3,'_',I1,'", DATAPACKING = POINT')
7785 FORMAT(I3,2X,F10.6,2X,7(G14.6,2X))
7784 FORMAT(I3,2X,F10.6,2X,F10.6,2X,14(G14.6,2X))
IF (IREAC == 2 .AND. IMF .ne. 0 .AND. MNRE <= 2.AND. IMFS == 1 .AND. MNRE>0) THEN
  !
  !----- Distribution of Et and Er------------------------
  !
  OPEN (3,FILE='IMF_ETR.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Nomralize by T: ",FTMP0,"K, Current T: ",VAR(8,NSPDF),"K "
  WRITE(3,"(5A)") 'VARIABLES = "E/kT",',&
    '"ET0_N","ET0_P","ET_N","ET_P","ETR_N","ETR_P",',&
    '"ER0_N","ER0_P","ER_N","ER_P","ERR_N","ERR_P",',&
    '"ER01_N","ER01_P","ER1_N","ER1_P","ERR1_N","ERR1_P",',&
    '"PEQ","PCOL","PBOLTZ"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      WRITE(3,7788) 1000,L,FTIME
      K = IMFpair(IREA(1,L),IREA(2,L))

      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE IF (IREA(1,L) > IREA(2,L)) THEN
        IJ = 2; JJ = 1
      END IF

      AA = SUM(NMFET0(:,K))
      BB = SUM(NMFER0(:,IJ,K)); BB2 = SUM(NMFER0(:,JJ,K))

      A = SUM(NMFET(:,K))
      C = SUM(NMFER(:,IJ,K));   C2 = SUM(NMFER(:,JJ,K))

      AS = SUM(NMFETR(:,L))
      AU = SUM(NMFERR(:,1,L)); AT = SUM(NMFERR(:,2,L))


      CC = 2.5d0-SPM(3,IREA(1,L),IREA(2,L))
      DO N = 1,1000
        B = DFLOAT(N)*0.01d0
        D = 0.5d0*SPI*(ERF(dsqrt(B))-ERF(dsqrt(B-0.01d0)))
        !WRITE(3,7789) B,NMFET0(N),NMFET0(N)/AA,&
        !NMFER0(N), NMFER0(N)/BB,&
        !NMFET(L,N),DFLOAT(NMFET(L,N))/A,&
        !NMFER(L,N),DFLOAT(NMFER(L,N))/C, &
        !NMFETR(L,N),DFLOAT(NMFETR(L,N))/AS, &
        !NMFERR(L,N),DFLOAT(NMFERR(L,N))/SRR,&
        !(0.99d0+B)*dexp(-B+0.01d0)-(1.0d0+B)*dexp(-B),&
        !dcosh(B-0.01d0)-dcosh(B)-dsinh(B-0.01d0)+dsinh(B)
        WRITE(3,7789) B,&
          NMFET0(N,K),   NMFET0(N,K)/AA,  NMFET(N,K),    NMFET(N,K)/A,&
          NMFETR(N,L),   NMFETR(N,L)/AS,  NMFER0(N,IJ,K),NMFER0(N,IJ,K)/BB, &
          NMFER(N,IJ,K), NMFER(N,IJ,K)/C, NMFERR(N,1,L), NMFERR(N,1,L)/AU, &
          NMFER0(N,JJ,K),NMFER0(N,JJ,K)/BB2, &
          NMFER(N,JJ,K), NMFER(N,JJ,K)/C2,NMFERR(N,2,L), NMFERR(N,2,L)/AT, &
          gamain(B,1.5d0,NN) - gamain(B-0.01d0,1.5d0,NN),&
          gamain(B,CC,NN) - gamain(B-0.01d0,CC,NN),&
          dcosh(B-0.01d0)-dcosh(B)-dsinh(B-0.01d0)+dsinh(B)
      END DO
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
  !
  ! -- Distribution of Ev
  !
  OPEN (3,FILE='IMF_EV.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Initial T: ",FVTMP0,"K Current T: ",VAR(10,NSPDF),"K"
  WRITE(3,"(3A)")'VARIABLES = "v","Ev (eV)","Ev2 (eV)",', &
    '"Ev0_N","Ev0_P","Ev_N","Ev_P","EvR_N","EvR_P","P0",',&
    '"Ev01_N","Ev01_P","Ev1_N","Ev1_P","EvR1_N","EvR1_P","P1"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      K = IMFpair(IREA(1,L),IREA(2,L))
      ! write chemical reaction
      WRITE(3,*)
      WRITE(3,'("#",I3,"+",I3,"->",I3,"+",I3,"+",I3)') LE(L),ME(L),KP(L),LP(L),MP(L)

      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE
        IJ = 2; JJ = 1
      END IF

      BB = SUM(NMFEV0(:,IJ,K)); AA = SUM(NMFEV(:,IJ,K)); C = SUM(NMFEVR(:,1,L))
      IF (BB .ge. 1.0d0) THEN
        BB = 1.0d0/BB
      ELSE
        BB = 0.0d0
      END IF
      IF (AA .ge. 1.0d0) THEN
        AA = 1.0d0/AA
      ELSE
        AA = 0.0d0
      END IF
      IF (C  .ge. 1.0d0) THEN
        C = 1.0d0/C
      ELSE
        C = 0.0d0
      END IF

      KK = 100
      IF (ISPV(IREA(2,L)) == 0) THEN
        ! this is an atom-diatom dissociation
        ! KK is the maximum vibrational level
        IF (IVMODEL(IREA(1,L),1) == 1) KK = IVMODEL(IREA(1,L),2)
        !--- write zone header
        WRITE(3,7788) KK+1,L,FTIME  !write zone header
        WRITE(3,"(A)") 'PASSIVEVARLIST = [3,11-17]'
        !--- write data
        DO N = 0,KK
          CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
          WRITE(3,7785) N,EVIB/EVOLT, NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB,&
            NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,NMFEVR(N,1,L),NMFEVR(N,1,L)*C,&
            DEXP(-EVIB/BOLTZ/VAR(10,NSPDF))
        END DO
      ELSE
        BB2 = SUM(NMFEV0(:,JJ,K));
        IF (BB2 .ge. 1.0d0) THEN
          BB2 = 1.0d0/BB2
        ELSE
          BB2 = 0.0d0
        END IF
        A = SUM(NMFEV(:,JJ,K))
        IF (A .ge. 1.0d0) THEN
          A = 1.0d0/A
        ELSE
          A = 0.0d0
        END IF
        CC = SUM(NMFEVR(:,2,L))
        IF (CC .ge. 1.0d0) THEN
          CC = 1.0d0/CC
        ELSE
          CC = 0.0d0
        END IF

        KR = 100; JR = 100;
        IF (IVMODEL(IREA(1,L),1) == 1)    KR = IVMODEL(IREA(1,L),2)
        IF (IVMODEL(IREA(2,L),1) == 1)    JR = IVMODEL(IREA(2,L),2)
        KK = MAX(KR,JR)

        !--- write zone header
        WRITE(3,7788) KK+1,L,FTIME
        !--- write data
        DO N=0,MIN(KR,JR)
          CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
          CALL VIB_ENERGY(AS,N,1,IREA(2,L))
          WRITE(3,7784) N,EVIB/EVOLT,AS/EVOLT,&
            NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB, NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,&
            NMFEVR(N,1,L),NMFEVR(N,1,L)*C, DEXP(-EVIB/BOLTZ/VAR(10,NSPDF)),&
            NMFEV0(N,JJ,K),NMFEV0(N,JJ,K)*BB2, NMFEV(N,JJ,K),NMFEV(N,JJ,K)*A,&
            NMFEVR(N,2,L),NMFEVR(N,2,L)*CC, DEXP(-AS/BOLTZ/VAR(10,NSPDF))
        END DO

        IF (KR < JR) THEN
          DO N=KR+1,JR
            CALL VIB_ENERGY(AS,N,1,IREA(2,L))
            WRITE(3,7784) N,0.0d0,AS/EVOLT,&
              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0 ,&
              NMFEV0(N,JJ,K),NMFEV0(N,JJ,K)*BB2, NMFEV(N,JJ,K),NMFEV(N,JJ,K)*A,&
              NMFEVR(N,2,L),NMFEVR(N,2,L)*CC, DEXP(-AS/BOLTZ/VAR(10,NSPDF))
          END DO
        ELSE IF (KR > JR) THEN
          DO N=JR+1,KR
            CALL VIB_ENERGY(EVIB,N,1,IREA(1,L))
            WRITE(3,7784) N,EVIB/EVOLT,0.0d0,&
              NMFEV0(N,IJ,K),NMFEV0(N,IJ,K)*BB, NMFEV(N,IJ,K),NMFEV(N,IJ,K)*AA,&
              NMFEVR(N,1,L),NMFEVR(N,1,L)*C, DEXP(-EVIB/BOLTZ/VAR(10,NSPDF)),&
              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0
          END DO
        END IF
      END IF
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
  !
  ! -- Reaction probability
  OPEN (3,FILE='IMF_Prob.DAT')
  WRITE(3,"(A,F10.3,A,F10.3,A)") "# Nomralize by T: ",FTMP0,"K, Current T:",VAR(8,NSPDF),"K"
  WRITE(3,"(A)")'VARIABLES = "Et/kT", "Ev (eV)","N","NR1","NR2","P1","P2"'
  DO L = 1,MNRE
    IF (KP(L) > 0) THEN ! dissociation alone
      K = IMFpair(IREA(1,L),IREA(2,L))
      IF (IREA(1,L) <= IREA(2,L)) THEN
        IJ = 1; JJ = 2
      ELSE
        IJ = 2; JJ = 1
      END IF

      ! write data for the dissociating particle
      KK = 100
      IF (IVMODEL(IREA(1,L),1) == 1)     KK = IVMODEL(IREA(1,L),2)
      WRITE(3,7790) 1000, KK+1,L,1

      DO N = 0,KK
        DO M = 1,1000
          CALL VIB_ENERGY(C,N,1,IREA(1,L))
          A = 0.0D0;       B = 0.0D0
          IF (NMFVT0(N,M,IJ,K) .NE. 0.0d0)  A = NMFVT(N,M,IJ,K)/NMFVT0(N,M,IJ,K)
          IF (NMFVT(N,M,IJ,K) .NE. 0.0d0)   B = NMFVTR(N,M,1,L) / NMFVT(N,M,IJ,K)

          WRITE(3,7789) DFLOAT(M)*0.01d0, C/EVOLT,&
            NMFVT0(N,M,IJ,K),NMFVT(N,M,IJ,K),NMFVTR(N,M,1,L),A,B
        END DO
      END DO
      WRITE(3,*)

      ! write data for the collider
      KK = 100
      IF (IVMODEL(IREA(2,L),1) == 1)  KK = IVMODEL(IREA(2,L),2)
      WRITE(3,7790) 1000, KK+1,L,2

      DO N = 0,KK
        DO M = 1,1000
          CALL VIB_ENERGY(C,N,1,IREA(2,L))
          A = 0.0D0;       B = 0.0D0
          IF (NMFVT0(N,M,JJ,K) .NE. 0.0d0)  A = NMFVT(N,M,JJ,K)/NMFVT0(N,M,JJ,K)
          IF (NMFVT(N,M,JJ,K) .NE. 0.0d0)   B = NMFVTR(N,M,2,L) / NMFVT(N,M,JJ,K)

          WRITE(3,7789) DFLOAT(M)*0.01d0, C/EVOLT,&
            NMFVT0(N,M,JJ,K),NMFVT(N,M,JJ,K),NMFVTR(N,M,2,L),A,B
        END DO
      END DO
      WRITE(3,*)
    END IF
  END DO
  CLOSE(3)
END IF


!
!----------------------------------------------------------------------------
!--write flowfield overall properties
!
!IF (ISF >= 1) OPEN (3,FILE='DS1FP00_'//E//'.DAT')
!IF (ISF == 0) OPEN (3,FILE='DS1FP00.DAT')
OPEN (3,FILE='DS1FP00.DAT') !to generate only one output file
!
WRITE(TNAME,7791) NOUT
7791 FORMAT('ZONE T = "'i4.4,'"')
!
WRITE (3,994) '# Flowfield Overall Properties at FTIME:',FTIME
WRITE (3,994) '# NCELLS:', NCELLS
WRITE (3,994) '# NSAMP: ', NSAMP
WRITE (3,994) '# X-coord.     Cell   Sample     Number Dens.   Density   u velocity &
            &  v velocity   w velocity   Trans. Temp.   Rot. Temp.   Vib. Temp. &
            &  Temperature  Mach no.     Mols/cell    m.c.t        m.f.p     &
            &  mcs/mfp        speed      Pressure        TTX          TTY          TTZ     &
            &  dtm/mct     <dx/mfp>      Fx          Fy            Fz         Qtransfer &
            &  Species Fractions'
WRITE (3,994) TRIM(TNAME)
WRITE (3,994) 'SOLUTIONTIME =',FTIME
WRITE (3,992) 'STRANDID = ',QCTMODEL
!
DO N=1,NCELLS
  WRITE (3,996) VAR(1,N),N,VAR(2:21,N),DTM/VAR(14,N),(CELL(3,N)-CELL(2,N))/DFLOAT(NCIS)/VAR(15,N),&
                CST(1:4,N)/CST(0,N),VARSP(1,N,1:MSP)
END DO
CLOSE(3)
!
!----------------------------------------------------------------------------
!--write flowfield properties per species
!
DO L=1,MSP
    WRITE(FILENAME,775) L
775 FORMAT('DS1FP',i2.2)
    OPEN (3,FILE=TRIM(FILENAME)//'.DAT')
!
  WRITE (3,994) '# Flowfield Properties Species:',L
  WRITE (3,994) '# NCELLS:', NCELLS
  WRITE (3,994) '# NSAMP: ', NSAMP
  WRITE (3,994) '# X-coord.     Cell   Sample    Fraction     Species TTx   Species TTy  Species TTz &
             & Trans. Temp.  Rot. Temp.   Vib. Temp.   Spec. Temp  u Diff Vel   v Diff Vel   w Diff Vel.'
  DO N=1,NCELLS
    WRITE (3,996) VAR(1,N),N,VARSP(0,N,L),VARSP(1,N,L),VARSP(2:11,N,L)
  END DO
CLOSE(3)
END DO
!
!----------------------------------------------------------------------------
!--write composition of a reacting gas as a function of time
!
OPEN (10,FILE='COMPOSITION.DAT',ACCESS='APPEND')
A=0.d0
AS=NM
IF (IENERS > 0) THEN
  IF (NOUT <= 1) ENERS0=ENERS
  A=ENERS/ENERS0-1.d0
  ENERS0=ENERS
END IF
WRITE (10,993) FTIME,VAR(8:11,NSPDF),VAR(18,NSPDF),NMS(1:MSP)/AS,A,NM,VAR(5,1),VAR(5,NCELLS)
CLOSE (10)
!
!----------------------------------------------------------------------------
!--write pdf files
!
IF(IPDF > 0) THEN
!
!--compute normalized velocity and speed pdfs
!
  PDF=0.d0
  PDFS=0.d0
  DO N=1,NBINS
    DO L=1,MSP
      DO K=1,3
        PDF(N,K)=PDF(N,K)+BINS(N,K,L)/DBINV(3)/BIN(0,K)
        PDFS(N,K,L)=BINS(N,K,L)/DBINV(3)/BINS(0,K,L)
      END DO
      PDF(N,4)=PDF(N,4)+BINS(N,4,L)/DBINC(3)/BIN(0,4)
      PDF(N,5)=PDF(N,5)+BINS(N,5,L)/DBINE(3)/BIN(0,5)
      PDFS(N,4,L)=BINS(N,4,L)/DBINC(3)/BINS(0,4,L)
      PDFS(N,5,L)=BINS(N,5,L)/DBINE(3)/BINS(0,5,L)
    END DO
  END DO
!
!--velocity components
  OPEN (7,FILE='DS1DVEL.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Interval, Equilibrium pdf, U V W pdfs total, U V W pdfs spec1, U V W pdfs spec2 ... '
  DO N=1,NBINS
    A=DBINV(1)+DBINV(3)*(DFLOAT(N)-0.5d0)                !normalized U/VMP
    B=(1.d0/SPI)/DEXP(A*A)                               !Boltzmann distribution
    WRITE (7,993) A,B,PDF(N,1:3),(PDFS(N,1:3,L),L=1,MSP) !pdfs
  END DO
  CLOSE(7)
!
!--speed
  OPEN (7,FILE='DS1DSPD.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Interval, Equilibrium pdf, C pdf total, C pdf spec1, C pdf spec2 ... '
  DO N=1,NBINS
    A=DBINC(1)+DBINC(3)*(DFLOAT(N)-0.5d0)                !normalized C/VMP
    B=((4.d0/SPI)*A*A)/DEXP(A*A)                         !Boltzmann distribution; ftr- Eqn 15.25 or N.5 from Laurendeau
    WRITE (7,993) A,B,PDF(N,4),(PDFS(N,4,L),L=1,MSP)     !pdfs
  END DO
  CLOSE(7)
!
!--thermal (translational) energy
  OPEN (7,FILE='DS1DTEN.DAT',FORM='FORMATTED')
  WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
  WRITE (7,994) '# NSamples total, spec1, spec2 ...:     ', BIN(0,1),BINS(0,1,1:MSP)  !summation values
  WRITE (7,994) '# Energy (eV), U V W pdfs total, U V W pdfs spec1, U V W pdfs spec2 ... '
  DO N=1,NBINS
    A=DBINE(1)+DBINE(3)*(DFLOAT(N)-0.5d0)                !energy in eV
    WRITE (7,993) A,PDF(N,5),(PDFS(N,5,L),L=1,MSP) !pdfs
  END DO
  CLOSE(7)
!
!--rotational energy
  IF (MMRM > 0) THEN
    DO L=1,MSP
      IF (ISPR(1,L) > 0) THEN
        EQROT=1.d0
        DO M=1,100
          A=(DFLOAT(M)-0.5D00)*0.1D00   !--rotational values
          BB=.5*ISPR(1,L)
          !Boltzmann distribution; frot - Eqn 11.21 from Bird94
          EQROT(L,M)=(1.d0/GAM(BB))*(A**(BB-1.d0))*DEXP(-A)*0.1D00
        END DO
        WRITE(FILENAME,777) L
        777 FORMAT('DS1DROT',i2.2)
        OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
        WRITE (7,994) '# NSPDF, X-coord, Tt, Tr, Tv, To: ', NSPDF,VAR(1,NSPDF),VAR(8:11,NSPDF)
        WRITE (7,994) '# Equilibrium based on the rotational temperature'
        WRITE (7,994) '# Species and number of modes:',L,ISPR(1,L)
        WRITE (7,994) '# erot/kT       sample      f_rot    f_rot_equilib    ratio'
        DO M=1,100
          A=(DFLOAT(M)-0.5D00)*0.1D00
          B=DFLOAT(NDROT(L,M))/BINS(0,1,L) !NDROT and EQROT are fractions in interval 0.1, so multiply by 10 to obtain pdfs
          WRITE (7,993) A,DFLOAT(NDROT(L,M)),B*10.D00,EQROT(L,M)*10.D00,B/EQROT(L,M)
        END DO
        CLOSE (7)
      END IF
    END DO
  END IF
!
!--vibrational levels
  F =0.0d0
  IF (MMVM > 0) THEN
    DO I=1,NSCELLS
      J=NSVEC(I) !sampling cell
      DO L=1,MSP
        IF (ISPV(L) > 0 .AND. (L==1 .or. L==3 .or. L==5) ) THEN !plotting only O2 or N2 populations
          K=1 !single vibrational mode
!
          !calculate the vibrational partition function
          EQVIBVT=1.d0
          EQVIBOT=1.d0
          QVIBVT=0.d0
          QVIBOT=0.d0
          MAXLEVEL=IVMODEL(L,2)
          DO M=0,MAXLEVEL
            CALL VIB_ENERGY(EVIB,M,K,L)
            QVIBVT=QVIBVT+DEXP(-EVIB/(BOLTZ*VARSP(7,J,L)))       !based on Tv_mode
            QVIBOT=QVIBOT+DEXP(-EVIB/(BOLTZ*VARSP(8,J,L))) !based on Tov
          END DO
!
          !calculate Boltzmann distribution function
          DO M=0,MAXLEVEL
            CALL VIB_ENERGY(EVIB,M,K,L)
            EQVIBVT(K,L,M)=DEXP(-EVIB/(BOLTZ*VARSP(7,J,L)))/QVIBVT       !based on Tv_mode
            EQVIBOT(K,L,M)=DEXP(-EVIB/(BOLTZ*VARSP(8,J,L)))/QVIBOT !based on Tov
          END DO
!
          !write sampled values
          !WRITE(FILENAME,778) L,NOUT
          !778 FORMAT('DS1DVIB',i2.2,'_',i4.4)
          WRITE(FILENAME,879) L
          879 FORMAT('DS1DVIB',i2.2)
          WRITE(TNAME,779) L,J !NOUT
          779 FORMAT('ZONE T = "',i2.2,'_',i4.4,'"')
          OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
          WRITE (7,994) '# NS, X-coord, Tt, Tr, Tv, To, t, MCT: ', J,VAR(1,J),VAR(8:11,J),FTIME,MCT
          WRITE (7,994) '# Equilibrium based on both vibrational temperature and overall gas temperature'
          WRITE (7,994) '# Species and number of modes:',L,ISPV(L)
          WRITE (7,994) '# level      Ev       sample      f_vib       f_equil(Tvib) ratio Eq.   f_equil(Tov)  ratio'
          WRITE (7,994) TRIM(TNAME)
          WRITE (7,994) 'SOLUTIONTIME =',FTIME
          WRITE (7,992) 'STRANDID = ',QCTMODEL
          !DO M=0,MAXLEVEL
          DO M=0,IVMODEL(L,2)
            F(1:ISPV(L),M,L)=DFLOAT(NDVIB(J,1:ISPV(L),L,M))/DFLOAT(NDVIB(J,0,L,0))
            CALL VIB_ENERGY(EVIB,M,K,L)
            A=EVIB/EVOLT !energy in eV
            WRITE (7,993) M,A,DFLOAT(NDVIB(J,1:ISPV(L),L,M)),F(1:ISPV(L),M,L),&
                          EQVIBVT(1:ISPV(L),L,M),F(1:ISPV(L),M,L)/EQVIBVT(1:ISPV(L),L,M),&
                          EQVIBOT(1:ISPV(L),L,M),F(1:ISPV(L),M,L)/EQVIBOT(1:ISPV(L),L,M)
          END DO
          CLOSE (7)
!
        END IF
      END DO
    END DO
  END IF
!
END IF

!-- vibrational state-specific rates
IF (IREAC == 2 .AND. IMF .ne. 0 .AND. MNRE <= 2.AND. IMFS == 1 .AND. MNRE>0 .AND. NSCELLS==1) THEN
  OPEN(3, FILE='IMF_Vrate.dat')
  WRITE(3,"(A,G14.6,A,F10.3)") '# time:',FTIME, ' VT: ',VAR(10,NSPDF)
  WRITE(3,"(A)") 'VARIABLES = "Evib(eV)","v","Nreac","Rate (cm3/mol/s)"'
  DO L = 1,MNRE
    NN = LE(L)
    IF (KP(L) > 0 .and. IVMODEL(NN,1) == 1) THEN ! dissociation alone
      WRITE(3,*)
      WRITE(3,'("#",I3,"+",I3,"->",I3,"+",I3,"+",I3)') LE(L),ME(L),KP(L),LP(L),MP(L)
      WRITE(3,7788) IVMODEL(NN,2)+1,L,FTIME  !write zone header
      DO N=0,IVMODEL(NN,2)
        CALL VIB_ENERGY(EVIB,N,1,NN)
        A=NMFEVR(N,1,L)*(FNUM/CELL(4,NSPDF))/(FTIME-TISAMP) !# of reaction/V/t
        B=(VARSP(1,NSPDF,LE(L))*F(1,N,NN)*VARSP(1,NSPDF,ME(L))*VAR(3,NSPDF)**2.d0)/(AVOG*1.d3)    !number density

        IF (NMFEVR(N,1,L) == 0) THEN
          WRITE(3,'(F12.6,2X,I3,2X,E14.6,2X,E14.6)') EVIB/EVOLT,N,NMFEVR(N,1,L),0.0d0
        ELSE
          WRITE(3,'(F12.6,2X,I3,2X,E14.6,2X,E14.6)') EVIB/EVOLT,N,NMFEVR(N,1,L),A/B
        END IF
        ! IF (VAR(10,NSPDF) < 6.d3)  SRR(K)=SRR(K)/1.d3        !trick to sample low T rates
      END DO
      write(3,*)
    END IF
  END DO
  CLOSE(3)
END IF
!
!--sampled pre- and post-reaction vibrational levels
IF (MNRE > 0) THEN
!
  J=1 !pre-reaction
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
    B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
    C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
    IF(A == 0) A=1.d0
    IF(B == 0) B=1.d0
    IF(C == 0) C=1.d0
    WRITE(FILENAME,780) KK
    780 FORMAT('DS1VIBL1_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled pre-reaction vibrational levels'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,B,C,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# I, Molec1 Sample and Fraction, Molec2 Sample and Fraction, , Molec3 Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NPVIB(J,KK,1,1:3,I),NPVIB(J,KK,1,1:3,I)/A,&
                      NPVIB(J,KK,2,1:3,I),NPVIB(J,KK,2,1:3,I)/B,&
                      NPVIB(J,KK,3,1:3,I),NPVIB(J,KK,3,1:3,I)/C
    END DO
    CLOSE (7)
  END DO
!
  J=2 !post-reaction
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
    B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
    C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
    IF(A == 0) A=1.d0
    IF(B == 0) B=1.d0
    IF(C == 0) C=1.d0
    WRITE(FILENAME,782) KK
    782 FORMAT('DS1VIBL2_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled post-reaction vibrational levels'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,B,C,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# I, Molec1 Sample and Fraction, Molec2 Sample and Fraction, , Molec3 Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NPVIB(J,KK,1,1:3,I),NPVIB(J,KK,1,1:3,I)/A,&
                      NPVIB(J,KK,2,1:3,I),NPVIB(J,KK,2,1:3,I)/B,&
                      NPVIB(J,KK,3,1:3,I),NPVIB(J,KK,3,1:3,I)/C
    END DO
    CLOSE (7)
  END DO
!
  DO KK=1,MNRE
    KR=IREV(KK)
    A=SUM(NEVIB(KK,1,:))
    IF(A == 0) A=1.d0
    WRITE(FILENAME,783) KK
    783 FORMAT('DS1VIBE1_',i2.2)
    OPEN (7,FILE=TRIM(FILENAME)//'.DAT')
    WRITE (7,994) '# Sampled pre-reaction total vibrational energiy'
    WRITE (7,994) '# Depleting and replenishing reactions, IKA-IREV(IKA)-Total:',KK,KR,A,REAC(KK)
    WRITE (7,994) '# Pre-reaction species:',IREA(1,KK),IREA(2,KK),JREA(2,KK,1)
    WRITE (7,994) '# Temperatures:',VAR(8:11,NSPDF)
    WRITE (7,994) '# [Ev/(BOLTZ*1e3)], Sample and Fraction'
    DO I=0,100
      WRITE (7,993) I,NEVIB(KK,1,I),NEVIB(KK,1,I)/A  !pre-reaction total vibrational energy
    END DO
    CLOSE (7)
  END DO
!
  IF (IREAC == 2) THEN
    FEVIB=0.d0
    DO KK=1,MNRE
      A=SUM(NEVIB(KK,1,:))
      IF(A == 0) A=1.d0
      DO I=0,100
        FEVIB(KK,1,1,I)=NEVIB(KK,1,I)/A  !pre-reaction total vibrational energies
      END DO
    END DO
!
    J=1
    FPVIB=0.d0
    DO KK=1,MNRE
      A=SUM(NPVIB(J,KK,1,1,:)) !molecule 1
      B=SUM(NPVIB(J,KK,2,1,:)) !molecule 2
      C=SUM(NPVIB(J,KK,3,1,:)) !molecule 3
      IF(A == 0) A=1.d0
      IF(B == 0) B=1.d0
      IF(C == 0) C=1.d0
      DO I=0,100
        FPVIB(KK,1,1,1:3,I)=NPVIB(J,KK,1,1:3,I)/A  !pre-reaction vibrational levels pdfs
        FPVIB(KK,1,2,1:3,I)=NPVIB(J,KK,2,1:3,I)/B
        FPVIB(KK,1,3,1:3,I)=NPVIB(J,KK,3,1:3,I)/C
      END DO
      ZCHECK=1234567
      WRITE(FILENAME,786) KK,J
      786 FORMAT('DS1VIBF1_',i2.2,'_',i2.2)
      787 CONTINUE
      OPEN (7,FILE=TRIM(FILENAME)//'.BIN',FORM='UNFORMATTED',ERR=787)
      WRITE (7) VAR(11,NSPDF),FEVIB(KK,J,:,:),FPVIB(KK,J,:,:,:),ZCHECK
      CLOSE (7)
    END DO
  END IF
END IF
!
!----------------------------------------------------------------------------
!--sample molecule properties for SNAPSHOT files
!
NSS=0
DO K=1,NSNAP
  CALL ZGF(RANF,IDT)
  N=INT(RANF*DFLOAT(NM))+1
  NSS=NSS+1
  ISNAP(1,NSS)=IPSP(N)
  PSNAP(1,NSS)=PX(N)
  VSNAP(:,NSS)=PV(:,N)
END DO
!
!----------------------------------------------------------------------------
!--write binary files for post processing of multiple run cases and snapshot

ZCHECK=1234567
!
102 CONTINUE
  !OPEN (7,FILE='DS1OUT_'//E//'.BIN',FORM='UNFORMATTED',ERR=102)
!trick  !OPEN (7,FILE='DS1OUT.BIN',FORM='UNFORMATTED',ERR=102)
  !WRITE (7)IFX,NOUT,FTIME,NSAMP,TISAMP,NM,TOTMOV,TOTCOL,PCOLLS,TOTDUP,TCOL,CSCR,&
  !         CTIME,TOTCOLI,TOTMOVI,NDISSOC,NRECOMB,&
  !         ITYPE,XB,VARS,NCELLS,VARSP,&  !end of general properties
  !         JCD,TREACG,TREACL,ITCV,IEAA,IZV,REAC,VAR,DTM,CELL,NCIS,CST,UVFX,NMS,& !end of reactions,flow properties, and composition
  !         IPDF,IPRS,BIN,BINS,DBINV,SPI,PDF,PDFS,NSPDF,NDROT,NDVIB,DBINC,DBINE,ISNAP,PSNAP,VSNAP,SP,ZCHECK    !end of sample pdf and snapshot
!trick  !CLOSE(7)
!
WRITE (*,*) 'Output and binary files are written'
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
992 FORMAT (A,I2.2)
!
!----------------------------------------------------------------------------
!--reset collision and transit times etc.
!
DTMI=DTM
DTM=DTM*2.
!--this makes it possible for DTM to increase, it will be reduced as necessary
DO N=1,NCCELLS
!
  NN=ICCELL(3,N)
  B=(CELL(3,NN)-CELL(2,NN))/DFLOAT(NCIS)     !--collision cell width
!
  IF ((NN > 0).AND.(NN <= NCELLS)) THEN
!
    IF (VAR(13,NN) > 20.D00) THEN
!--consider the local collision rate
      CCELL(3,N)=0.5D00*VAR(14,NN)*CPDTM
!--look also at collision cell transit time based on the local flow speed
      A=0.5D00*(B/(ABS(VAR(5,NN))))*TPDTM
      IF (A < CCELL(3,N)) CCELL(3,N)=A
      IF (2.D00*CCELL(3,N) < DTM) DTM=2.D00*CCELL(3,N)
    ELSE
!-- base the time step on a collision cell transit time at the refence vmp
      A=TPDTM*B/VMPM
      IF (A < CCELL(3,N)) CCELL(3,N)=A
    END IF
    IF (1.9D00*CCELL(3,N) < DTM) DTM=1.9D00*CCELL(3,N)
!
  END IF
END DO
!
WRITE (9,*) 'DTM changes  from',DTMI,' to',DTM
WRITE (*,*) 'DTM changes  from',DTMI,' to',DTM
!
!----------------------------------------------------------------------------
!--update output interval
!
TPOUT=OUTRAT
IF (ISF > 0) THEN
  IF ((NOUT >= 0  ).AND.(NOUT < 100)) TPOUT=OUTRAT*.1d0
  IF ((NOUT >= 0  ).AND.(NOUT < 100)) TPOUT=OUTRAT*.1d0
  IF ((NOUT >= 100).AND.(NOUT < 150)) TPOUT=OUTRAT*.2d0
!  IF ((NOUT >= 150).AND.(NOUT < 170)) TPOUT=OUTRAT*.5d0
!  IF ((NOUT >= 200)) TPOUT=OUTRAT*INT(.9999999+(NOUT-190)/10) !comment for RELAX.DAT
!
  IF (NOUT >= 50 ) CPDTM=0.005
  IF (NOUT >= 100) CPDTM=0.010
  IF (NOUT >= 150) CPDTM=0.050
  IF (NOUT >= 200) CPDTM=0.100
!  IF ((CPDTM < 0.02  ).AND.(NOUT >= 200)) CPDTM=0.02
!  IF ((CPDTM < 0.2   ).AND.(NOUT >= 250)) CPDTM=0.2
END IF
!
IF (ISF==0) THEN !special coding for shockwave sampling
  TPOUT=OUTRAT*2.d0
END IF
!
IF (MNRE > 0) THEN
  IF ((IRM == 201).AND.(FTIME > 4.d-5)) ISF=0 !steady state is achieved
  IF ((IRM == 202).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 206).AND.(FTIME > 4.d-5)) ISF=0
  IF ((IRM == 213).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 214).AND.(FTIME > 8.d-5)) ISF=0
  IF ((IRM == 218).AND.(FTIME > 4.d-5)) ISF=0
  IF ((IRM == 227).AND.(FTIME > 2.d-6)) ISF=0
  IF ((IRM == 229).AND.(FTIME > 2.d-6)) ISF=0
  IF ((IRM == 231).AND.(FTIME > 1.d-1)) ISF=0
END IF
IF ((IGS == 3).AND.(FTIME > 3.d-3)) ISF=0
!
DTSAMP=DTSAMP*DTM/DTMI
DTOUT=DTSAMP*TPOUT
TOUT=FTIME
!
WRITE (9,*) 'NOUT:',NOUT,' OUTRAT:',TPOUT
WRITE (*,*) 'NOUT:',NOUT,' OUTRAT:',TPOUT
!
!----------------------------------------------------------------------------
!
RETURN
!
END SUBROUTINE OUTPUT_RESULTS
!
!**************************************************************************************
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
SUBROUTINE DISSOCIATION
!
!--dissociate molecules that have been marked with IPVIB(0,:) < 0
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
!
IMPLICIT NONE
!
INTEGER ::K,KK,L,M,LS,MS,KS,JS,KV,IKA,NC,NS,NSP,MAXLEV,II,IV,IDT=0,ITEST=0 !--isebasti: included IDT
REAL(KIND=8) :: A,B,C,ECT,ECC,EVIB,VRR,VR,RMM,RML,RANF,SVDOF(2),PROB,ERM !--isebasti: included RANF
REAL(KIND=8), DIMENSION(3) :: VRC,VCM,VRCP
!
L=0
DO WHILE (L < NM)
  L=L+1
  LS=IPSP(L)
  IF (ISPV(LS) > 0) THEN
    IF (IPVIB(0,L) < 0) THEN  !dissociate through reaction IKA
      IKA=-IPVIB(0,L)
      ECT=PROT(L)  !--ECT is the energy available for post-reaction states
      DO KV=1,ISPV(LS)
        CALL VIB_ENERGY(EVIB,IPVIB(KV,L),KV,LS)
        ECT=ECT+EVIB !IPVIB(KV,L)*BOLTZ*SPVM(1,KV,LS)
      END DO
      IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
      NM=NM+1
      M=NM
!--update species
      IPSP(L)=JREA(1,IKA,1)
      IPSP(M)=JREA(1,IKA,2)
      LS=IPSP(L)
      MS=IPSP(M)
!--set center of mass velocity as that of molecule
      VCM(1:3)=PV(1:3,L)
!--new born molecule
      PX(M)=PX(L)
      IPCELL(M)=IPCELL(L)
      PTIM(M)=PTIM(L)
!--set any internal mode to the ground state and reset IPVIB(0,L) and IPVIB(0,M) to 0
      PROT(L)=0.d0
      PROT(M)=0.d0
      IPVIB(:,L)=0
      IPVIB(:,M)=0
!--effective vibrational degrees of freedom
      NC=IPCELL(L)        !collision cell
      NS=ICCELL(3,NC)     !sampling cell
      A=VAR(10,NS)        !sampling cell vibrational temperature
      SVDOF(:)=0.d0
      IF (ISPV(LS) > 0) THEN
        DO KV=1,ISPV(LS)
          SVDOF(1)=SVDOF(1)+2.d0 !*(SPVM(1,KV,LS)/A)/(DEXP(SPVM(1,KV,LS)/A)-1.d0)
        END DO
      END IF
      IF (ISPV(MS) > 0) THEN
        DO KV=1,ISPV(MS)
          SVDOF(2)=SVDOF(2)+2.d0 !*(SPVM(1,KV,MS)/A)/(DEXP(SPVM(1,KV,MS)/A)-1.d0)
        END DO
      END IF
IF (ITEST == 1) THEN
!--Larsen-Borgnakke serial vibrational energy redistribution
      DO NSP=1,2
        IF (NSP == 1) THEN
          K=L ; KS=LS ; JS=MS
        ELSE
          K=M ; KS=MS ; JS=LS
        END IF
        IF (ISPV(KS) > 0) THEN
          DO KV=1,ISPV(KS)
            CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
            !EVIB=DFLOAT(IPVIB(KV,K))*BOLTZ*SPVM(1,KV,KS)
            ECC=ECT+EVIB
            CALL VIB_LEVEL(ECC,MAXLEV,KV,KS)
            !MAXLEV=ECC/(BOLTZ*SPVM(1,KV,KS))
            CALL ZGF(RANF,IDT)
            II=0
            DO WHILE (II == 0)
              CALL ZGF(RANF,IDT)
              IV=RANF*(MAXLEV+0.99999999D00)
              IPVIB(KV,K)=IV
              CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
              !EVIB=DFLOAT(IV)*BOLTZ*SPVM(1,KV,KS)
              IF (EVIB < ECC) THEN
                IF (SVDOF(NSP) > 0.d0) SVDOF(NSP)=SVDOF(NSP)-2.d0 !*(SPVM(1,KV,KS)/VAR(10,NN))/(DEXP(SPVM(1,KV,KS)/VAR(10,NN))-1.d0)
                IF (NSP == 1) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,KS)+ISPR(1,JS)+SVDOF(2)+SVDOF(1)
                IF (NSP == 2) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,KS)+ISPR(1,JS)+SVDOF(2)
                PROB=(1.D00-EVIB/ECC)**(A*0.5d0-1.d0)   !--see eqn (5.61) and (5.45)
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) II=1
              END IF
            END DO
            ECT=ECC-EVIB
          END DO
        END IF
      END DO
!-------------------------------------------------
!--Larsen-Borgnakke serial rotational energy redistribution
      DO NSP=1,2
        IF (NSP == 1) THEN
          K=L ; KS=LS ; JS=MS
        ELSE
          K=M ; KS=MS ; JS=LS
        END IF
        IF (ISPR(1,KS) > 0) THEN
          ECC=ECT+PROT(K)
          IF (NSP == 1) A=(5.d0-2.d0*SPM(3,KS,JS))+ISPR(1,JS)
          IF (NSP == 2) A=(5.d0-2.d0*SPM(3,KS,JS))
          CALL LBS(ISPR(1,KS)*0.5d0-1.d0,A*0.5d0-1.d0,ERM,IDT)
          PROT(K)=ERM*ECC
          ECT=ECC-PROT(K)
        END IF
      END DO
END IF
!-------------------------------------------------
!--calculate new velocities
      VRR=2.d0*ECT/SPM(1,LS,MS)
      VR=DSQRT(VRR)
      VRC=0.d0                  !initial relative velocity
      RML=SPM(1,LS,MS)/SP(5,MS)
      RMM=SPM(1,LS,MS)/SP(5,LS)
      CALL VHSS(LS,MS,VR,VRC,VRCP,IDT)
      PV(1:3,L)=VCM(1:3)+RMM*VRCP(1:3)
      PV(1:3,M)=VCM(1:3)-RML*VRCP(1:3)
      IPCP(L)=M
      IPCP(M)=L
    END IF
  END IF
END DO
!
RETURN
!
END SUBROUTINE DISSOCIATION
!
!************************************************************************************
!
SUBROUTINE CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,VRR,ECT,SXSECTION,IVDC,IDT)
!
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
USE OMP_LIB
USE MFDSMC,only : IMF, IMFdia, NMFANG,&
  MF_SAMPLE_PHASE,MF_CALC_COLL,MF_SAMPLE_ANGLE,MF_EVAL_F
!
IMPLICIT NONE
!
!
REAL(8),INTENT(IN) :: VRR,ECT,SXSECTION
INTEGER,INTENT(IN) :: N,L,M,LS,MS

INTEGER :: J,K,NRE,KK,KS,MK,KA,KAA,MKK,KV,IA,ISTE(MNRE),&
           IDT,NS,NPM,I,IVDC(MNRE),II,MTYPE,IVIB,JROT,IMAX,JMAX
REAL(8),EXTERNAL :: GAM
REAL(KIND=8) :: A,B,ECR,ECV,EC,THBCELL,RANF,&
                VDOF1(MMVM),VDOF2(MMVM),EV1(MMVM),EV2(MMVM),EV,ECV1,ECV2,SVDOF1,SVDOF2,ECR1,ECR2,AL,&
                ATDOF,AIDOF,ANDOF,AVDOF,X,C,D,AA,BB,XI,PSI,PHI,CV,EE,SEV(2),TEMP,&
                RXSECTION(MNRE),STER(MNRE),ECA(MNRE),CF,&
                ALPHA1,ALPHA2,EA,CC,DD
REAL(KIND=8),ALLOCATABLE :: VEC_I(:),VEC_J(:),VALUES(:),ARRAY_TEMP(:,:)
REAL(KIND=8) :: MFANG(8),MFV1,MFV2,MFR1,MFR2,MFCtheta
REAL(KIND=8) :: MFF(2),MFcoll(2)
INTEGER      :: NMFCALL, viblevel(2), IINPM
LOGICAL      :: ISAME
!
!--A,B,C working variables
!--J,K,KK,MK,KA,KAA,MKK,JR,KR,JS,KV,IA working integers
!--N the collision cell
!--L,S  the molecule numbers
!--LS,MS,KS species
!--ECR rotational energy
!--ECV vibrational energy
!--ECA available energy
!--XSECTION cross-section
!--ISTE third body species
!--NRE number of reactions that may occur for this pair
!--TXSECTION total cross-section
!--THBCELL number of third body molecules in cell
!--WF weighting factor
!--VCM center-of-mass velocity components
!--VRR square of rel speed
!--RML,RMM molecule mass parameters
!--MFANG angles sampled by MF-MC model
!--MFF threshold function maximum 4 kinds of configuration
!--SXSECTION total collision cross sections
!
NMFCALL = 0           !count of MFmodel using
NS=ICCELL(3,N)        !sampling cell
TEMP=VAR(10,NS)       !sampling cell vibrational temperature
NRE=NRSP(LS,MS)       !number of reaction involving species LS and MS
ISAME = .false.             !check if two species are the same
IF (LS == MS) ISAME = .true.
!
IF (NRE > 0) THEN
!
  MFF = ECT*2.0d0  !initialization
  ECR=0.D00; ECV=0.D00
  EV1=0.D00; EV2=0.D00
  ECR1=0.D00; ECR2=0.D00
  ECV1=0.D00; ECV2=0.D00
  VDOF1=0.D00; VDOF2=0.D00
  SVDOF1=0.D00; SVDOF2=0.D00
  IF (ISPR(1,LS) > 0) ECR1=PROT(L)
  IF (ISPR(1,MS) > 0) ECR2=PROT(M)
  ECR=ECR1+ECR2
  IF (MMVM > 0) THEN
    IF (ISPV(LS) > 0) THEN
      DO KV=1,ISPV(LS)
        CALL VIB_ENERGY(EV1(KV),IPVIB(KV,L),KV,LS)
        VDOF1(KV)=2.d0*(SPVM(1,KV,LS)/TEMP)/(DEXP(SPVM(1,KV,LS)/TEMP)-1.d0)
        ECV1=ECV1+EV1(KV)
        SVDOF1=SVDOF1+VDOF1(KV)
      END DO
    END IF
    IF (ISPV(MS) > 0) THEN
      DO KV=1,ISPV(MS)
        CALL VIB_ENERGY(EV2(KV),IPVIB(KV,M),KV,MS)
        VDOF2(KV)=2.d0*(SPVM(1,KV,MS)/TEMP)/(DEXP(SPVM(1,KV,MS)/TEMP)-1.d0)
        ECV2=ECV2+EV2(KV)
        SVDOF2=SVDOF2+VDOF2(KV)
      END DO
    END IF
    ECV=ECV1+ECV2
  END IF
!
  EC=ECT+ECR+ECV  !total collision energy
  RXSECTION=0.d0  !initializing reaction cross-sections
  STER=0.d0       !initializing steric factors used in TCE model
!
  AL=SPM(3,LS,MS)-0.5d0
  ATDOF=(ISPR(1,LS)+ISPR(1,MS)+SVDOF1+SVDOF2+(4.d0-2.d0*AL))*0.5d0 !average number of dofs
  AIDOF=(ISPR(1,LS)+ISPR(1,MS)+SVDOF1+SVDOF2)*0.5d0                !average number of internal dofs
!
  DO J=1,NRE
    K=IRCD(J,LS,MS)
    NPM=NREA(1,K)+NREA(2,K)  !number of post-reaction species
    X=REA(4,K)-0.5d0+AL
    ECA(J)=EC
    CF=1.d0
!
!--collision energy is above reaction energy barrier
    IF ((EC >= REA(2,K)).AND.(EC+REA(5,K) > 0)) THEN !2nd condition prevents a negative post-reaction total energy
!
!--check if molecules follow the order given by IREA
      IVDC(K) = 1
      IF (ISAME) THEN
        IF (ISPV(LS) > 0)THEN
         ! choose the molecule with higher dissociation energy
         ! some model like MFDSMC may change this
          IF (ECV1 .ge. ECV2) THEN
            IVDC(K) = 1
          ELSE
            IVDC(K) = 2
          END IF
        END IF
      ELSE
        IF (LS == IREA(1,K)) IVDC(K) = 1  !corrected order
        IF (MS == IREA(1,K)) IVDC(K) = 2  !inverted order
      END IF
!
!--check reaction model
      MTYPE=0 !TCE is the standard model
      IF (QCTMODEL==2.AND.GASCODE==8)THEN !QCT model
        IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) MTYPE=1 !N2-O collision
        IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) MTYPE=3 !O2-O collision
      END IF
!     MTYPE=0 !to enforce the use of TCE model for all collisions
!
!--------------------------------------------------------------------------------
      IF (MTYPE == 0) THEN !use TCE/MFDSMC/VDC model
!
!--possible recombination reaction
        IF (NPM == 1) THEN
          MKK=0
          CALL ZGF(RANF,IDT)
          IF (IRECOM == 1) THEN! .AND. RANF < 0.001d0) THEN  !only one in 1000 possible recombinations are considered and
                                                             !the probability is increased by a factor of 1000
            IF (JREA(2,K,1) < 0) THEN
!--third-body molecule may be any species in cell
              IF (ICCELL(2,N) >= 3) THEN
                MK=L
                DO WHILE ((MK == L).OR.(MK == M))
                  CALL ZGF(RANF,IDT)
                  KK=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                  MK=ICREF(KK)
                END DO
                KS=IPSP(MK)
                ISTE(J)=MK
                A=REA(6,K)
                IA=-JREA(2,K,1)
                A=A*THBP(IA,KS)
              ELSE
                A=0.D00
              END IF
              THBCELL=ICCELL(2,N)
            ELSE
              KS=JREA(2,K,1)
!--the third-body molecule must be species KS
              THBCELL=0.
              MKK=0
              DO KAA=1,ICCELL(2,N)
                KA=ICCELL(1,N)+KAA
                MK=ICREF(KA)
                IF ((IPSP(MK) == KS).AND.(IPCELL(MK) > 0).AND.(IPVIB(0,MK) >= 0)) THEN
                  THBCELL=THBCELL+1.
                  IF ((MK /= L).AND.(MK /= M)) MKK=MK
                END IF
              END DO
              IF (MKK > 0) THEN
                MK=MKK
                ISTE(J)=MK
                A=REA(6,K)
              ELSE
                A=0.
              END IF
            END IF
            B=THBCELL*FNUM/CCELL(1,N)  !number density of third body species in collision cell
            CC=3.5d0-AL !3.5d0-AIDOF-AL
            DD=3.5d0-AL+X !3.5d0-AIDOF-AL+X
            C=GAM(CC)
            D=GAM(DD)
            STER(J)=CF*B*A*(C/D)*ECA(J)**X  !*1000.d0  !note the factor
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          END IF
        END IF
!
!--possible exchange reaction
        IF (NPM == 2) THEN
          CC=ATDOF
          DD=AIDOF+1.5d0+REA(4,K)
          C=GAM(CC)
          D=GAM(DD)
          AA=X+ATDOF-1.d0
          BB=ATDOF-1.d0
          !--note measures to prevent underflows
          STER(J)=CF*(REA(6,K)*(C/D))*(((ECA(J)-REA(2,K))*1.d10)**AA/(ECA(J)*1.d10)**BB)*1.d10**(BB-AA)
          RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
        END IF
!
!--possible dissociation reaction
        IF (NPM == 3) THEN
          IF(REA(1,K) < 0.) THEN
            IF (IMF == 0 .or. (IMF .ne. 0 .and. NMFANG(K) == 0))THEN
              !TCE model (same as in exchange reactions)
              CC=ATDOF
              DD=AIDOF+1.5d0+REA(4,K)
              C=GAM(CC)
              D=GAM(DD)
              AA=X+ATDOF-1.d0
              BB=ATDOF-1.d0
              !--testing nonequilibrium reaction rate corrections
              IF (IVDC(K) == 1) KS=LS
              IF (IVDC(K) == 2) KS=MS
              !
              !--Kuznetsov (KSS)
              !E=(1.d0-DEXP(-SPVM(1,1,KS)/VAR(10,NS)))/(1.d0-DEXP(-SPVM(1,1,KS)/VAR(8,NS)))
              !A=1.d0/(BOLTZ*VAR(10,NS))
              !B=1.d0/(BOLTZ*VAR(8,NS))
              !CF=E*DEXP(-0.7d0*REA(2,K)*(A-B)) !note that the tecplot exact curves are considering a 0.78 factor
              !
              !--testing modified (0.78) Macharet-Fridam correction
              !IF (LS==3.AND.MS==3) THEN !TCE+MF only for O2-O2 collisions
              !  E=(1.d0-DEXP(-SPVM(1,1,KS)/VAR(10,NS)))/(1.d0-DEXP(-SPVM(1,1,KS)/VAR(8,NS)))
              !  A=1.d0/(BOLTZ*VAR(10,NS))
              !  B=1.d0/(BOLTZ*VAR(8,NS))
              !  ALP=(1.d0/2.d0)**2.d0                  !alpha; for homonuclear collisions
              !  TA=ALP*VAR(10,NS)+(1-ALP)*VAR(8,NS)    !T_a
              !  IF (ISPV(LS) >= 1 .AND. ISPV(MS) >= 1) THEN
              !    LF=2.d0*(1-ALP)/(PI**2.d0*ALP**0.75d0) * (BOLTZ*VAR(8,NS)/REA(2,K))**(1.5d0-REA(4,K)) &
              !                                           * (1+3.5*(1-ALP)*(1+DSQRT(ALP))*(BOLTZ*VAR(8,NS)/REA(2,K)))
              !  ELSE
              !    LF=(9.d0/64.d0)*DSQRT(PI*(1-ALP)) * (BOLTZ*VAR(8,NS)/REA(2,K))**(1.d0-REA(4,K)) &
              !                                      * (1+2.5*(1-ALP)*(BOLTZ*VAR(8,NS)/REA(2,K)))
              !  END IF
              !  PHI2=1.373274d6*VAR(8,NS)**(-1.309734d0)
              !  CF=E*(1-LF)*DEXP(-0.78d0*REA(2,K)*(A-B)) + PHI2*LF*DEXP(-0.78d0*REA(2,K)*(1.d0/(BOLTZ*TA)-B))
              !  IF (CF > 1.d0) CF=1.d0 !to correct only vibrationally cold conditions
              !END IF
              !
              !--note measures to prevent underflows
              STER(J)=CF*(REA(6,K)*(C/D))*(((ECA(J)-REA(2,K))*1.d10)**AA/(ECA(J)*1.d10)**BB)*1.d10**(BB-AA)
            ELSE
              ! Macheret-Fridman-MC mode
              !
              !---get energy of each one particle
              ! here I assume that if the molecule doesn't have vibrational energy
              ! IPVIB will be 0
              IF (IVDC(K) == 1)THEN
                MFV1 = ECV1; MFR1 = ECR1
                MFV2 = ECV2; MFR2 = ECR2
                viblevel(1) = IPVIB(1,L); viblevel(2) = IPVIB(1,M)
              ELSE
                MFV1 = ECV2; MFR1 = ECR2
                MFV2 = ECV1; MFR2 = ECR1
                viblevel(1) = IPVIB(1,M); viblevel(2) = IPVIB(1,L)
              END IF
              !
              NMFCALL = NMFCALL + 1
              IF (NMFCALL > 2) THEN
                WRITE (*,"(A,I3,'+',I3)") "There are more than 2 MF-diss for ",LS,MS
                STOP
              END IF
              !---set mass of mb (MFcol(1)) and mc (MFcol(2))
              CALL MF_CALC_COLL(K,NMFCALL,MFcoll,IDT)
              !
              !---generate geometries of collision
              !
              ! atom-diatom: gamma1, gamma2, theta, phi
              ! diatom-diatom: gamma1, gamma2, theta, phi0, phi1, beta1, beta2, delta, beta
              ! MFANG:
              !  1. gamma_1
              !  2. gamma_2
              !  3. gamma_
              !  4. phi_0
              !  5. phi_1
              !  6. beta_1
              !  7. beta_2
              !  8. delta
              CALL MF_SAMPLE_ANGLE(K,NMFCALL,MFANG,MFCtheta,viblevel,MFV1,MFV2,IDT)
                ! log angles for binning
                !IF (IREAC == 2)THEN
                  !DO II = 1,NMFANG(K)
                    !JJ = FLOOR(MFANG(II)*180.0d0)+1
                    !JJ = MIN(JJ,180)
                    !!              NCANGLE(II,JJ) = NCANGLE(II,JJ)+1.d0
                  !END DO
                !END IF
              CALL MF_EVAL_F(K,NMFCALL,MFV1,MFV2,MFR1,MFR2,MFcoll,MFANG,MFCtheta,MFF(NMFCALL),IDT)
              STER(J) = 0.0d0
              IF (ECT >= MFF(NMFCALL)) STER(J) = 1.0D0

              IF (ISAME .and. IMFdia == 1) THEN
                ! check the possibility of dissociate the molecule with lower vibrational energy
                NMFCALL = NMFCALL + 1
                ! exchange of viblevel
                ! this is actually not needed, but we have it here to keep consistency
                II = viblevel(1)
                viblevel(1) = viblevel(2)
                viblevel(2) = II
                ! exchange MFcoll
                CALL MF_CALC_COLL(K,NMFCALL,MFcoll,IDT)
                ! exchange some angles
                CALL MF_SAMPLE_ANGLE(K,NMFCALL,MFANG,MFCtheta,viblevel,MFV2,MFV1,IDT)
                ! exchange internal energy
                CALL MF_EVAL_F(K,NMFCALL,MFV2,MFV1,MFR2,MFR1,MFcoll,MFANG,MFCtheta,MFF(NMFCALL),IDT)

                IF (MFF(2) < MFF(1) .and. MFF(2) <= ECT) THEN
                  STER(J) = 1.0d00
                  IVDC(K) = 3-IVDC(K)
                  ! the configuration IVDC(K) is prefered
                END IF
              ELSE
                IF ((NMFCALL == 2) .and. (STER(J) > 0.9d0) )THEN
                  IF (MFF(2) < MFF(1))THEN
                    DO II=1,J-1
                      IINPM = NREA(1,IRCD(II,LS,MS))+NREA(2,IRCD(II,LS,MS))
                      IF (IINPM == 3 .and. NMFANG(II) > 0) then
                        STER(II) =0.0d0 ! set previous dissociation to 0
                        RXSECTION(II) = 0.0d0
                      END IF
                    END DO
                  ELSE
                    STER(J) = 0.0d0
                  END IF
                END IF
              !ELSE  ! MF probability version
                !STER(J)=0.0D0
                !IF (MFV1 .LE. 1.0D-5) MFV1 = 0.5d0*BOLTZ*SPVM(1,1,IREA(1,K))
                !C = MFALPHA(K)*MFDSTAR
                !CC = DSQRT(MFV1 / MFDSTAR)
                !DD = DSQRT(MFALPHA(K))
                !IF (MFV1 .LT. C) THEN
                  !MFF = MFDSTAR/(1.0d0-MFALPHA(K))*(1.0d0-DD*CC)**2
                  !IF ( ECT > MFF) THEN
                    !STER(J) = 4.0D0*(ECT-MFF)**1.5d0 /(3.0d0*PI*PI*MFF*DSQRT(MFDSTAR/(1.0d0-MFALPHA(K))))
                    !STER(J) = STER(J) / DSQRT((1.0d0-CC*DD)*(1.0d0-(2.0d0-DD)*CC))
                  !END IF
                !ELSE
                  !MFF = MFDSTAR - MFV1
                  !IF ( ECT > MFF) THEN
                    !STER(J) = (ECT-MFF)**2/(2.0d0*PI*PI*MFF*MFDSTAR*DD/DSQRT((1.0d0+DD)/(1.0d0-DD)))
                    !STER(J) = STER(J) / DSQRT((1.0d0+MFV1/MFDSTAR)*(MFV1/MFDSTAR/MFALPHA(K)-1.0D0))
                  !END IF
                !END IF
                !IF ( DABS(MFV1-MFDSTAR) < 1.0D0 ) STER(J) = 1.0D0
                !IF (STER(J) > 0.999D0)  STER(J) = 1.0D0
              !END IF
              END IF
            END IF
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          ELSE
            !VDC model
            IF (IVDC(K) == 1) AVDOF=SVDOF1*0.5d0    !internal dofs contributing to vibrational favoring
            IF (IVDC(K) == 2) AVDOF=SVDOF2*0.5d0
            ANDOF=ATDOF-AVDOF                    !internal dofs not contributing to vibrational favoring
            XI=ANDOF-1.d0
            PSI=REA(4,K)+ATDOF+0.5d0-(4.d0-2.d0*AL)*0.5d0
            PHI=REA(1,K)
            CC=ANDOF
            DD=PSI+1.d0
            C=GAM(CC)
            D=GAM(DD)
            CV=1.d0
            BB=1.d0
            EE=1.d0
            EV=0.d0
            IF (MMVM > 0) THEN
              IF (IVDC(K) == 1) THEN
                DO KV=1,ISPV(LS)
                  AA=(TEMP/SPVM(1,KV,LS))**(VDOF1(KV)*0.5d0)
                  SEV=0.d0
                  DO I=1,IDL(KV,LS)
                    SEV(1)=SEV(1)+(DFLOAT(I)*SPVM(1,KV,LS))**PHI
                    SEV(2)=SEV(2)+DEXP(-DFLOAT(I)*SPVM(1,KV,LS)/TEMP)
                  END DO
                  CV=CV*AA*SEV(1)/SEV(2)
                  BB=BB*SPVM(1,KV,LS)**(VDOF1(KV)*0.5d0)
                  EE=EE*(EV1(KV)/BOLTZ)**PHI
                END DO
                EV=ECV1
              ELSE
                DO KV=1,ISPV(MS)
                  AA=(TEMP/SPVM(1,KV,MS))**(VDOF2(KV)*0.5d0)
                  SEV=0.d0
                  DO I=1,IDL(KV,MS)
                    SEV(1)=SEV(1)+(DFLOAT(I)*SPVM(1,KV,MS))**PHI
                    SEV(2)=SEV(2)+DEXP(-DFLOAT(I)*SPVM(1,KV,MS)/TEMP)
                  END DO
                  CV=CV*AA*SEV(1)/SEV(2)
                  BB=BB*SPVM(1,KV,MS)**(VDOF2(KV)*0.5d0)
                  EE=EE*(EV2(KV)/BOLTZ)**PHI
                END DO
                EV=ECV2
              END IF
            END IF
            STER(J)=CF*(REA(6,K)*(C/D)/(CV*BB))*(((ECA(J)-REA(2,K))/BOLTZ)**PSI/((ECA(J)-EV)/BOLTZ)**XI)*EE
            RXSECTION(J)=STER(J)*SXSECTION  !in Angstrons^2
          END IF
        END IF
!
      ELSE
!--------------------------------------------------------------------------------
!
!--possible exchange reaction with QCT model
        IF (NPM == 2) THEN
          IF (MTYPE==1) THEN
            ! data only for N2+O=>NO+N reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 50) IVIB=50
!
            IMAX=60
            ALLOCATE (VEC_I(0:IMAX))
            VEC_I=(/2.454482d-4,2.427769d-4,2.401402d-4,2.375291d-4,2.348600d-4,2.322024d-4,&
                    2.296469d-4,2.270147d-4,2.243892d-4,2.217688d-4,2.191515d-4,2.166443d-4,&
                    2.140356d-4,2.114257d-4,2.088131d-4,2.061956d-4,2.036905d-4,2.010655d-4,&
                    1.984311d-4,1.957851d-4,1.932510d-4,1.905840d-4,1.878997d-4,1.851937d-4,&
                    1.826079d-4,1.798666d-4,1.772361d-4,1.745758d-4,1.717679d-4,1.689237d-4,&
                    1.661840d-4,1.632698d-4,1.604614d-4,1.576103d-4,1.547141d-4,1.516307d-4,&
                    1.486383d-4,1.455892d-4,1.424788d-4,1.393016d-4,1.360511d-4,1.327185d-4,&
                    1.294563d-4,1.259450d-4,1.224778d-4,1.187449d-4,1.150340d-4,1.111769d-4,&
                    1.071570d-4,1.029527d-4,9.853521d-5,9.402503d-5,8.922437d-5,8.408293d-5,&
                    7.852777d-5,7.259131d-5,6.598628d-5,5.858944d-5,4.994199d-5,3.901624d-5,&
                    1.683039d-5/)
            B=VEC_I(IVIB)  !ro-vibrational energy ladder alpha_1 fitting parameter for N2
            JROT=-0.5d0+DSQRT(0.25d0+(ECR/EVOLT)/B)
            B=2.88d0    !characteristic rot temperature for N2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 0) JROT=0
            IF (JROT > 150) JROT=150
            DEALLOCATE (VEC_I)
!
            IMAX=5
            JMAX=10
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,20,50,100,150/)         !rot levels
            VEC_J=(/0,1,5,7,10,15,20,25,30,50/) !vib levels
!
            VALUES=(/ 106.6040, 104.3597, 108.0698, 351.2132,  0.    ,&
                       87.8864,  70.9021, 197.9063, 359.9827, 46.6272,&
                      266.9474, 259.8274, 206.4693, 463.3182, 53.9792,&
                      417.9959, 385.7538, 305.3864, 409.3383, 48.7572,&
                      800.0000, 672.3581, 240.2589, 196.8796, 42.5553,&
                      800.0000, 622.1744, 260.5081, 160.8333, 37.3285,&
                      321.1994, 269.5515,  58.4387,  84.0821, 36.1900,&
                      213.5262, 100.3110,  44.1216,  73.1788, 21.1884,&
                      107.6799,  76.4448,  93.3640,  88.2306, 17.3448,&
                      80.4866,   94.9898,  71.1551,   0.    ,  0.     /)
            VALUES=VALUES/3.d0
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),A)
!
            VALUES=(/-2.1495, -2.1276, -2.0766, -2.5219, -2.3018,&
                     -2.0134, -1.8961, -2.3555, -2.5117, -2.3504,&
                     -2.4319, -2.4155, -2.2638, -2.7787, -2.5659,&
                     -2.6048, -2.5582, -2.4071, -2.8969, -2.5814,&
                     -2.8641, -2.7715, -2.2370, -2.6994, -2.6589,&
                     -3.1547, -3.0446, -2.6986, -2.9661, -2.8068,&
                     -3.0753, -2.9859, -1.9479, -2.7587, -3.1888,&
                     -3.2309, -2.5819, -1.9158, -3.0487, -2.9818,&
                     -3.0438, -2.6741, -3.1186, -3.9082, -2.9116,&
                     -6.6570, -7.7534, -9.3717,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA1)
!
            VALUES=(/ 3.2628, 3.2630, 3.2628, 3.2629, 6.4555,&
                      3.2628, 3.2637, 3.2628, 3.2632, 6.5782,&
                      3.2629, 3.2628, 3.2629, 3.7551, 6.9618,&
                      3.2628, 3.2635, 3.2629, 4.2386, 7.2907,&
                      3.2628, 3.2629, 3.3366, 4.9327, 7.8277,&
                      3.9871, 4.0736, 4.5092, 6.0078, 8.5064,&
                      5.0910, 5.1721, 5.5801, 6.9790, 9.1292,&
                      6.0918, 6.1672, 6.5467, 7.8421, 9.9821,&
                      6.9871, 7.0568, 7.4069, 8.5925,10.2970,&
                      9.4269, 9.4681, 9.6689, 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),EA)
!
            D=9.8216d0                                                              !dissociation energy in eV
            ALPHA2=ALPHA1*(1.d0-D/EA)                                               !2nd exponential coefficient
            IF ((EC/EVOLT) <= EA) A=0.d0                                            !to satisfy the fitting constraint
            RXSECTION(J)=A*(((EC/EVOLT)/EA)**ALPHA1)*((1.d0-EA/(EC/EVOLT))**ALPHA2) !SIGMA_EXCHANGE in Angstrons^2
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
        END IF
!
!--possible dissociation reaction with QCT model
        IF (NPM == 3) THEN
          IF (MTYPE==1) THEN
            ! data only for N2+O=>2N+O reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 55) IVIB=55
!
            IMAX=60
            ALLOCATE (VEC_I(0:IMAX))
            VEC_I=(/2.454482d-4,2.427769d-4,2.401402d-4,2.375291d-4,2.348600d-4,2.322024d-4,&
                    2.296469d-4,2.270147d-4,2.243892d-4,2.217688d-4,2.191515d-4,2.166443d-4,&
                    2.140356d-4,2.114257d-4,2.088131d-4,2.061956d-4,2.036905d-4,2.010655d-4,&
                    1.984311d-4,1.957851d-4,1.932510d-4,1.905840d-4,1.878997d-4,1.851937d-4,&
                    1.826079d-4,1.798666d-4,1.772361d-4,1.745758d-4,1.717679d-4,1.689237d-4,&
                    1.661840d-4,1.632698d-4,1.604614d-4,1.576103d-4,1.547141d-4,1.516307d-4,&
                    1.486383d-4,1.455892d-4,1.424788d-4,1.393016d-4,1.360511d-4,1.327185d-4,&
                    1.294563d-4,1.259450d-4,1.224778d-4,1.187449d-4,1.150340d-4,1.111769d-4,&
                    1.071570d-4,1.029527d-4,9.853521d-5,9.402503d-5,8.922437d-5,8.408293d-5,&
                    7.852777d-5,7.259131d-5,6.598628d-5,5.858944d-5,4.994199d-5,3.901624d-5,&
                    1.683039d-5/)
            B=VEC_I(IVIB)  !ro-vibrational energy ladder alpha_1 fitting parameter for N2
            JROT=-0.5d0+DSQRT(0.25d0+(ECR/EVOLT)/B)
            B=2.88d0    !characteristic rot temperature for N2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 0) JROT=0
            IF (JROT > 150) JROT=150
            DEALLOCATE (VEC_I)
!
            IMAX=5
            JMAX=11
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,20,50,100,150/)              !rot levels
            VEC_J=(/0,1,5,7,10,15,20,25,30,50,55/) !vib levels
!
            VALUES=(/ 13.6954,  3.8612,  8.8992, 32.0071, 22.7469,&
                      15.9766, 16.6928, 10.1433, 18.1889, 20.5395,&
                      21.7207, 28.8407, 14.1536, 43.3566, 33.4931,&
                      17.5769, 15.7801, 27.7189, 29.9199, 45.3146,&
                      21.7649, 23.5739, 29.0426, 33.5678, 50.3760,&
                      46.4568, 38.8535, 47.1537, 44.6373, 66.0602,&
                      62.2184, 53.2935, 55.2718, 71.2303, 67.1038,&
                      58.4835, 67.2176, 59.4007, 65.7182, 76.9170,&
                      76.7714, 70.4267, 46.9692, 55.4398, 62.7420,&
                      52.5293, 52.5709, 56.9886,  0.    ,  0.    ,&
                      34.5883, 33.3195, 40.4456,  0.    ,  0.     /)
            VALUES=VALUES/3.d0
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),A)
!
            VALUES=(/-0.8562,  0.1535, -0.4418, -1.1983, -0.7820,&
                     -0.9275, -0.9514, -0.5159, -0.7844, -0.6867,&
                     -0.9779, -1.1816, -0.6282, -1.3260, -0.9121,&
                     -0.7563, -0.6698, -1.0632, -1.0196, -1.0701,&
                     -0.8428, -0.9005, -1.0460, -1.0323, -1.0748,&
                     -1.2716, -1.1245, -1.2475, -1.0968, -1.1749,&
                     -1.3609, -1.2548, -1.2502, -1.3362, -1.1218,&
                     -1.2437, -1.3355, -1.2301, -1.2168, -1.1619,&
                     -1.3797, -1.3113, -1.0134, -1.0469, -0.9061,&
                     -0.9389, -0.9274, -0.9389,  0.    ,  0.    ,&
                     -0.4158, -0.3554, -0.2485,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA1)
!
            VALUES=(/ 2.7421, 1.9663, 2.4130, 3.1189, 2.4395,&
                      2.8119, 2.8250, 2.4530, 2.6669, 2.3048,&
                      2.8314, 3.0095, 2.4938, 3.0070, 2.4202,&
                      2.6106, 2.5168, 2.8200, 2.5957, 2.5188,&
                      2.5600, 2.5917, 2.6326, 2.4992, 2.4386,&
                      2.7926, 2.6603, 2.7038, 2.4336, 2.3362,&
                      2.7358, 2.6093, 2.5622, 2.4924, 2.0266,&
                      2.4224, 2.5002, 2.3485, 2.1483, 1.7466,&
                      2.3265, 2.2610, 1.9071, 1.7467, 1.3238,&
                      0.9529, 0.9338, 0.8428, 0.    , 0.    ,&
                      0.3998, 0.3449, 0.1476, 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(JROT),DFLOAT(IVIB),ALPHA2)
!
            D=9.8216d0                                                            !dissociation energy in eV
            RXSECTION(J)=A*(((EC/EVOLT)/D)**ALPHA1)*((1.d0-D/(EC/EVOLT))**ALPHA2) !SIGMA_DISS in Angstrons^2
!    IF (TEMP < 6.d3)  RXSECTION(J)=RXSECTION(J)*1.d3                      !trick to sample low T rates
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
!
          IF (MTYPE==3) THEN
            ! data only for O2+O=>3O reaction
!
            IF (IVDC(K) == 1) THEN
              IVIB = IPVIB(1,L)
            ELSE
              IVIB = IPVIB(1,M)
            END IF
            IF (IVIB < 0) IVIB=0
            IF (IVIB > 40) IVIB=40
!
            B=2.07d0    !characteristic rot temperature for O2 (from Bird94)
            JROT=-0.5d0+DSQRT(0.25d0+ECR/(BOLTZ*B))
            IF (JROT < 11) JROT=11
            IF (JROT > 231) JROT=231
!
            IMAX=5
            JMAX=6
            ALLOCATE (VEC_I(IMAX),VEC_J(JMAX),VALUES(IMAX*JMAX),ARRAY_TEMP(IMAX,JMAX))
            VEC_I=(/0,10,20,30,40/)         !vib levels
            VEC_J=(/11,51,101,151,201,231/) !rot levels
!
            VALUES=(/-0.2146,-1.4838,-0.9152,-0.6273, 0.3143,&
                     -0.4933,-1.2553,-0.6653,-0.8075, 0.    ,&
                     -0.7825,-1.1472,-0.6812, 0.    , 0.    ,&
                      0.3910, 0.1276, 0.    , 0.    , 0.    ,&
                      0.4847, 0.    , 0.    , 0.    , 0.    ,&
                      0.6097, 0.    , 0.    , 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),ALPHA1)
!
            VALUES=(/ 2.7549, 3.0977, 2.2660, 1.4933, 0.2179,&
                      2.8395, 2.8527, 1.9294, 1.4879, 0.    ,&
                      2.6707, 2.4818, 1.7044, 0.    , 0.    ,&
                      1.9083, 1.8382, 0.    , 0.    , 0.    ,&
                      2.0407, 0.    , 0.    , 0.    , 0.    ,&
                      0.6837, 0.    , 0.    , 0.    , 0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),ALPHA2)
!
            VALUES=(/ 13.3864, 182.0998,115.1628,215.0688, 93.7781,&
                      20.2472, 165.1482, 80.5651,300.0000,  0.    ,&
                      53.4346, 160.0991,252.3550,  0.    ,  0.    ,&
                      28.6315, 198.8150,  0.    ,  0.    ,  0.    ,&
                     106.1007,   0.    ,  0.    ,  0.    ,  0.    ,&
                     128.6496,   0.    ,  0.    ,  0.    ,  0.     /)
            ARRAY_TEMP=RESHAPE(VALUES,(/IMAX,JMAX/))
            CALL INTERPOLATE(VEC_I,VEC_J,IMAX,JMAX,ARRAY_TEMP,DFLOAT(IVIB),DFLOAT(JROT),A)
!
            D=5.21275d0                                                           !dissociation energy in eV
            RXSECTION(J)=A*(((EC/EVOLT)/D)**ALPHA1)*((1.d0-D/(EC/EVOLT))**ALPHA2) !SIGMA_DISS in Angstrons^2
            DEALLOCATE (VEC_I,VEC_J,VALUES,ARRAY_TEMP)
          END IF
        END IF
      END IF
!--------------------------------------------------------------------------------
    END IF
!
  END DO
END IF
!
RETURN
END SUBROUTINE CHECK_RXSECTION
!
!************************************************************************************
!
SUBROUTINE CHECK_REACTION(IKA,N,L,LM,M,LS,LMS,MS,VRR,VR,VRC,VRI,VCM,RML,RMM,IVDC)
!
!
USE MOLECS
USE GAS
USE CALC
USE OUTPUT
USE GEOM
USE OMP_LIB
USE MFDSMC,only:IMF,IMFS,NMFETR,NMFERR,NMFEVR,NMFVTR
!
IMPLICIT NONE
!
!
INTEGER :: J,K,L,LM,M,N,LS,LMS,MS,IKA,JR,KV,ISTE(MNRE),NS,NPM,I,IVDC(MNRE)
REAL(8),EXTERNAL :: GAM
REAL(KIND=8) :: A,ECT,ECR,ECV,EC,VR,VRR,RML,RMM,&
                VDOF1(MMVM),VDOF2(MMVM),EV1(MMVM),EV2(MMVM),ECV1,ECV2,SVDOF1,SVDOF2,ECR1,ECR2,&
                EP,TEMP,ECT2,VRC(3),VCM(3),VRI
INTEGER :: II,IETDX,IERDX(2)
              !
!--A,B,C working variables
!--J,K,KK,MK,KA,KAA,MKK,JR,KR,JS,KV,IA working integers
!--N the collision cell
!--L,S  the molecule numbers
!--LS,MS,KS species
!--IKA reaction indicator
!--ECR rotational energy
!--ECV vibrational energy
!--ECA available energy
!--STER steric factor
!--ISTE third body species
!--NRE number of reactions that may occur for this pair
!--STERT cumulative steric factor
!--THBCELL number of third body molecules in cell
!--WF weighting factor
!--PSTERT probability
!--VRC relative velocity components: can only be changed by recombination
!--VCM center-of-mass velocity components: can only be changed by recombination
!--VRI magnitude of VRC: can only be changed by recombination
!--VRR new value of u^2+v^2+w^2 containg all the energy
!--VR  new value of relative speed
!--RML,RMM molecule mass parameters
!
NS=ICCELL(3,N)     !sampling cell
TEMP=VAR(10,NS)    !sampling cell vibrational temperature
!
ECT=0.5D00*SPM(1,LS,MS)*VRR
ECR=0.D00; ECV=0.D00
EV1=0.D00; EV2=0.D00
ECR1=0.D00; ECR2=0.D00
ECV1=0.D00; ECV2=0.D00
VDOF1=0.D00; VDOF2=0.D00
SVDOF1=0.D00; SVDOF2=0.D00
IF (ISPR(1,LS) > 0) ECR1=PROT(L)
IF (ISPR(1,MS) > 0) ECR2=PROT(M)
ECR=ECR1+ECR2
IF (MMVM > 0) THEN
  IF (ISPV(LS) > 0) THEN
    DO KV=1,ISPV(LS)
      CALL VIB_ENERGY(EV1(KV),IPVIB(KV,L),KV,LS)
      VDOF1(KV)=2.d0*(SPVM(1,KV,LS)/TEMP)/(DEXP(SPVM(1,KV,LS)/TEMP)-1.d0)
      ECV1=ECV1+EV1(KV)
      SVDOF1=SVDOF1+VDOF1(KV)
    END DO
  END IF
  IF (ISPV(MS) > 0) THEN
    DO KV=1,ISPV(MS)
     CALL VIB_ENERGY(EV2(KV),IPVIB(KV,M),KV,MS)
      VDOF2(KV)=2.d0*(SPVM(1,KV,MS)/TEMP)/(DEXP(SPVM(1,KV,MS)/TEMP)-1.d0)
      ECV2=ECV2+EV2(KV)
      SVDOF2=SVDOF2+VDOF2(KV)
    END DO
  END IF
  ECV=ECV1+ECV2
END IF
!
EC=ECT+ECR+ECV
!
!--check if the molecules follow the order given by IREA
! IF (LS /= MS) THEN
  ! IF (LS == IREA(1,IKA)) IVDC=1  !corrected order
  ! IF (MS == IREA(1,IKA)) IVDC=2  !inverted order
! ELSE
  ! IF (ECV1 >= ECV2) IVDC=1  !choose the high energy molecule for dissociation
  ! IF (ECV2 >  ECV1) IVDC=2
! END IF
!-- Han changed this on Oct 15th 2018
! IVDC(IKA) = 1: LS = IREA(1,IKA), MS = IREA(2,IKA)
! IVDC(IKA) = 2: LS = IREA(2,IKA), MS = IREA(1,IKA)
! For dissociation reaction, IVDC(IKA) tell which molecule is the one to be dissociated
!
!
IF (IKA > 0) THEN
  NPM=NREA(1,IKA)+NREA(2,IKA)
!
!--sample pre-reaction vibrational levels
  DO KV=1,MMVM
    I=0; J=0; K=0
    IF (IVDC(IKA) == 1) THEN
      SELECT CASE(NPM)
        CASE(1) !recombination
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=IPVIB(KV,ISTE(JR)) !third body in recombination
        CASE(2) !exchange
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=0
        CASE(3) !dissociation
          I=IPVIB(KV,L)
          J=IPVIB(KV,M)
          K=0
      END SELECT
    ELSE
      SELECT CASE(NPM)
        CASE(1) !recombination
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=IPVIB(KV,ISTE(JR)) !third body in recombination
        CASE(2) !exchange
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=0
        CASE(3) !dissociation
          I=IPVIB(KV,M)
          J=IPVIB(KV,L)
          K=0
      END SELECT
    END IF
    IF (I > 100) I=100
    IF (J > 100) J=100
    IF (K > 100) K=100
    !$omp atomic
    NPVIB(1,IKA,1,KV,I)=NPVIB(1,IKA,1,KV,I)+1  !bin counter
    !$omp atomic
    NPVIB(1,IKA,2,KV,J)=NPVIB(1,IKA,2,KV,J)+1
    !$omp atomic
    NPVIB(1,IKA,3,KV,K)=NPVIB(1,IKA,3,KV,K)+1
    !$omp critical
    IF (IREAC > 0 .and. NPM == 3 .and. IMF .ne. 0 .and. KV == 1 .and. IMFS == 1) THEN
      IETDX = FLOOR(ECT/BOLTZ/FTMP0/0.01D0)+1;
      IETDX=MIN(IETDX,1000)
      NMFETR(IETDX,IKA) = NMFETR(IETDX,IKA) + 1.0d0

      ! IREDX(1) is the one dissociated
      IF (IVDC(IKA) == 1) THEN
        IERDX(1) = FLOOR(ECR1/BOLTZ/FTMP0/0.01D0)+1
        IERDX(2) = FLOOR(ECR2/BOLTZ/FTMP0/0.01D0)+1
      ELSE
        IERDX(2) = FLOOR(ECR1/BOLTZ/FTMP0/0.01D0)+1
        IERDX(1) = FLOOR(ECR2/BOLTZ/FTMP0/0.01D0)+1
      END IF

      DO II = 1,2
        IERDX(II) = MIN(IERDX(II),1000)
        NMFERR(IERDX(II),II,IKA) = NMFERR(IERDX(II),II,IKA)+1.0d0
      END DO

      NMFEVR(I,1,IKA) = NMFEVR(I,1,IKA) + 1.0d0
      NMFEVR(J,2,IKA) = NMFEVR(J,2,IKA) + 1.0d0

      NMFVTR(I,IETDX,1,IKA) = NMFVTR(I,IETDX,1,IKA) + 1.0d0
      NMFVTR(J,IETDX,2,IKA) = NMFVTR(J,IETDX,2,IKA) + 1.0d0
    ENDIF

    IF (NPM == 3 .and. KV == 1 .and. NPM ==3) THEN
      IF (IVDC(K) == 1) THEN
        EVREM(IKA) = EVREM(IKA) + ECV1  ! EVREM in Joule
      ELSE
        EVREM(IKA) = EVREM(IKA) + ECV2
      END IF
    END IF
    !$omp end critical
  END DO
!
!--sample pre-reaction vibrational energies (normalized by 1000*BOLTZ)
  I=ECV/(1.d3*BOLTZ)         !bin index
  IF (I < 100) THEN
    !$omp atomic
    NEVIB(IKA,1,I)=NEVIB(IKA,1,I)+1 !bin accumulation
  END IF
!
!   -----------------------------------------------------------
!--modify molecule states according to post-reaction conditions
  IF (IREAC <= 1) THEN
    EP=EC+REA(5,IKA)         !post-reaction total energy (not accounting for third body yet)
!
    IF (NPM == 1) THEN
!--one post-collision molecule (recombination)
      IPSP(L)=JREA(1,IKA,1)  !molecule L becomes the recombined molecule
      LS=IPSP(L)
      PV(1:3,L)=VCM(1:3)
      IPCELL(M)=-IPCELL(M)   !molecule M is marked for removal
      M=ISTE(JR)             !calculate collision of molecule L with third body
      MS=IPSP(M)
      VRC(1:3)=PV(1:3,L)-PV(1:3,M)
      VRR=VRC(1)**2+VRC(2)**2+VRC(3)**2
      VRI=DSQRT(VRR)
      ECT2=0.5D00*SPM(1,LS,MS)*VRR
      ECR2=0.d0; ECV2=0.d0
      IF (ISPR(1,MS) > 0) ECR2=PROT(M)
      IF (ISPV(MS) > 0) THEN
          DO KV=1,ISPV(MS)
            CALL VIB_ENERGY(A,IPVIB(KV,M),KV,MS)
            ECV2=ECV2+A !DFLOAT(IPVIB(KV,M))*BOLTZ*SPVM(1,KV,MS)
          END DO
      END IF
      EP=EP+ECT2+ECR2+ECV2   !add third body Etra, Erot, and Evib
      IVDC(IKA)=1  !update order of the molecules
    END IF
!
    IF (NPM == 2) THEN
!--two post-collision molecules (exchange)
      IPSP(L)=JREA(1,IKA,1)
      IPSP(M)=JREA(2,IKA,1)
      LS=IPSP(L)
      MS=IPSP(M)
      IVDC(IKA)=1  !update order of the molecules
    END IF
!
    IF (NPM == 3) THEN
!--three post-collision molecules (dissociation)
!$omp critical
      IF (NM >= MNM) CALL EXTEND_MNM(1.1d0)
      NM=NM+1
      LM=NM                    !new born molecule
!$omp end critical
      IF (IVDC(IKA) == 1) THEN
        IPSP(L)=JREA(1,IKA,1)  !update species
        IPSP(LM)=JREA(1,IKA,2)
        LS=IPSP(L)
        LMS=IPSP(LM)
        PX(LM)=PX(L)
        IPCELL(LM)=IPCELL(L)
        PTIM(LM)=PTIM(L)
      ELSE
        IPSP(M)=JREA(1,IKA,1)
        IPSP(LM)=JREA(1,IKA,2)
        MS=IPSP(M)
        LMS=IPSP(LM)
        PX(LM)=PX(M)
        IPCELL(LM)=IPCELL(M)
        PTIM(LM)=PTIM(M)
      END IF
    PROT(LM)=0.d0
    IPVIB(0:MMVM,LM)=0
    END IF
!
!--set all internal energies to zero
    PROT(L)=0.d0
    PROT(M)=0.d0
    IPVIB(0:MMVM,L)=0
    IPVIB(0:MMVM,M)=0
!
!--update collision pair data
    VRR=2.D00*EP/SPM(1,LS,MS)
    VR=DSQRT(VRR)
    IF (NPM < 3) THEN
      RML=SPM(1,LS,MS)/SP(5,MS)
      RMM=SPM(1,LS,MS)/SP(5,LS)
    ELSE
      IF (IVDC(IKA) == 1) THEN
        RML=SPM(1,IREA(1,IKA),MS)/SP(5,MS)
        RMM=SPM(1,IREA(1,IKA),MS)/SP(5,IREA(1,IKA))
      ELSE
        RML=SPM(1,LS,IREA(1,IKA))/SP(5,IREA(1,IKA))
        RMM=SPM(1,LS,IREA(1,IKA))/SP(5,LS)
      END IF
    END IF
    IF (NPM == 1) THEN       !update VCM for recombination reactions
      VCM(1:3)=RML*PV(1:3,L)+RMM*PV(1:3,M)
    END IF
  END IF
!   -----------------------------------------------------------
END IF
!
RETURN
END SUBROUTINE CHECK_REACTION
!
!************************************************************************************
!
SUBROUTINE ADAPT_CELLS
!
!--adapt the sampling cells through the splitting of the divisions into successive levels
!--the collision cells are divisions of the sampling cells
!
USE MOLECS
USE GAS
USE GEOM
USE CALC
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: M,N,L,K,KK,I,J,JJ,MSEG,NSEG,NSEG1,NSEG2,IDT=0 !included IDT
REAL(KIND=8) :: A,B,DDE,DCRIT,DENNR,RANF !--isebasti: included RANF,DENNR
INTEGER, ALLOCATABLE, DIMENSION(:) :: KDIV,NC
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISD
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: XMIN,XMAX,DRAT
!
!--DCRIT  the number density ratio that causes a cell to be subdivided
!--KDIV(N) the number of divisions/subdivisions (cells or further subdivisions) at level N
!--DRAT(N) the contriburion to the density ratio of element N
!--NC(I) the number of sampling cells at level I
!--DDE the width of an element
!--MSEG the maximum number of segments (a segment is the size of the smallest subdivision
!--NSEG1 the (first segment-1) in the subdivision
!--NSEG2 the final segment in the subdivision
!--ISD(N,M) 0,1 for cell,subdivided for level N subdivision
!--DENNR the refence number density
!
DCRIT=4.0d0    !--may be altered !--isebasti: original value is 1.5d0
!
!--determine the level to which the divisions are to be subdivided
!
A=1.D00
IF (IADAPT == 1) DENNR=FND(1)               !--isebasti: included
IF (IADAPT == 2) DENNR=MINVAL(VAR(3,:))     !--isebasti: included
DO N=1,NCELLS
  IF (VAR(3,N)/DENNR > A) A=VAR(3,N)/DENNR  !--isebasti: FND(1) replaced by DENNR
END DO
ILEVEL=0
DO WHILE (A > DCRIT)
  ILEVEL=ILEVEL+1
  A=A/2.D00
END DO
WRITE (*,*) 'ILEVEL =',ILEVEL   !--isebasti: included
WRITE (9,*) 'ILEVEL =',ILEVEL
NSEG=2**ILEVEL
MSEG=NDIV*NSEG
!
ALLOCATE (KDIV(0:ILEVEL),DRAT(MSEG),NC(0:ILEVEL),ISD(0:ILEVEL,MSEG),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR KDIV ARRAY',ERROR
ENDIF
!
DDE=(XB(2)-XB(1))/DFLOAT(MSEG)
DO N=1,MSEG
  A=XB(1)+(DFLOAT(N)-0.5D00)*DDE
  CALL FIND_CELL(A,M,L)
  DRAT(N)=VAR(3,L)/(DENNR*DFLOAT(NSEG)) !--isebasti: FND(1) replaced by DENNR
END DO
!
!--calculate the number of subdivisions at the various levels of subdivision
KDIV=0
!--also the number of sampling cells at each level
NC=0
!
DO N=1,NDIV    !--divisions
  ISD=0
  ISD(0,1)=1
  KDIV(0)=KDIV(0)+1
!  WRITE (9,*) 'DIVISION',N
  DO I=0,ILEVEL  !--level of subdivision
!    WRITE (9,*) 'LEVEL',I
    J=2**I  !--number of possible subdivisions at this level
    JJ=NSEG/J  !--number of segments in a subdivision
    DO M=1,J
!      WRITE (9,*) 'SUBDIVISION',M
      IF (ISD(I,M) == 1) THEN
        NSEG1=(N-1)*NSEG+(M-1)*JJ+1
        NSEG2=NSEG1+JJ-1
        A=0.D00
!        WRITE (9,*) 'NSEG RANGE',NSEG1,NSEG2
        DO L=NSEG1,NSEG2
          A=A+DRAT(L)
        END DO
!        WRITE (9,*) 'DENS CONTRIB',A
        IF (A < DCRIT) THEN
          NC(I)=NC(I)+1
!          WRITE (9,*) 'LEVEL',I,' CELLS TO', NC(I)
        ELSE
          KDIV(I+1)=KDIV(I+1)+2
!          WRITE (9,*) 'LEVEL',I+1,' SUBDIVISIONS TO',KDIV(I+1)
          DO L=NSEG1-(N-1)*NSEG,NSEG2-(N-1)*NSEG
            ISD(I+1,L)=1
          END DO
        END IF
      END IF
    END DO
  END DO
END DO
!
WRITE (9,*) 'KDIV',KDIV
!
WRITE (9,*) 'NC',NC
WRITE (*,*)
WRITE (9,*) 'Number of divisions',NDIV
A=0
NCELLS=0
DO N=0,ILEVEL
  A=A+DFLOAT(NC(N))/(2.D00**N)
  NCELLS=NCELLS+NC(N)
END DO
WRITE (9,*) 'Total divisions from sampling cells',A
WRITE (9,*) 'Adapted sampling cells',NCELLS
NCCELLS=NCELLS*NCIS
WRITE (9,*) 'Adapted collision cells',NCCELLS
!
DEALLOCATE (JDIV,CELL,ICELL,CCELL,ICCELL,COLLS,WCOLLS,CLSEP,VAR,VARSP,CS,CST,STAT=ERROR) !--isebasti: CST included
IF (ERROR /= 0) THEN
  WRITE (*,*)'PROGRAM COULD NOT DEALLOCATE ARRAYS IN ADAPT',ERROR
END IF
!
DO N=0,ILEVEL
  IF (KDIV(N) > MDIV) MDIV=KDIV(N)
END DO
!
ALLOCATE (JDIV(0:ILEVEL,MDIV),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR JDIV ARRAY IN ADAPT',ERROR
ENDIF
!
ALLOCATE (CELL(4,NCELLS),ICELL(NCELLS),CCELL(5,NCCELLS),ICCELL(3,NCCELLS),XMIN(NCCELLS),XMAX(NCCELLS),STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR CELL ARRAYS IN ADAPT',ERROR
ENDIF
!--isebasti: CST included bellow
ALLOCATE (COLLS(NCELLS),WCOLLS(NCELLS),CLSEP(NCELLS),VAR(21,NCELLS),VARSP(0:11,NCELLS,MSP),&
          CS(0:8+MMVM,NCELLS,MSP),CST(0:4,NCELLS),STAT=ERROR)  !--isebasti: correcting CS allocation
IF (ERROR /= 0) THEN
  WRITE (*,*) 'PROGRAM COULD NOT ALLOCATE SPACE FOR SAMPLING ARRAYS IN ADAPT',ERROR
ENDIF
!
NCCELLS=0
NCELLS=0
!
!--set the JDIV arrays and the sampling cells at the various levels of subdivision
KDIV=0
JDIV=0
!
DO N=1,NDIV    !--divisions
  ISD=0
  ISD(0,1)=1
  KDIV(0)=KDIV(0)+1
  DO I=0,ILEVEL  !--level of subdivision
    J=2**I  !--number of possible subdivisions at this level
    JJ=NSEG/J  !--number of segments in a subdivision
    DO M=1,J
      IF (ISD(I,M) == 1) THEN
        NSEG1=(N-1)*NSEG+(M-1)*JJ+1
        NSEG2=NSEG1+JJ-1
        A=0.D00
        DO L=NSEG1,NSEG2
          A=A+DRAT(L)
        END DO
        IF (A < DCRIT) THEN
          NCELLS=NCELLS+1
          XMIN(NCELLS)=XB(1)+DFLOAT(NSEG1-1)*DDE
          XMAX(NCELLS)=XMIN(NCELLS)+DFLOAT(NSEG2-NSEG1+1)*DDE
          WRITE (9,*) NCELLS,I,' XMIN,XMAX',XMIN(NCELLS),XMAX(NCELLS)
          JDIV(I,KDIV(I)-(J-M))=-NCELLS
!          WRITE (9,*) 'JDIV(',I,',',KDIV(I)-(J-M),')=',-NCELLS
        ELSE
          JDIV(I,KDIV(I)-(J-M))=KDIV(I+1)
!          WRITE (9,*) 'JDIV(',I,',',KDIV(I)-(J-M),')=',KDIV(I+1)
          KDIV(I+1)=KDIV(I+1)+2
          DO L=NSEG1-(N-1)*NSEG,NSEG2-(N-1)*NSEG
            ISD(I+1,L)=1
          END DO
        END IF
      END IF
    END DO
  END DO
END DO
!
!--set the other quantities associated with the sampling cells and the collision cells
!
NCCELLS=0
DO N=1,NCELLS
  CELL(1,N)=(XMIN(N)+XMAX(N))/2.D00
  CELL(2,N)=XMIN(N)
  CELL(3,N)=XMAX(N)
  IF (IFX == 0) CELL(4,N)=XMAX(N)-XMIN(N)    !--calculation assumes unit cross-section
  IF (IFX == 1) CELL(4,N)=PI*(XMAX(N)**2-XMIN(N)**2)
  IF (IFX == 2) CELL(4,N)=1.33333333333333333333D00*PI*(XMAX(N)**3-XMIN(N)**3)
  ICELL(N)=NCCELLS
  DO M=1,NCIS
    NCCELLS=NCCELLS+1
    ICCELL(3,NCCELLS)=N
    CCELL(1,NCCELLS)=CELL(4,N)/DFLOAT(NCIS)
    CCELL(3,NCCELLS)=DTM/2.D00
    CCELL(4,NCCELLS)=2.D00*VMPM*SPM(2,1,1)
    CALL ZGF(RANF,IDT)
    CCELL(2,NCCELLS)=RANF
    CCELL(5,NCCELLS)=FTIME
  END DO
END DO

IF (nonVHS .ne. 0) CCELL(4,:) = CCELL(4,:)*1.2D0
!
!--assign the molecules to the cells
!
DO N=1,NM
  CALL FIND_CELL(PX(N),IPCELL(N),JJ)
  M=IPCELL(N)
END DO
!
!--deallocate the local variables
DEALLOCATE (KDIV,NC,ISD,XMIN,XMAX,DRAT,STAT=ERROR)
IF (ERROR /= 0) THEN
  WRITE (*,*)'PROGRAM COULD NOT DEALLOCATE LOCAL ARRAYS IN ADAPT',ERROR
END IF
!
RETURN
!
END SUBROUTINE ADAPT_CELLS
!
!*****************************************************************************
!
SUBROUTINE ENERGY(M,I,TOTEN)
!
!--calculate the total energy (all molecules if I=0, otherwise molecule I)
!--used for diagnostic purposes only
!
USE MOLECS
USE GAS
USE CALC
!
IMPLICIT NONE
!
INTEGER :: K,L,N,I,II,M,IV,KV,J
REAL(KIND=8) :: A,TOTEN,TOTENI,HF(MSP)
!
TOTEN=0.
!
IF ((M == 6).OR.(M == 4)) THEN
  IF (I == 0) THEN
    DO N=1,NM
      IF (IPCELL(N) > 0) THEN
        L=IPSP(N)
        TOTENI=TOTEN
        IF (M == 6) THEN
          IF (L==3) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
          IF (L==4) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,2)
          IF (L==5) TOTEN=TOTEN+1.49D-19
        END IF
        IF (M == 4) THEN
          IF ((L==2).OR.(L == 3)) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
        END IF
        TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
        IF (ISPV(L) > 0) THEN
          DO KV=1,ISPV(L)
          J=IPVIB(KV,N)
          IF (J <0) THEN
            J=-J
            IF (J == 99999) J=0
          END IF
            TOTEN=TOTEN+DFLOAT(J)*BOLTZ*SPVM(1,KV,L)
          END DO
        END IF
      END IF
      IF ((TOTEN-TOTENI) > 1.D-16) WRITE (*,*) 'MOL',N,' ENERGY',TOTEN-TOTENI
    END DO
!
!    WRITE (9,*) 'Total Energy =',TOTEN,NM
!    WRITE (*,*) 'Total Energy =',TOTEN,NM
  ELSE
    N=I
    IF (IPCELL(N) > 0) THEN
      L=IPSP(N)
      IF (M == 6) THEN
        IF (L==3) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
        IF (L==4) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,2)
        IF (L==5) TOTEN=TOTEN+1.49D-19
      END IF
      IF ((M == 4).OR.(M > 8)) THEN
!        IF ((L==2).OR.(L == 3)) TOTEN=TOTEN+0.5D00*BOLTZ*SPVM(4,1,1)
      END IF
      TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
      IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
      IF (ISPV(L) > 0) THEN
        DO KV=1,ISPV(L)
          J=IPVIB(KV,N)
          IF (J <0) THEN
            J=-J
            IF (J == 99999) J=0
          END IF
          TOTEN=TOTEN+DFLOAT(J)*BOLTZ*SPVM(1,KV,L)
        END DO
      END IF
    END IF
  END IF
END IF
IF (M == 8) THEN
  !HF(1)=0.d0         !H2 entalpy of formation [kcal/mol]; from O Conaire et al (2004)
  !HF(2)=52.098d0     !H
  HF(1)=0.d0         !N2 entalpy of formation [kcal/mol]; from O Conaire et al (2004)
  HF(2)=112.954d0    !N
  HF(3)=0.d0         !O2
  HF(4)=59.56d0      !O
  HF(5)=8.91d0       !OH
  HF(6)=-57.77d0     !H2O
  HF(7)=3.d0         !HO2
  HF(8)=-94.05d0     !CO2   !from wikipedia = 393509 J/mol
!  HF(8)=0.d0         !Ar
  HF(:)=HF(:)*4184.d3/AVOG  !converting to [J/molec]
  IF (I == 0) THEN
    DO N=1,NM
      IF (IPCELL(N) > 0) THEN  !energies are normalized by BOLTZ
        TOTENI=TOTEN
        L=IPSP(N)
        TOTEN=TOTEN+HF(L)/BOLTZ
        TOTEN=TOTEN+0.5D00*(SP(5,L)/BOLTZ)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)/BOLTZ
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (IPVIB(K,N) < 0) THEN
              WRITE (*,*) 'Dissociation marked molecule still in flow',N,IPVIB(K,N)
!              STOP
            END IF
            CALL VIB_ENERGY(A,IPVIB(K,N),K,L)
            TOTEN=TOTEN+A/BOLTZ !IPVIB(K,N)*SPVM(1,K,L)
         END DO
        END IF
        IF ((TOTEN-TOTENI) > 1.d-16/BOLTZ) WRITE (*,*) 'MOL',N,' SPECIES',L,' ENERGY/BOLTZ',TOTEN-TOTENI
      END IF
    END DO
!
!    WRITE (9,*) 'Total Energy =',TOTEN,NM
!    WRITE (*,*) 'Total Energy =',TOTEN,NM
  ELSE
    N=I
    IF (IPCELL(N) > 0) THEN
      L=IPSP(N)
      IF (IPCELL(N) > 0) THEN
        L=IPSP(N)
        IF (L==2) TOTEN=TOTEN+3.62D-19
        IF (L==4) TOTEN=TOTEN+4.14D-19
        IF (L==5) TOTEN=TOTEN+0.65D-19
        IF (L==6) TOTEN=TOTEN-4.02D-19
        IF (L==7) TOTEN=TOTEN+0.17D-19
        TOTEN=TOTEN+0.5D00*SP(5,L)*(PV(1,N)**2+PV(2,N)**2+PV(3,N)**2)
        IF (ISPR(1,L) > 0) TOTEN=TOTEN+PROT(N)
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IV=IPVIB(K,N)
            IF (IV < 0) THEN
              IF (IV == -99999) THEN
                IV=0
              ELSE
                IV=-IV
              END IF
              IF (L == 1) TOTEN=TOTEN+7.24E-19
              IF (L == 3) TOTEN=TOTEN+8.28E-19
              IF (L == 5) TOTEN=TOTEN+7.76E-19
              IF (L == 6) TOTEN=TOTEN+8.29E-19
              IF (L == 7) TOTEN=TOTEN+3.45E-19
            END IF
            TOTEN=TOTEN+IV*BOLTZ*SPVM(1,K,L)
          END DO
        END IF
        IF (ABS(TOTEN) > 1.D-16) THEN
          CONTINUE
        END IF
      END IF
    END IF
  END IF
END IF
!
RETURN
!
END SUBROUTINE ENERGY
!
!************************************************************************************
!
SUBROUTINE COLLISIONS
!
!--calculate the collisions
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
USE MFDSMC,only: NMFET0,NMFER0,NMFEV0,NMFVT0,NMFEV,NMFET,NMFER,NMFVT, &
  IMF, IMFS, IMFpair
!
IMPLICIT NONE
!
INTEGER :: N,NS,NN,M,MM,L,LL,K,KK,KT,J,I,II,III,NSP,MAXLEV,IV,IVP,NSEL,KV,LS,MS,KS,JS,IKA,LZ,KL,IS,IREC,&
           NLOOP,IA,IDISS,IE,IEX,JJ,NPRI,LIMLEV,KVV,KW,ISP1,ISP2,ISP3,INIL,INIM,JI,LV,IVM,NMC,NVM,LSI,&
           JX,MOLA,KR,JKV,NSC,KKV,NUMEXR,IAX,NSTEP,ISWITCH,NPM,LM,LMS,IVDC(MNRE),ISHUF(3),IT,IDT=0,IVDC0
REAL(KIND=8) :: A,AA,AAA,AB,B,BB,BBB,ASEL,DTC,SEP,VRI,VR,VRR,ECT,ET0,EVIB,ECC,ZV,COLT,ERM,C,OC,SD,D,CVR,PROB,&
                RML,RMM,ECTOT,ETI,EREC,ET2,XMIN,XMAX,WFC,CENI,CENF,VRRT,TCOLT,EA,DEN,SVDOF(3,3),RANF,DT(ITMAX),& !--isebasti: included RANF,DT
                RXSECTION(MNRE),SXSECTION,TXSECTION,&
                QNU,A1,A2,B1,B2,C1,C2,EROT,EV,SIGMA_REF,EV_POST,SUMF,E,F,EL,ED,EF,S !ME-QCT variables
REAL(KIND=8),DIMENSION(0:100) :: VTXSECTION
REAL(KIND=8),DIMENSION(3) :: VRC,VCM,VRCP,VRCT
REAL(KIND=8) :: ECR(2),EVIBEV,BMAX,REST_DOF
INTEGER :: IETDX,IERDX(2),IEVDX(2),IVPS(2)
logical :: IREACSP
!
!--N,M,K working integer
!--LS,MS,KS,JS molecular species
!--VRC components of the relative velocity
!--RML,RMM molecule mass parameters
!--VCM components of the center of mass velocity
!--VRCP post-collision components of the relative velocity
!--SEP the collision partner separation
!--VRR the square of the relative speed
!--VR the relative speed
!--ECT relative translational energy (update after different process)
!--ET0 relative intitial translational energy in eV
!--EVIB vibrational energy
!--ECC collision energy (rel trans +vib)
!--MAXLEV maximum vibrational level
!--ZV vibration collision number
!--COLT collision temperature (rel trans + vib)
!--TCOLT collision temperature based rel trans
!--SDF the number of degrees of freedom associated with the collision
!--ERM rotational energy
!--NSEL integer number of selections
!--CVR product of collision cross-section and relative speed
!--PROB a probability
!--IKA reaction indicator  ( 0 for no reaction, >0 for reaction )
!--KT third body molecule code
!--ECTOT energy added at recmbination
!--IREC initially 0, becomes 1 if a recombination occurs
!--NPRI the priority collision molecule (the last one rejected as a previous collision molecule)
!--WFC weighting factor in the cell
!--IEX is the reaction that occurs (1 if only one is possible)
!--EA activation energy
!--IVDC order of the two colliding particles
!
! REAL(8) :: MVELC(3), MSVELC(3), VARVELC(3)
!--MVELC mean velocity of the cell, mass averaged
!--MSVELC mean square of the velocity
!--VARIANCE

NUMEXR=0

!
!$omp parallel &
!$omp private(idt) & !always =0 outside parallel regions
!$omp private(n,ns,dtc,nn,wfc,aaa,c,xmin,xmax,asel,nsel,i,ie,kl,npri,iii,k,l,lm,m,ranf,ll) &
!$omp private(sep,j,mm,a,kk,vrc,vrr,vri,vr,cvr,ect,nsp,kv,evib,ecc,maxlev,colt,b,zv,ii,iv,ivp,prob,svdof) &
!$omp private(erm,vcm,vrcp,oc,sd,d,ls,lms,ms,ivdc,ivdc0,ishuf,rml,rmm,iswitch,idiss,ks,js,limlev,lz,iex,irec) &
!$omp private(tcolt,kt,aa,bb,kvv,ji,lsi,ivm,nmc,nvm,ectot,psf,jj,ea,den,iax,jx,ika,npm,nstep,it,dt) &
!$omp private(SXSECTION,RXSECTION,VTXSECTION,TXSECTION) &
!$omp private(QNU,A1,A2,B1,B2,C1,C2,E,F,EL,ED,ET0,EROT,EV,SIGMA_REF,EV_POST,SUMF,EF,S) &
!$omp private(ECR,EVIBEV,IVPS,BMAX,REST_DOF)&
!$omp private(IETDX,IEVDX,IERDX,IREACSP) &
!$omp reduction(+:ndissoc,ndissl,trecomb,nrecomb,treacl,treacg,tnex,tforex,trevex) & !Q-K
!$omp reduction(+:totdup,totcol,pcolls,tcol,cscr,colls,wcolls,clsep,reac,npvib) &
!$omp reduction(+:NMFET0,NMFER0,NMFEV0,NMFVT0,NMFEV,NMFET,NMFER,NMFVT)
!$    idt=omp_get_thread_num()  !thread id
!$omp do !schedule (dynamic,NCCELLS/32)

! OPEN(19, FILE="VelocityVAR.txt", ACCESS="APPEND")
DO N=1,NCCELLS
!
  IF (FTIME-CCELL(5,N) > CCELL(3,N)) THEN
!
!--calculate collisions appropriate to  time DTC
    DTC=2.D00*CCELL(3,N)
    IF (ICCELL(2,N) > 1) THEN  !--no collisions calculated if there are less than two molecules in collision cell
      NN=ICCELL(3,N)
      WFC=1.D00
      IF ((IWF == 1).AND.(IVB == 0)) WFC=1.D00+WFM*CELL(1,NN)**IFX
      CCELL(5,N)=CCELL(5,N)+DTC
      IF (IVB == 0) AAA=CCELL(1,N)
      IF (IVB == 1) THEN
        C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV*NCIS)
        XMIN=XB(1)+DFLOAT(N-1)*C
        XMAX=XMIN+C
        WFC=1.D00+WFM*(0.5D00*(XMIN+XMAX))**IFX
        IF (IFX == 0) AAA=XMAX-XMIN
        IF (IFX == 1) AAA=PI*(XMAX**2-XMIN**2)  !--assumes unit length of full cylinder
        IF (IFX == 2) AAA=1.33333333333333333333D00*PI*(XMAX**3-XMIN**3)    !--flow is in the full sphere
      END IF
!
!--these statements implement the N(N-1) scheme
      ASEL=0.5D00*ICCELL(2,N)*(ICCELL(2,N)-1)*WFC*FNUM*CCELL(4,N)*DTC/AAA+CCELL(2,N)
      NSEL=ASEL
      CCELL(2,N)=ASEL-DFLOAT(NSEL)
!
      IF (NSEL > 0) THEN
!
        ! MVELC = 0.0d0
        ! MSVELC = 0.0d0
        ! DO K = 1, ICCELL(2,N)
          ! M = ICREF(ICCELL(1,N) + K)
          ! LS = ABS(IPSP(M))
          ! DO J = 1,3
            ! MVELC(J) = MVELC(J) + PV(J,M)
            ! MSVELC(J) = MSVELC(J) + PV(J,M)**2
          ! END DO
        ! END DO
        ! DO J=1,3
          ! MVELC(J) = MVELC(J) / ICCELL(2,N)
          ! MSVELC(J) = MSVELC(J) / ICCELL(2,N)
          ! VARVELC(J) = MSVELC(J) - MVELC(J)**2
        ! END DO
        ! WRITE(19,'(I3,1X,20G13.5)') N, (MVELC(J),J=1,3), (MSVELC(J),J=1,3), (VARVELC(J),J=1,3)
!
        NS=0   !--counts the number of selections
        KL=0   !--becomes 1 if it is the last selection
        NPRI=0
!
        DO KL=1,NSEL
          NS=NS+1
          IE=0  !--becomes >0 if an inelastic collision occurs
!
          IF (ICCELL(2,N) == 2) THEN
            K=1+ICCELL(1,N)
            L=ICREF(K)
            K=2+ICCELL(1,N)
            M=ICREF(K)
            IF (M == IPCP(L)) THEN
              III = 1
              CCELL(5,N)=CCELL(5,N)-DTC
            END IF
          ELSE
            K=NPRI
            IF (K > 0) THEN
              L=K
              NPRI=0
            ELSE
              CALL ZGF(RANF,IDT)
              K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
              L=ICREF(K)
            END IF
!--one molecule has been selected at random
            IF (NNC == 0) THEN
!--select the collision partner at random
              M=L
              DO WHILE (M == L)
                CALL ZGF(RANF,IDT)
                K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                M=ICREF(K)
              END DO
            ELSE
!--select the nearest from the total number (< 30) or a random 30
              IF (ICCELL(2,N) < 30) THEN
                LL=ICCELL(2,N)
              ELSE
                LL=30
              END IF
              SEP=1.D10
              M=0
              DO J=1,LL
                IF (LL < 30) THEN
                  K=J+ICCELL(1,N)
                ELSE
                  CALL ZGF(RANF,IDT)
                  K=INT(RANF*DFLOAT(ICCELL(2,N)))+ICCELL(1,N)+1
                END IF
                MM=ICREF(K)
                IF (MM /= L) THEN   !exclude the already selected molecule
                  IF (MM /= IPCP(L)) THEN   !exclude the previous collision partner
                    A=DABS(PX(L)-PX(MM))
                    IF ((A < SEP).AND.(A >1.D-8*DDIV)) THEN
                      M=MM
                      SEP=A
                    END IF
                  ELSE
                    NPRI=MM    !MM is made the priority molecule
                  END IF
                END IF
              END DO
            END IF
          END IF
!
!--Calculate the relative velocity
          DO KK=1,3
            VRC(KK)=PV(KK,L)-PV(KK,M)
          END DO
          VRR=VRC(1)*VRC(1)+VRC(2)*VRC(2)+VRC(3)*VRC(3)
          VR=DSQRT(VRR)
          VRI=VR
!
!--Prepare variables
          LS=ABS(IPSP(L))
          MS=ABS(IPSP(M))
          RML=SPM(1,LS,MS)/SP(5,MS)
          RMM=SPM(1,LS,MS)/SP(5,LS)
          DO KK=1,3
            VCM(KK)=RML*PV(KK,L)+RMM*PV(KK,M)
          END DO
!
          IF (IREAC .EQ. 2) THEN
            DO I = 1,MNRE
              IF ( (LS == IREA(1,I)) .and. (MS == IREA(2,I)) ) THEN
                IREACSP = .true.
                EXIT
              ELSE IF ( (MS == IREA(1,I)) .and. (LS == IREA(2,I))) THEN
                IREACSP = .true.
                EXIT
              ELSE
                IREACSP = .false.
              END IF
            END DO
          ELSE
            IREACSP = .true.
          END IF

!--Simple gas
          IF (MSP == 1) THEN
            CVR=VR*CXSS*((2.D00*BOLTZ*SP(2,1)/(RMAS*VRR))**(SP(3,1)-0.5D00))*RGFS
            IF (CVR > CCELL(4,N)) CCELL(4,N)=CVR
            CALL ZGF(RANF,IDT)
            IF (RANF < CVR/CCELL(4,N)) THEN
!--the collision occurs
              IF ((M == IPCP(L)).AND.(L == IPCP(M))) TOTDUP=TOTDUP+1.D00
              TOTCOL=TOTCOL+1.D00
              TCOL(1,1)=TCOL(1,1)+2.D00
              CSCR(1,1)=CSCR(1,1)+2.d0*CVR
              COLLS(NN)=COLLS(NN)+1.D000
              WCOLLS(NN)=WCOLLS(NN)+WFC
              IF(NN == NSPDF) PCOLLS=PCOLLS+WFC
              SEP=DABS(PX(L)-PX(M))
              CLSEP(NN)=CLSEP(NN)+SEP
              IF ((ISPR(1,1) > 0)) THEN
!--Larsen-Borgnakke serial redistribution
                ECT=0.5D00*RMAS*VRR
                DO NSP=1,2
!--consider the molecules in turn
                  IF (NSP == 1) THEN
                    K=L
                  ELSE
                    K=M
                  END IF
                  IF (MMVM > 0) THEN
                    IF (ISPV(1) > 0) THEN
                      DO KV=1,ISPV(1)
                        CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,1)
                        ECC=ECT+EVIB
                        CALL VIB_LEVEL(ECC,MAXLEV,KV,1)
                        IF (SPVM(3,KV,1) > 0.) THEN
!--note quantizing of COLT in the following statements
                          IF (IZV == 1) THEN
                            COLT=(DFLOAT(MAXLEV)*SPVM(1,KV,1))/(3.5D00-SPM(3,1,1))        !--quantized collision temperature
                          ELSE
                            COLT=VAR(8,NN)        !--isebasti: Ttrans, according to DSMC.f90 (2013)
                          END IF
                          B=SPVM(4,KV,1)/SPVM(3,KV,1)    !Tdiss/Tref
                          A=SPVM(4,KV,1)/COLT            !Tdiss/T
       ZV=(A**SPM(3,1,1))*(SPVM(2,KV,1)*(B**(-SPM(3,1,1))))**(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                        ELSE
                          ZV=SPVM(2,KV,1)
                        END IF
                        CALL ZGF(RANF,IDT)
                        IF (1.D00/ZV > RANF) THEN
                          II=0
                          DO WHILE (II == 0)
                            CALL ZGF(RANF,IDT)
                            IV=RANF*(MAXLEV+0.99999999D00)
                            IPVIB(KV,K)=IV
                            CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,1)
                            IF (EVIB < ECC) THEN
                              CALL INELASTIC_VT(1,1,EVIB,ECC,5.0d0-2.0d0*SPM(3,1,1),PROB)
!--PROB is the probability ratio of eqn (5.61)
                              CALL ZGF(RANF,IDT)
                              IF (PROB > RANF) II=1
                            END IF
                          END DO
                          ECT=ECC-EVIB
                        END IF
                      END DO
                    END IF
                  END IF
!--now rotation of this molecule
                  IF (ISPR(1,1) > 0) THEN
                    IF (ISPR(2,1) == 0) THEN
                      B=1.D00/SPR(1,1)
                    ELSE    !--use molecule rather than mean value
                      COLT=ECC/((2.5D00-SP(3,1))*BOLTZ)
                      B=1.D00/(SPR(1,1)+SPR(2,1)*COLT+SPR(3,1)*COLT*COLT)
                    END IF
                    CALL ZGF(RANF,IDT)
                    IF (B > RANF) THEN
                      ECC=ECT+PROT(K)
                      CALL INELASTIC_RT(1,1,ECC,ERM,IDT)
                      PROT(K)=ERM*ECC
                      ECT=ECC-PROT(K)
                    END IF
                  END IF
                END DO
!--adjust VR for the change in energy
                VR=SQRT(2.D00*ECT/SPM(1,1,1))
              END IF
!--end of L-B redistribution
              DO KK=1,3
                VCM(KK)=0.5d0*(PV(KK,L)+PV(KK,M))
              END DO
!--calculate new velocities
              CALL VHSS(1,1,VR,VRC,VRCP,IDT)
              PV(1:3,L)=VCM(1:3)+0.5d0*VRCP(1:3)
              PV(1:3,M)=VCM(1:3)-0.5d0*VRCP(1:3)
              IPCP(L)=M
              IPCP(M)=L
            END IF   !--collision occurrence
!
!--Gas Mixture
          ELSE IF (IPCELL(L)>0.AND.IPCELL(M)>0.AND.IPVIB(0,L)>=0.AND.IPVIB(0,M)>=0.AND. IREACSP) THEN !trick
            !IPCELL<0 for recombined molecule marked for removal
            !IPVIB(0,L)<0 for molecule marked for dissociation
!
!--Calculate the total scattering (VHS) cross-section
            ECT=0.5D00*SPM(1,LS,MS)*VRR         !collision translational energy in J
            ET0=ECT/EVOLT                        !convert to eV
            CALL CALC_TOTXSEC(LS, MS, VR, VRR, ET0, -1.0d0, SXSECTION, BMAX, CVR)


! get rotational energy and vibrational energy, ECR1 is always the one with lower sp number
! 1 is always the one with lower sp number  or the one with higher Ev
            ECR = 0.0d0; IVPS = 0
            IF (LS < MS ) THEN
              IF (ISPV(LS) > 0) THEN
                ECR(1) = PROT(L); IVPS(1) = IPVIB(1,L)
              END IF
              IF (ISPV(MS) > 0) THEN
                ECR(2) = PROT(M); IVPS(2) = IPVIB(1,M)
              END IF
            ELSE IF (LS > MS) THEN
              IF (ISPV(LS) > 0) THEN
                ECR(2) = PROT(L); IVPS(2) = IPVIB(1,L)
              END IF
              IF (ISPV(MS) > 0) THEN
                ECR(1) = PROT(M); IVPS(1) = IPVIB(1,M)
              END IF
            ELSE ! same specie
              IF (ISPV(LS) > 0) THEN
                IF (IPVIB(1,L) >= IPVIB(1,M)) THEN ! select one with higher Ev
                  ECR(1) = PROT(L); IVPS(1) = IPVIB(1,L)
                  ECR(2) = PROT(M); IVPS(2) = IPVIB(1,M)
                ELSE
                  ECR(1) = PROT(M); IVPS(1) = IPVIB(1,M)
                  ECR(2) = PROT(L); IVPS(2) = IPVIB(1,L)
                END IF
              END IF
            END IF
            ! Sample rotational, translational and vibrational energy
            IF ((IREAC .eq. 2) .and. (IMF .ne. 0) .and. (MNRE <= 2) .and. (IMFS == 1)) THEN
              IETDX = FLOOR(ECT/BOLTZ/FTMP0/0.010D0)+1
              IETDX = MIN(IETDX,1000)
              NMFET0(IETDX,IMFpair(LS,MS)) = NMFET0(IETDX,IMFpair(LS,MS)) + 1.0D0
              ! KK=1 the one with lower smaller sp number
              DO KK = 1,2
                K = KK
                ! for same specie, X(:,2,:) is always zero
                ! rotational energy
                IERDX(KK) = FLOOR(ECR(KK)/BOLTZ/FTMP0/0.01D00)+1
                IERDX(KK) = MIN(IERDX(KK),1000)
                NMFER0(IERDX(KK),K,IMFpair(LS,MS)) = NMFER0(IERDX(KK),K,IMFpair(LS,MS)) + 1.0D0
                ! vibrational energy
                IEVDX(KK) = IVPS(KK)
                IF (IEVDX(KK) > 100)    IEVDX(KK) = 100
                NMFEV0(IEVDX(KK),K,IMFpair(LS,MS)) = NMFEV0(IEVDX(KK),K,IMFpair(LS,MS)) + 1.0d0
                NMFVT0(IEVDX(KK),IETDX,K,IMFpair(LS,MS)) = NMFVT0(IEVDX(KK),IETDX,K,IMFpair(LS,MS))+1.0D0
              END DO
            END IF

!
!--Calculate the VT (ME-QCT) cross-sections (only for O2+O and N2+O collisions)
            J=0
            IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) J=1 !N2-O collision
            IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) J=3 !O2-O collision
            ! J is only larger than 0 for the above two pairs

            ! Calcuate pre-collision vibrational energy
            ! IVP is the vibrational level of the molecule which is more
            ! likely to be dissociated
            ! For J>0, predefined QCT model *MIGHT* be used later
            ! For J<0, QCT model will never be used
            IVP = -1
            EVIB = 0.0d0
            IF (LS == J) THEN
              IVP=IPVIB(1,L);
            ELSE IF  (MS == J) THEN
              IVP=IPVIB(1,M)
            ELSE IF ((ISPV(LS) == 1) .and. (ISPV(MS) == 0))THEN
              IVP = IPVIB(1,L)
              J = -LS
            ELSE IF ((ISPV(MS) == 1) .and. (ISPV(LS) == 0))THEN
              IVP = IPVIB(1,M)
              J = -MS
            ELSE ! both are diatom, select the one with higher vibrational energy to be dissociated
              IF (IPVIB(1,L) >= IPVIB(1,M)) THEN
                J = -LS
                IVP = IPVIB(1,L)
              ELSE
                J = -MS
                IVP = IPVIB(1,M)
              ENDIF
            END IF
            IF (IVP > 0)THEN
              CALL VIB_ENERGY(EVIB,IVP,1,ABS(J))
            END IF
            ! the following values are used in ME-QCT-VT model
            ECC=(ECT+EVIB)/EVOLT                !available collision energy in eV

            ! Han: add nonVHS cross section model for N2+O
            ! This model has dependency on vibrational level
            IF (nonVHS == 1) THEN
              IF ((LS == 1 .and. MS == 4) .or. (LS ==4 .and. MS == 1)) THEN
                CALL CALC_TOTXSEC(LS, MS, VR, VRR, ET0, EVIB, SXSECTION,BMAX,CVR)
              END IF
            END IF

!--Calculate the reaction cross-sections
            IVDC=0                 !to track the order of reacting species
            RXSECTION=0.d0
            IF (MNRE>0 .and. IMF .ne.  0 .and. QCTMODEL == 2) THEN
              ! precalculate Reaction cross section before collision occur
              ! This is for the case when QCTMODEL is used
              ! Reaction cross sections might become larger than others
              CALL CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,VRR,ECT,SXSECTION,IVDC,IDT)
            END IF
!


            ! QCT VT cross sections
            VTXSECTION=0.d0
            IF (QCTMODEL==2.AND.GASCODE==8.AND.J>0) THEN
!
              IF(J==1) THEN  !N2+O collisions
                SIGMA_REF=14.5351d0/3.d0 !24.2788d0/3.d0
                QNU=0.7074d0 !0.3994d0
                A1=1.0570d0  !0.1253d0
                A2=-0.0503d0 !0.0221d0
                B1=8.2246d0  !5.8506d0
                B2=0.1639d0  !0.1398d0
                C1=0.4170d0  !0.2493d0
                C2=-0.3942d0 !-0.2075d0
                D=0.6974d0   !0.6988d0
                E=-7.3212d0  !-5.0131d0
                F=-0.0273d0  !-0.0868d0
                ED=9.8216d0  !E_dissociation in eV
              END IF
!
              IF(J==3) THEN  !O2+O collisions
                SIGMA_REF=6.0108d0/27.d0
                QNU=0.6390d0
                A1=0.0331d0
                A2=0.0016d0
                B1=2.0921d0
                B2=0.2459d0
                C1=-0.1357d0
                C2=0.4272d0
                D=0.5859d0
                E=-3.1009d0
                F=23.3814d0
                ED=5.21275d0  !E_dissociation in eV
              END IF
!
!
              A=A1+A2*ECC
              B=B1+B2*ECC
              C=C1+C2*ECC
!
              SUMF=0.d0
              EL=DMIN1(ECC,ED)                    !in eV
!
              IF(J==1) THEN !N2+O collisions
                AB=ECT+EVIB
                CALL VIB_LEVEL(AB,MAXLEV,1,J)              !max level based on the collision energy
                DO IV=0,MAXLEV
                  CALL VIB_ENERGY(EV_POST,IV,1,J)
                  AB=DABS(EVIB-EV_POST)/EVOLT              !energy difference in eV
                  BB=DABS(AB/C)
                  EF=(EV_POST/EVOLT)/ECC                   !energy ratio
                  IF(EF > 1.d0) EF=1.d0                    !to avoid round off errors when EV_POST=ECC
                  AAA=EL-(EVIB/EVOLT)
                  BBB=EL-(EV_POST/EVOLT)
                  S=A*AB+B*DEXP(-BB**D)+E*(AAA*BBB)**F     !surprisal function
                  VTXSECTION(IV)=((1.d0-EF)**QNU)*DEXP(S)
                  SUMF=SUMF+(1.d0-EF)**QNU                 !normalization constant
                END DO
              END IF
!
              IF(J==3) THEN !O2+O collisions
                MAXLEV=IVMODEL(J,2)
                DO IV=0,MAXLEV
                  CALL VIB_ENERGY(EV_POST,IV,1,J)
                  AB=DABS(EVIB-EV_POST)/EVOLT             !energy difference in eV
                  BB=DABS(AB/C)
                  EF=(EV_POST/EVOLT)/ECC                  !energy ratio
                  IF(EF > 1.d0) EF=1.d0                   !to avoid round off errors when EV_POST=ECC
                  AAA=(EV_POST/EVOLT)/ED-1.d0
                  BBB=(EVIB/EVOLT)/ED-1.d0
                  S=A*AB+B*DEXP(-BB**D)+E*DEXP(-F*AAA*BBB) !surprisal function
                  VTXSECTION(IV)=((1.d0-EF)**QNU)*DEXP(S)
                  SUMF=SUMF+(1.d0-EF)**QNU                 !normalization constant
                END DO
              END IF
!
              VTXSECTION(:)=VTXSECTION(:)/SUMF                     !normalized
              VTXSECTION(:)=SIGMA_REF*ET0**(QNU-1.d0)*VTXSECTION(:) !SIGMA_ME-QCT-VT in Angstrons^2
            END IF
!
!--Define the total cross-section
            TXSECTION=DMAX1(SXSECTION,SUM(VTXSECTION(:)))      !trick
            IF (MNRE>0 .and. IMF .ne.  0 .and. QCTMODEL == 2) THEN
              TXSECTION=DMAX1(TXSECTION,SUM(RXSECTION(:)))
            END IF
            !TXSECTION=DMAX1(TXSECTION,SUM(VTXSECTION(:)))     !total cross-section in Angstrons^2

            IF (TXSECTION*VRI/1.d20 > CCELL(4,N)) THEN
              !IF (((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) .and. nonVHS ==1) THEN
                !WRITE(*,*) TXSECTION*VRI/1.d20/CCELL(4,N)
              !ENDIF
              ! WRITE(*,*) TXSECTION*VRI/1.d20/CCELL(4,N)
              CCELL(4,N)=TXSECTION*VRI/1.d20 !update maximum value in (m2)*(m/s)
            END IF
!
!--NTC approach based on total cross-section
            PROB=(TXSECTION*VRI/1.d20)/CCELL(4,N)

            CALL ZGF(RANF,IDT)
            IF (PROB > RANF) THEN !collision occurs
!
!--Sample and initialize collision data
!-------------------------------------------------
              IF ((M==IPCP(L)).AND.(L==IPCP(M))) TOTDUP=TOTDUP+1.D00
              TOTCOL=TOTCOL+1.D00
              TCOL(LS,MS)=TCOL(LS,MS)+1.D00
              TCOL(MS,LS)=TCOL(MS,LS)+1.D00
              CSCR(LS,MS)=CSCR(LS,MS)+TXSECTION*VRI/1.d20 !accumulate species-dependent SigmaTotal*Cr (in m3/s)
              CSCR(MS,LS)=CSCR(MS,LS)+TXSECTION*VRI/1.d20
              COLLS(NN)=COLLS(NN)+1.D00
              WCOLLS(NN)=WCOLLS(NN)+WFC
              IF(NN == NSPDF) PCOLLS=PCOLLS+WFC
              SEP=DABS(PX(L)-PX(M))
              CLSEP(NN)=CLSEP(NN)+SEP

              !
              !--Sample internal states and energy
              !-------------------------------------------------
              ! sample t,r,v energy for particles that collide
              IF ((IREAC .eq. 2) .and. (IMF .ne. 0) .and. (MNRE <= 2) .and. (IMFS == 1)) THEN
                NMFET(IETDX,IMFpair(LS,MS)) = NMFET(IETDX,IMFpair(LS,MS)) + 1.0D0
                ! KK=1 the one with lower smaller sp number
                DO KK = 1,2
                  K = KK
                  ! for same specie, X(:,2,:) is always zero
                  ! rotational energy
                  NMFER(IERDX(KK),K,IMFpair(LS,MS)) = NMFER(IERDX(KK),K,IMFpair(LS,MS)) + 1.0D0
                  ! vibrational energy
                  NMFEV(IEVDX(KK),K,IMFpair(LS,MS)) = NMFEV(IEVDX(KK),K,IMFpair(LS,MS)) + 1.0d0
                  NMFVT(IEVDX(KK),IETDX,K,IMFpair(LS,MS)) = NMFVT(IEVDX(KK),IETDX,K,IMFpair(LS,MS))+1.0D0
                END DO
              END IF
!
!--Check chemical reaction
!-------------------------------------------------
              IKA=0                  !becomes >0 if reaction occurs
              JX=2                   !standard number of post reaction species
              LMS=1                  !initialize
              IF (MNRE>0) THEN
                IF (QCTMODEL .ne. 2) THEN
                  RXSECTION = 0.0d0
                  CALL CHECK_RXSECTION(RXSECTION,N,L,M,LS,MS,VRR,ECT,SXSECTION,IVDC,IDT)
                  ! for a chemical reaction K
                  ! IVDC(K) = 1:  LS = IREA(1,K), MS = IREA(2,K)
                  ! IVDC(K) = 2:  LS = IREA(2,K), MS = IREA(1,K)
                  ! For the case IREA(1,K) = IREA(2,K)
                  ! IVDC(K) = 1:  Molecule L is the one to be dissociated
                  ! IVDC(K) = 2:  Molecule M is the one to be dissociated
                END IF
                SUMF=SUM(RXSECTION(:)) !total reaction cross-section in Angstron^2

                PROB=SUMF/TXSECTION
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) THEN  !reaction occurs
                  IF (PROB > 1.d0) write(*,*) 'Prob_RXN above unity:',PROB,SUMF,TXSECTION
                  J=0
                  A=0.d0
                  CALL ZGF(RANF,IDT)
                  DO WHILE (A < RANF)
                    J=J+1
                    A=A+RXSECTION(J)/SUMF
                  END DO
                  IKA=IRCD(J,LS,MS) !reaction IKA occurs
                  CALL CHECK_REACTION(IKA,N,L,LM,M,LS,LMS,MS,VRR,VR,VRC,VRI,VCM,RML,RMM,IVDC)
                  ! VRR, VR become the velocites make 0.5*mr*VRR = Etot
                  ! For recombination, VRC is changed from A-B relative velocity to AB-C, the same applies to VCM
                END IF
              END IF
!
!--Update variables for VRT energy redistribution
!-------------------------------------------------
              IF (IKA /= 0) THEN
                NPM=NREA(1,IKA)+NREA(2,IKA)
                IF (NPM == 3) JX=3          !dissociation
                IF (IREAC == 0) THEN
                  REAC(IKA)=REAC(IKA)+1.D00 !count reactions in entire domain
                ELSE
                  IF (NN == NSPDF) REAC(IKA)=REAC(IKA)+1.D00 !count reactions only in NSPDF cell
                END IF
                IF (IPRS == 1) THEN
                  DO IT=1,ITMAX
                    DT(IT)=DABS(VAR(8,NN)-FPTEMP(IT))
                  END DO
                  IT=MINLOC(DT,1)
                END IF
                IF (IREAC == 2) THEN
                  IKA=0
                  JX=2
                  IVDC = 0
                  IVDC0 = 0
                ELSE
                  IVDC0 = IVDC(IKA)
                END IF
              ELSE
                IVDC0 = 0
              END IF

!             !ECT is updated here
              ! This is the very tricky part, for the case without chemical reaction, VRR is the original one
              ! otherwise ECT is the total energy containg all modes
              ECT=(0.5D00*SPM(1,LS,MS)*VRR)-ELACK
              IF (ECT > 0.d0) THEN
                ELACK=0.d0
              ELSE
                ELACK=-ECT+1.d-26
                ECT=1.d-26
              END IF
!
              CALL ZGF(RANF,IDT)
              IF (RANF > 0.5d0) THEN  !shuffle the sequence of molecules in energy redistribution
                IAX=1;  KK=JX; IVM=1
              ELSE
                IAX=JX; KK=1;  IVM=-1
              END IF
!
!--Check VT energy redistribution
!-------------------------------------------------
              J=0
              IF ((LS==1.AND.MS==4).OR.(LS==4.AND.MS==1)) J=1 !N2-O collision
              IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) J=3 !O2-O collision
!
              IF (IREAC<=1.AND.QCTMODEL==2.AND.IKA==0.AND.GASCODE==8.AND.J>0) THEN
                ! The following code is run under the condition that
                !  1. Reaction can happend (IREA <=1)
                !  2. QCTMODEL is on
                !  3. Reaction doesn't occur (IKA = 0)
                !  4. GASCODE == 8
                !  5. the species is either N2-O or O2-O
              !--Vibrational relaxation based on ME-QCT model
!
                SUMF=SUM(VTXSECTION(:))      !total VT cross-section (SIGMA_ME-QCT-VT,T)
                PROB=SUMF/TXSECTION          !SIGMA_ME-QCT-VT,T to SIGMA_TOTAL ratio
                IF (PROB > 1.d0) write(*,*) 'Prob_VT above unity:',PROB
!
                CALL ZGF(RANF,IDT)
                IF (PROB > RANF) THEN        !vibrational relaxation occurs
                  II=0
                  MAXLEV=IVMODEL(J,2)
                  DO WHILE (II == 0)
                    CALL ZGF(RANF,IDT)
                    IV=RANF*(MAXLEV+0.99999999D00)
                    PROB=VTXSECTION(IV)/SUMF !SIGMA_ME-QCT-VT to SIGMA_ME-QCT-VT,T ratio
                    CALL ZGF(RANF,IDT)
                    IF (PROB > RANF) THEN    !a vibrational level is accepted
                      CALL VIB_ENERGY(EV_POST,IV,1,J)
                      II=1
                    END IF
                  END DO
!
                  IF (LS == J) IPVIB(1,L)=IV
                  IF (MS == J) IPVIB(1,M)=IV
!
                  IF (IRELAX /= 0) ECT=(ECC*EVOLT)-EV_POST  !remaining energy available for redistribution; comment this line for isothermal relaxations
                END IF
!
              ELSE
!-------------------------------------------------
                !--Larsen-Borgnakke serial vibrational energy redistribution
                REST_DOF = 0.0d0
                IF ((IKA /= 0).AND.(IPRS == 2)) THEN
                      A=VAR(8,NN)  !sampling cell translational temperature
                      SVDOF=0.d0
                      IF (ISPV(LS) > 0) THEN
                        DO KV=1,ISPV(LS)
                          SVDOF(1,KV)=2.d0*(SPVM(1,KV,LS)/A)/(DEXP(SPVM(1,KV,LS)/A)-1.d0)
                        END DO
                      END IF
                      IF (ISPV(MS) > 0) THEN
                        DO KV=1,ISPV(MS)
                          SVDOF(2,KV)=2.d0*(SPVM(1,KV,MS)/A)/(DEXP(SPVM(1,KV,MS)/A)-1.d0)
                        END DO
                      END IF
                      IF (JX == 3) THEN
                        IF (ISPV(LMS) > 0) THEN
                          DO KV=1,ISPV(LMS)
                            SVDOF(3,KV)=2.d0*(SPVM(1,KV,LMS)/A)/(DEXP(SPVM(1,KV,LMS)/A)-1.d0)
                          END DO
                        END IF
                      END IF
                      REST_DOF=5.d0-2.d0*SPM(3,LS,MS)   !translational dofs
                      REST_DOF=REST_DOF+ISPR(1,LS)+ISPR(1,MS) !rotational dofs
                      REST_DOF=REST_DOF+SUM(SVDOF)            !vibrational dofs
                END IF
!
                DO NSP=IAX,KK,IVM
                      IF (IVDC0 < 2) THEN
                        IF (NSP == 1) THEN
                          K=L ; KS=LS ; JS=MS
                        ELSE IF (NSP == 2) THEN
                          K=M ; KS=MS ; JS=LS
                        ELSE
                          K=LM ; KS=LMS ; JS=MS  !new born molecule from dissociation reaction
                        END IF
                      ELSE
                        IF (NSP == 1) THEN
                          K=M ; KS=MS ; JS=LS
                        ELSE IF (NSP == 2) THEN
                          K=L ; KS=LS ; JS=MS
                        ELSE
                          K=LM ; KS=LMS ; JS=LS  !new born molecule from dissociation reaction
                        END IF
                      END IF
!
                      IF (MMVM > 0) THEN
                        IF (ISPV(KS) > 0) THEN
                          IF (ISPV(KS) == 3) THEN  !shuffle the sequence of vibrational modes over serial LB
                            CALL ZGF(RANF,IDT)
                            JI=RANF*(6.d0+0.9999999d0)
                            SELECT CASE (JI)
                              CASE DEFAULT; ISHUF(1)=1; ISHUF(2)=2; ISHUF(3)=3;
                              CASE (1); ISHUF(1)=1; ISHUF(2)=2; ISHUF(3)=3;
                              CASE (2); ISHUF(1)=1; ISHUF(2)=3; ISHUF(3)=2;
                              CASE (3); ISHUF(1)=2; ISHUF(2)=3; ISHUF(3)=1;
                              CASE (4); ISHUF(1)=2; ISHUF(2)=1; ISHUF(3)=3;
                              CASE (5); ISHUF(1)=3; ISHUF(2)=1; ISHUF(3)=2;
                              CASE (6); ISHUF(1)=3; ISHUF(2)=2; ISHUF(3)=1;
                            END SELECT
                          END IF
                          DO I=1,ISPV(KS)
                            KV=I
                            IF (ISPV(KS) == 3) KV=ISHUF(I)
                            IF (IPVIB(KV,K) >= 0) THEN
!--do not redistribute to a dissociating molecule marked for removal or for a recombined molecule
                              IF (IKA .eq. 0) THEN
                                ! IF IKA == 0, ECC is the total energy before chemical reaction + reaction heat
                                CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
                                ECC=ECT+EVIB
                              END IF
                              CALL VIB_LEVEL(ECC,MAXLEV,KV,KS)
                              IF (SPVM(3,KV,KS) > 0.) THEN
                                IF (IZV == 1) THEN
                                  COLT=(DFLOAT(MAXLEV)*SPVM(1,KV,KS))/(3.5D00-SPM(3,KS,JS))  !--quantized collision temperature
                                ELSE
                                  COLT=VAR(8,NN)   !--isebasti: Ttrans, according to DSMC.f90 (2013)
                                END IF
                                B=SPVM(4,KV,KS)/SPVM(3,KV,KS)    !Tdiss/Tref
                                A=SPVM(4,KV,KS)/COLT             !Tdiss/T
       ZV=(A**SPM(3,KS,JS))*(SPVM(2,KV,KS)*(B**(-SPM(3,KS,JS))))**(((A**0.3333333D00)-1.D00)/((B**0.33333D00)-1.D00))
                              ELSE
                                IF (SPVM(3,1,KS) == -1.) THEN
                                  ZV=SPVM(2,KV,KS)
                                ELSE IF (SPVM(3,1,KS) == -2.) THEN
                                  IF (KS == JS) THEN
                                    ZV=SPVM(2,1,KS)
                                  ELSE
                                    ZV=SPVM(2,2,KS)
                                  END IF
                                END IF
                              END IF
                              IF (IREAC<=1) THEN
!VAR(8,NN)=0.154*SPVM(1,KV,KS)
!do II=0,20
!VAR(10,NN)=3.5d0*SPVM(1,KV,KS)*DFLOAT(ii)/20.d0
                                IF (VAR(8,NN) >= VAR(11,NN)) THEN
                                  A=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  A=0.5d0*A*A*DEXP(SPVM(1,KV,KS)/VAR(8,NN))
                                  B=5.d0-2.d0*SPM(3,LS,MS)
                                  C=B/(B+A) !C is the Zcorrection
                                ELSE
                                  A=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  B=2.d0*(SPVM(1,KV,KS)/VAR(10,NN))/(DEXP(SPVM(1,KV,KS)/VAR(10,NN))-1.d0)
                                  IF(VAR(10,NN) < 1.d-3) B=0.d0
                                  C=5.d0-2.d0*SPM(3,LS,MS)
                                  D=B*VAR(10,NN)+C*VAR(8,NN)
                                  CALL ROOTF_SECANT(2,KS,LS,MS,VAR(11,NN),D,E)
                                  C=(C*(VAR(8,NN)-E))/(A*VAR(8,NN)-B*VAR(10,NN)) !C is the Zcorrection
                                END IF
!write(*,*) VAR(10,NN)/SPVM(1,KV,KS),VAR(8,NN)/SPVM(1,KV,KS),E/SPVM(1,KV,KS),C
!end do
!stop
                                IF ((LS==3.AND.MS==4).OR.(LS==4.AND.MS==3)) THEN
                                  !calibrated with O2-O MeQct model
                                  IF (VAR(8,NN)<= 1.d3) ZV=C*(60.01d0+0.04268d0*VAR(8,NN)-2.29d-5*VAR(8,NN)**2.d0)
                                  IF (VAR(8,NN) > 1.d3) ZV=C*DEXP(-0.1372d0*DLOG(VAR(8,NN))+5.328d0)
                                END IF
                                !IF (LS==3.AND.MS==3) THEN
                                !  !calibrated with O2-O2 ibraguimova2013 data
                                !  A=VAR(8,NN)**(-1.d0/3.d0)                          !T^(-1/3)
                                !  B=7.206d0-289.4d0*A+4919.d0*A**2.d0-2.23d4*A**3.d0 !log10(Zv)
                                !  ZV=C*10.d0**B
                                !END IF
                              END IF
                              CALL ZGF(RANF,IDT)
                              IF ((1.D00/ZV > RANF).OR.(IKA /= 0)) THEN
!-- a vibrational redistribution is always made after a chemical reaction
                                II=0
                                nstep=0
                                IF ((IKA /= 0).AND.(IPRS == 2)) THEN  !modified LB for reacting collisions
                                  B=2.d0*(SPVM(1,KV,KS)/VAR(8,NN))/(DEXP(SPVM(1,KV,KS)/VAR(8,NN))-1.d0)
                                  REST_DOF=REST_DOF-B
                                  A=REST_DOF    !remaining degrees of freedom
                                ELSE
                                  A=5.d0-2.d0*SPM(3,KS,JS)              !nonreacting collisions
                                END IF

                                DO WHILE ((II == 0).and.(nstep < 100000))
                                  nstep=nstep+1
                                  if (nstep > 99000) then
                                    write (*,*) nstep,ecc,maxlev
                                    stop
                                  end if
                                  CALL ZGF(RANF,IDT)
                                  IV=RANF*(MAXLEV+0.99999999D00)
                                  IPVIB(KV,K)=IV
                                  CALL VIB_ENERGY(EVIB,IPVIB(KV,K),KV,KS)
                                  !EVIB=DFLOAT(IV)*BOLTZ*SPVM(1,KV,KS)
                                  IF (EVIB < ECC) THEN
                                    CALL INELASTIC_VT(KS, JS, EVIB, ECC, A, PROB)
                                    ! PROB=(1.D00-EVIB/ECC)**(A*0.5d0-1.d0)
 !--PROB is the probability ratio of eqn (5.61) and (5.45)
                                    CALL ZGF(RANF,IDT)
                                    IF (PROB > RANF) II=1
                                  END IF
                                END DO
                                IF (IRELAX /= 0) ECT=ECC-EVIB !remaining energy available for redistribution; comment this line for isothermal relaxations
                              END IF
                            END IF
                          END DO
                        END IF
                      END IF
                END DO
              END IF  !--end of VT exchange
!
!--sample post-collisional vibrational levels
!-------------------------------------------------
              IF (IKA /= 0) THEN
                    DO KV=1,MMVM
                      IF (IVDC0 == 1) THEN
                        SELECT CASE(JX)
                        CASE(2) !recombination and exchange
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,M)
                          K=0
                        CASE(3) !dissociation
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,LM) !new born
                          K=IPVIB(KV,M)  !third body
                        END SELECT
                      ELSE
                        SELECT CASE(JX)
                        CASE(2) !recombination and exchange
                          I=IPVIB(KV,L)
                          J=IPVIB(KV,M)
                          K=0
                        CASE(3) !dissociation
                          I=IPVIB(KV,M)
                          J=IPVIB(KV,LM) !new born
                          K=IPVIB(KV,L)  !third body
                        END SELECT
                      END IF
                      IF (I > 100) I=100
                      IF (J > 100) J=100
                      IF (K > 100) K=100
                      NPVIB(2,IKA,1,KV,I)=NPVIB(2,IKA,1,KV,I)+1  !bin counter
                      NPVIB(2,IKA,2,KV,J)=NPVIB(2,IKA,2,KV,J)+1
                      NPVIB(2,IKA,3,KV,K)=NPVIB(2,IKA,3,KV,K)+1
                    END DO
              END IF
!
!--Check RT energy redistribution (LB model)
!-------------------------------------------------
              DO NSP=IAX,KK,IVM
                    IF  (IVDC0 < 2) THEN
                      IF (NSP == 1) THEN
                        K=L ; KS=LS ; JS=MS
                      ELSE IF (NSP == 2) THEN
                        K=M ; KS=MS ; JS=LS
                      ELSE
                        K=LM ; KS=LMS ; JS=MS  !new born molecule from dissociation reaction
                      END IF
                    ELSE
                      IF (NSP == 1) THEN
                        K=M ; KS=MS ; JS=LS
                      ELSE IF (NSP == 2) THEN
                        K=L ; KS=LS ; JS=MS
                      ELSE
                        K=LM ; KS=LMS ; JS=LS  !new born molecule from dissociation reaction
                      END IF
                    END IF
                    IF (ISPR(1,KS) > 0) THEN
                      IF (ISPR(2,KS) == 0) THEN
                        B=SPM(7,KS,JS) !Zr
                      ELSE
                        COLT=100.d0/VAR(8,NN) !T*/Ttrans; using SMILE data
                        B=20.d0/(1.d0+(0.5d0*PI**(3.d0/2.d0))*DSQRT(COLT)+PI*(1.d0+PI/4.d0)*COLT) !Zr
                      END IF
                      A=5.d0-2.d0*SPM(3,KS,LS)
                      C=A/(A+ISPR(1,KS)+ISPR(1,JS)) !Zr correction
                      A=B*C                         !Corrected Zr
                      CALL ZGF(RANF,IDT)
                      IF ((1.d0/A > RANF).OR.(IKA /= 0)) THEN
                        ECC=ECT+PROT(K)
                        CALL INELASTIC_RT(KS, JS, ECC, ERM, IDT)
                        PROT(K)=ERM*ECC
                        ECT=ECC-PROT(K)
                      END IF
                    END IF
              END DO
!
!--adjust VR after energy redistribution
!-------------------------------------------------
              IF (JX == 2) THEN
                VR=DSQRT(2.D00*ECT/SPM(1,LS,MS))
              ELSE
                ! Reaction off
                VR=DSQRT(2.D00*ECT/SPM(1,IREA(1,IKA),IREA(2,IKA))) !same as in SPARTA
              END IF
!
!--calculate new velocities
!-------------------------------------------------
              ! At this point for exchange or dissocation reaction
              ! VRC and VRI are still the relative velocity of the two molecule
              ! For recombination, it's the relative velocity of the
              ! new molecule and the third body
              ! The values are changed in CHECK_REACTION
              VRC = VRC/VRI*VR             ! rescale to match current magnitude

              ! The current implementation may not handle recombination well for EXP scattering
              IF (JX == 2) THEN
                CALL SCATTER_MOL(LS, MS, BMAX, ET0, VRI, VR, VRC, VRCP, IDT)
              ELSE
                CALL SCATTER_MOL(IREA(1,IKA), IREA(2,IKA), BMAX,ET0, VRI, VR, VRC, VRCP, IDT)
              END IF

              PV(1:3,L)=VCM(1:3)+RMM*VRCP(1:3)
              PV(1:3,M)=VCM(1:3)-RML*VRCP(1:3)
              IF (JX == 3) THEN
                IF (IVDC0 == 1) THEN
                  PV(1:3,LM)=PV(1:3,L)
                  IPCP(LM)=L
                ELSE
                  PV(1:3,LM)=PV(1:3,M)
                  IPCP(LM)=M
                END IF
              END IF
              IPCP(L)=M
              IPCP(M)=L
!
              !--to sustain a bimodal translational-rotatioanal distribution
              !IF (IGS==4) THEN
              !  DO NSP=1,2
              !    IF (NSP == 1) THEN
              !      K=L ; KS=LS
              !    ELSE
              !      K=M ; KS=MS
              !    END IF
              !    CALL ZGF(RANF,IDT)
              !    IF (RANF > 0.5d0) THEN
              !      A=FTMP(1)   !based on translational reference temperature
              !    ELSE
              !      A=FRTMP(1)  !based on rotational reference temperature
              !    END IF
              !    B=SQRT(2.D00*BOLTZ*A/SP(5,KS))
              !    DO J=1,3
              !      CALL RVELC(PV(J,K),C,B,IDT) !for bimodal distribution
              !    END DO
              !    PV(1,K)=PV(1,K)+VFX(1)
              !    PV(2,K)=PV(2,K)+VFY(1)
              !    IF (ISPR(1,KS) > 0) CALL SROT(KS,A,PROT(K),IDT) !for bimodal distribution
              !  END DO
              !END IF
!
            END IF   !--end of collision
          END IF   !--separate simple gas / mixture cases
        END DO   !--do for all candidate collisions
        ! MVELC = 0.0d0
        ! MSVELC = 0.0d0
        ! DO K = 1, ICCELL(2,N)
          ! M = ICREF(ICCELL(1,N) + K)
          ! LS = ABS(IPSP(M))
          ! DO J = 1,3
            ! MVELC(J) = MVELC(J) + PV(J,M)
            ! MSVELC(J) = MSVELC(J) + PV(J,M)**2
          ! END DO
        ! END DO
        ! DO J=1,3
          ! MVELC(J) = MVELC(J) / ICCELL(2,N)
          ! MSVELC(J) = MSVELC(J) / ICCELL(2,N)
          ! VARVELC(J) = MSVELC(J) - MVELC(J)**2
        ! END DO
        ! WRITE(19,'(I3,1X,20G13.5)') N, (MVELC(J),J=1,3), (MSVELC(J),J=1,3), (VARVELC(J),J=1,3)
      END IF
    END IF
  END IF
END DO
! WRITE(9,*)
! WRITE(9,*)
! CLOSE(9)
!$omp end do
!$omp end parallel
!
!--remove any recombined atoms
N=0
DO WHILE (N < NM)
  N=N+1
  IF (IPCELL(N) < 0) THEN
    CALL REMOVE_MOL(N)
    N=N-1
  END IF
END DO
!
!--dissociate marked molecules
!IF (MNRE > 0) CALL DISSOCIATION
!
RETURN
END SUBROUTINE COLLISIONS
!
!*****************************************************************************
!
SUBROUTINE CALC_TOTXSEC(LS,MS,VR,VRR,ET0,EVIB,TOTXSEC,BMAX,CVR)
!
!--calculate the total cross sections, i.e. collision cross sections
! EVIB is the vibrational energy of N2 molecule for N2-O nonVHS model
! TOTXSEC in angstrom^2
! CVR,VR,VRR, EVIB in si
! ET0 in eV
USE GAS, only : INONVHS, SPM
USE CALC, only: BOLTZ, EVOLT, PI
USE EXPCOL, only: EXPCOL_TOTXSEC
IMPLICIT NONE
INTEGER,intent(in) :: LS,MS
REAL(8),intent(in) :: VRR, EVIB, VR, ET0
REAL(8) :: EVIBEV
REAL(8),intent(out) ::TOTXSEC, CVR, BMAX
REAL(8),PARAMETER :: CTOT(8) = (/46.039504155886434, 0.051004420857885, -0.584551057667247, &
                              &  -0.002969283806997, 0.281618756125794, 0.030181202283512, &
                              &  -0.436592532266083, 0.152224780739684/)
INTEGER :: IERROR

! ET = 0.5d0*SPM(1,LS,MS)*VRR/EVOLT  ! collisional energy in eV

IF (INONVHS(LS,MS) == 0 .or. (INONVHS(LS,MS) == 1 .and. EVIB < 0.0D0)) THEN
  TOTXSEC = SPM(2,LS,MS)*((2.D00*BOLTZ*SPM(5,LS,MS)/(SPM(1,LS,MS)*VRR))**(SPM(3,LS,MS)-0.5D00))*SPM(6,LS,MS)*1.0d20
ELSE IF (INONVHS(LS,MS) == 1) THEN
  EVIBEV = EVIB/EVOLT  ! convert to eV
  IF (EVIBEV >= 6.9d0  ) EVIBEV = 6.9d0  ! extrapolate
  TOTXSEC = DLOG(ET0) - (CTOT(8) + EVIBEV*(CTOT(7) + EVIBEV*CTOT(6)))
  TOTXSEC = (CTOT(3)+CTOT(4)*EVIBEV*EVIBEV)*DTANH(CTOT(5)*TOTXSEC) + CTOT(2)*EVIBEV
  TOTXSEC = CTOT(1)*DEXP(TOTXSEC)
ELSE IF (INONVHS(LS,MS) == 2) THEN
  CALL EXPCOL_TOTXSEC(LS, MS, ET0, TOTXSEC, IERROR)
END IF
BMAX = DSQRT(TOTXSEC/PI)  !angstrom
CVR = TOTXSEC/1.0d20*VR

END SUBROUTINE CALC_TOTXSEC
!
!*****************************************************************************
!
SUBROUTINE SCATTER_MOL(LS,MS,BMAX,ET0,VRI,VR,VRC,VRCP,IDT)
!-- a wrapper of different scattering model
USE GAS, only: INONVHS
USE EXPCOL, only: EXPCOL_Scatter
IMPLICIT NONE
INTEGER,INTENT(IN) :: LS, MS
REAL(8),INTENT(IN) :: BMAX,ET0,VRI, VR, VRC(3)
REAL(8) :: VRCP(3)
INTEGER :: IDT
! VRI original relative speed
! VR, VRC new speed without scattering
! all the following functions should only change VRCP
IF (INONVHS(LS,MS) == 0 ) THEN
  CALL VHSS(LS, MS, VR, VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 1) THEN
  CALL N2OScatter(LS, MS, BMAX,VRI, VR,VRC, VRCP, IDT)
ELSE IF (INONVHS(LS,MS) == 2) THEN
  CALL EXPCOL_Scatter(LS,MS,BMAX,ET0,VR,VRC,VRCP,IDT)
END IF

END SUBROUTINE SCATTER_MOL
!
!*****************************************************************************
!
SUBROUTINE VHSS(LS,MS,VR,VRC,VRCP,IDT)
!
!--calculate new scattering angles based on VHS/VSS model
!
USE GAS
USE CALC
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: LS,MS,IDT
REAL(KIND=8) :: A,B,C,D,OC,SD,VR,VRC(3),VRCP(3),RANF
!
IF (ABS(SPM(8,LS,MS)-1.) < 1.d-3) THEN
!--use the VHS logic
  CALL ZGF(RANF,IDT)
  B=2.D00*RANF-1.D00
!--B is the cosine of a random elevation angle
  A=DSQRT(1.D00-B*B)
  VRCP(1)=B*VR
  CALL ZGF(RANF,IDT)
  C=2.D00*PI*RANF
!--C is a random azimuth angle
  VRCP(2)=A*DCOS(C)*VR
  VRCP(3)=A*DSIN(C)*VR
ELSE
!--use the VSS logic
  CALL ZGF(RANF,IDT)
  B=2.D00*(RANF**SPM(8,LS,MS))-1.D00  !isebasti: SP(4,1) was used instead of SPM(8,LS,MS)
!--B is the cosine of the deflection angle for the VSS model (eqn (11.8)
  A=DSQRT(1.D00-B*B)
  CALL ZGF(RANF,IDT)
  C=2.D00*PI*RANF
  OC=DCOS(C)
  SD=DSIN(C)
  D=SQRT(VRC(2)**2+VRC(3)**2)
  VRCP(1)=B*VRC(1)+A*SD*D
  VRCP(2)=B*VRC(2)+A*(VR*VRC(3)*OC-VRC(1)*VRC(2)*SD)/D
  VRCP(3)=B*VRC(3)-A*(VR*VRC(2)*OC+VRC(1)*VRC(3)*SD)/D
!--the post-collision rel. velocity components are based on eqn (2.22)
END IF
RETURN
END SUBROUTINE VHSS
!
!*****************************************************************************
SUBROUTINE N2OScatter(LS,MS,BMAX, VRI,VR,VRC,VRCP,IDT)
!
! -- calculate new scattering angles based on non VHS model
! -- only for N2+O collision
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
INTEGER :: LS,MS,IDT
REAL(KIND=8) :: B,VR,VRC(3),VRCP(3),RANF,C,D,VRI
REAL(8) :: COF(6), BMAX, ET, RMASS, EVIBEV,C1,C2,CHI,CCHI,SCHI
REAL(8) :: EPSI,CEPSI,SEPSI

COF(1) = 1.783071942310433
COF(2) = 0.016980219000487/1000.0D0
COF(3) = 2.846002511624816
COF(4) = 0.294962630947434/1000.0D0
COF(5) = 2.189954664950628
COF(6) = 0.078701684994745


EVIBEV =  0.063577602025999  !Erv(v=0,J=0) for N2

CALL ZGF(RANF,IDT)
B = DSQRT(RANF)*BMAX
C1 = COF(3)*DEXP(-COF(4)*VRI) +COF(5)*DEXP(-COF(6)*EVIBEV)
C2 = COF(1) + COF(2)*VRI
CHI = C2*(1.0d0 - 1.0d0/(C1+DLOG(2.0D0))*(B-DLOG(DCOSH(B-C1))))
! Scattering angle should only be a function of precollision conditions
! this two condition shouldn't occur, for safety we set it here
IF (CHI > PI) CHI = PI
IF (CHI < 0) CHI = 0.0D0
CCHI = DCOS(CHI); SCHI = DSIN(CHI)

CALL ZGF(RANF,IDT)
EPSI = RANF*2.0D0*PI
CEPSI = DCOS(EPSI)
SEPSI = DSIN(EPSI)

D=DSQRT(VRC(2)**2+VRC(3)**2)
VRCP(1) = CCHI*VRC(1) + SCHI*SEPSI*D
VRCP(2) = CCHI*VRC(2) + SCHI*(VR*VRC(3)*CEPSI-VRC(1)*VRC(2)*SEPSI)/D
VRCP(3) = CCHI*VRC(3) - SCHI*(VR*VRC(2)*CEPSI+VRC(1)*VRC(3)*SEPSI)/D

RETURN
END SUBROUTINE N2OScatter
!
! -----------------------------------------------------------------
SUBROUTINE N2OScatter_old(LS,MS,BMAX,VR,VRC,VRCP,IDT)
!
! -- calculate new scattering angles based on non VHS model
! -- only for N2+O collision
! -- old version: which uses the new VR to calculate scattering angle
! -- this doesn't make sense, but it is the one Han used in PRF paper
!
USE GAS
USE CALC
!
IMPLICIT NONE
!
INTEGER :: LS,MS,IDT
REAL(KIND=8) :: B,VR,VRC(3),VRCP(3),RANF,C,D
REAL(8) :: COF(6), BMAX, ET, RMASS, EVIBEV,C1,C2,CHI,CCHI,SCHI
REAL(8) :: EPSI,CEPSI,SEPSI

COF(1) = 1.783071942310433
COF(2) = 0.016980219000487/1000.0D0
COF(3) = 2.846002511624816
COF(4) = 0.294962630947434/1000.0D0
COF(5) = 2.189954664950628
COF(6) = 0.078701684994745


EVIBEV =  0.063577602025999  !Erv(v=0,J=0) for N2

CALL ZGF(RANF,IDT)
B = DSQRT(RANF)*BMAX  !sampled impact parameter
C1 = COF(3)*DEXP(-COF(4)*VR) +COF(5)*DEXP(-COF(6)*EVIBEV)
C2 = COF(1) + COF(2)*VR
CHI = C2*(1.0d0 - 1.0d0/(C1+DLOG(2.0D0))*(B-DLOG(DCOSH(B-C1))))
! this two condition shouldn't occur, for safety we set it here
IF (CHI > PI) CHI = PI
IF (CHI < 0) CHI = 0.0D0
CCHI = DCOS(CHI); SCHI = DSIN(CHI)

CALL ZGF(RANF,IDT)
EPSI = RANF*2.0D0*PI
CEPSI = DCOS(EPSI)
SEPSI = DSIN(EPSI)

D=DSQRT(VRC(2)**2+VRC(3)**2)
VRCP(1) = CCHI*VRC(1) + SCHI*SEPSI*D
VRCP(2) = CCHI*VRC(2) + SCHI*(VR*VRC(3)*CEPSI-VRC(1)*VRC(2)*SEPSI)/D
VRCP(3) = CCHI*VRC(3) - SCHI*(VR*VRC(2)*CEPSI+VRC(1)*VRC(3)*SEPSI)/D

RETURN
END SUBROUTINE N2OScatter_old
!
!*****************************************************************************
SUBROUTINE THERMOPHORETIC
!
!--author: isebasti
!--calculate thermophoretic force based on Gimelshein et al, AIAA 2005-766;
!--assume one montionless particle of constant temperature per sampling cell;
!--weighting factors are not considered in the current version;
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: L,K
INTEGER(KIND=8) :: N,NC,NS
REAL(KIND=8) :: ALPHAM,ALPHAT,ALPHAR,ALPHAV,A,B,D,E,F,G,EVIB,FVIB,DEGENV,TH
REAL(KIND=8),DIMENSION(3) :: VR,TF
REAL(KIND=8),DIMENSION(NCELLS) :: PTEMP
REAL(KIND=8),DIMENSION(3,NCELLS) :: PPV
!
!--N molecule number
!--NC,NS,L,K working integer
!--A,B,E,F,FVIB auxiliar variable
!--ALPHAM accommodation coefficient for momentum
!--ALPHAT accommodation coefficient for translational energy
!--ALPHAR accommodation coefficient for rotational energy
!--ALPHAV accommodation coefficient for vibrational energy
!--PPV(1:3,N) velocity components of the macroscopic particle (up)
!--PTEMP(:) temperature of the macroscopic particle (Tp)
!--VR(1:3) components of relative velocity btw particle and molecule
!--G magnitude of the relative velocity
!--EVIB molecule vibrational energy
!--DEGENV degeneracy of a vibrational mode
!--TF(1:3,N) components of thermophoretic force per unit cross section area (pi*Rp^2)
!--TH(N) thermophoretic heat transfer per unit cross section area (pi*Rp^2)
!--CST(0:4,:) sum of thermoforetic properties in a sampling cell
!             0 number sum, 1:3 thermophoretic force compoments, 4 heat transfer
!
ALPHAM=1.d0
ALPHAT=1.d0
ALPHAR=1.d0
ALPHAV=1.d0
DEGENV=1.d0 !assuming all vibrational modes are non-degenerated
!
PPV(1:3,:)=0.d0
PTEMP(:)=273.d0
!
!------------------------------------------------------------------------------
!--loop over molecules (maybe it would be better a loop only over the macroscopic particles!)
!$omp parallel &
!$omp private(n,nc,ns,l,vr,g,a,b,d,e,tf,evib,fvib,th) &
!$omp reduction(+:cst)
!$omp do
DO N=1,NM
  NC=IPCELL(N)        !collision cell
  NS=ICCELL(3,NC)     !sampling cell
  L=IPSP(N)           !species of molecule n
!
  VR(1:3)=PV(1:3,N)-PPV(1:3,NS)
  G=SQRT(VR(1)*VR(1)+VR(2)*VR(2)+VR(3)*VR(3))
!
  A=SP(5,L)*FNUM/CELL(4,NS)
  B=1.d0+4.d0*ALPHAM*(1.d0-ALPHAT)/9.d0
  D=ALPHAM*ALPHAT*SPI*SQRT(2.d0*BOLTZ*PTEMP(NS)/SP(5,L))/3.d0
  TF(1:3)=A*VR(1:3)*(B*G+D) !thermophoretic force
!
  A=ALPHAM*FNUM*G/CELL(4,NS)
  B=ALPHAT*(0.5d0*SP(5,L)*G*G-2.d0*BOLTZ*PTEMP(NS))
!
  D=0.d0
  IF ((MMRM>0).AND.(ISPR(1,L)>0)) D=ALPHAR*(PROT(N)-0.5d0*ISPR(1,L)*BOLTZ*PTEMP(NS)) !to be validated
!
  E=0.d0
  IF ((MMVM>0).AND.(ISPV(L)>0)) THEN !to be validated
    EVIB=0.d0
    FVIB=0.d0
    DO K=1,ISPV(L)
      CALL VIB_ENERGY(EVIB,IPVIB(K,N),K,L)
      !EVIB=EVIB+DFLOAT(IPVIB(K,N))*BOLTZ*SPVM(1,K,L)
      FVIB=FVIB+DEGENV*BOLTZ*SPVM(1,K,L)/(DEXP(SPVM(1,K,L)/PTEMP(NS))-1.d0)
    END DO
    E=ALPHAV*(EVIB-FVIB)
  END IF
  TH=A*(B+D+E)
!
!$omp critical(thermoph)
  CST(1:3,NS)=CST(1:3,NS)+TF(1:3)
  CST(4,NS)=CST(4,NS)+TH
!$omp end critical(thermoph)
END DO
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
 CST(0,:)=CST(0,:)+1.d0
!Now I need to create new thermophoretic variables, e.g., VARTF so that I can reset CST samples
!
!--loop over collision cells (this an alternative approach)
!DO 500 NC=1,NCCELLS
!NM=ICCELL(2,NC)      !#of molecules in collision cell nc
!
!DO 400 J=1,NM
!K=J+ICCELL(1,NC)+1   !adress of molecule j in ICREF  !MUST BE DONE AFTER INDEX TO GET UPDATED ICREF
!N=ICREF(K)           !molecule number
!L=IPSP(N)            !species of l
!
!DO EVERYTHING...
!
!400 CONTINUE
!500 CONTINUE
!
RETURN
END SUBROUTINE THERMOPHORETIC
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
!$    write(*,*) 'Available threads     ', omp_get_max_threads( )
!$    write(*,*) 'Threads in use        ', nth
!$omp end master
!$omp end parallel
!
allocate(iseed(0:nth-1))
!
do idtl=0,nth-1
call ishr(ishr3)
iseed(idtl)=ishr3         !set a seed for each thread
!$    write(*,*) 'Thread',idtl,'Seed',iseed(idtl)
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
SUBROUTINE SAMPLE_PDF
!
!--author: isebasti
!--sample pdfs
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
USE OMP_LIB
!
IMPLICIT NONE
!
REAL(KIND=8) :: A,VMPL,SC(5)
INTEGER(KIND=8) :: NS,NCI,NCF,NC,NMC,J,K,N,L,M,I,JB,KV
!
!--NS cell to be sampled
!--NBINS number of bins
!--BINS(0,I,L) sampled sum for property I and species L
!--BIN(0,I) total sampled sum for property I
!
!--set auxiliar variables
!
DBINV(1)=-3.d0 ; DBINV(2)=3.d0  !normalized velocity bin lower/upper boundaries
DBINC(1)= 0.d0 ; DBINC(2)=3.d0  !normalized speed
DBINE(1)= 0.d0 ; DBINE(2)=10.d0 !thermal energy in eV
!
DBINV(3)=(DBINV(2)-DBINV(1))/DFLOAT(NBINS) !bin intervals
DBINC(3)=(DBINC(2)-DBINC(1))/DFLOAT(NBINS)
DBINE(3)=(DBINE(2)-DBINE(1))/DFLOAT(NBINS)
!
!--sum to corresponding bin
!
NSPDF=NCELLS/2.+0.9999999d0       !user defined; must be consistent with grid
NS=NSPDF                          !sampling cell
NCI=ICELL(NS)+1                   !initial collision cell
IF(NS < NCELLS)  NCF=ICELL(NS+1)  !final collision cell
IF(NS == NCELLS) NCF=NCCELLS
IF(NCELLS == 1) NCI=1             !debuging; it can be improved
!
DO NC=NCI,NCF          !do over collision cells
  NMC=ICCELL(2,NC)     !#of molecules in collision cell NC
!
  DO J=1,NMC           !do over molecules in collision cell NC
    K=J+ICCELL(1,NC)   !adress of molecule j in ICREF
    N=ICREF(K)         !molecule number
    L=IPSP(N)          !species of molecule n
!
    VMPL=SQRT(2.D00*BOLTZ*VARSP(5,NS,L)/SP(5,L)) !most probable velocity species L
!
    DO I=1,3
      SC(I)=(PV(I,N)-VAR(I+4,NS))/VMPL                      !normalized thermal velocity
      JB=1+(SC(I)-DBINV(1))/DBINV(3)                        !bin index
      IF (JB>0.AND.JB<NBINS) BINS(JB,I,L)=BINS(JB,I,L)+1.d0 !bin accumulation per species
    END DO
!
    SC(4)=DSQRT(SC(1)*SC(1)+SC(2)*SC(2)+SC(3)*SC(3))        !normalized thermal speed
    JB=1+(SC(4)-DBINC(1))/DBINC(3)                          !bin index
    IF (JB>0.AND.JB<NBINS) BINS(JB,4,L)=BINS(JB,4,L)+1.d0   !bin accumulation per species
!
    A=0.d0
    DO I=1,3
      A=A+(PV(I,N)-VAR(I+4,NS))**2.d0
    END DO
    SC(5)=0.5d0*SP(5,L)*A/EVOLT                    !translational thermal energy (eV)
    JB=1+(SC(5)-DBINE(1))/DBINE(3)                 !bin index
    IF (JB>NBINS) JB=NBINS
    BINS(JB,5,L)=BINS(JB,5,L)+1.d0                 !bin accumulation per species
!
    IF (ISPR(1,L) > 0) THEN                        !rotational energy
      M=1+(PROT(N)/(BOLTZ*VARSP(6,NS,L)))/0.1D00   !bin index
      IF (M < 101) NDROT(L,M)=NDROT(L,M)+1         !bin accumulation per species
    END IF
!
    BIN(0,:)=BIN(0,:)+1.d0          !total counter
    BINS(0,:,L)=BINS(0,:,L)+1.d0    !species counter
  END DO
END DO
!
!--sum to corresponding bin (vibrational levels)
!
DO I=1,NSCELLS
  NS=NSVEC(I)                       !sampling cell
  NCI=ICELL(NS)+1                   !initial collision cell
  IF(NS < NCELLS)  NCF=ICELL(NS+1)  !final collision cell
  IF(NS == NCELLS) NCF=NCCELLS
  IF(NCELLS == 1) NCI=1             !debuging; it can be improved
!
  DO NC=NCI,NCF          !do over collision cells
    NMC=ICCELL(2,NC)     !#of molecules in collision cell NC
!
    DO J=1,NMC           !do over molecules in collision cell NC
      K=J+ICCELL(1,NC)   !adress of molecule j in ICREF
      N=ICREF(K)         !molecule number
      L=IPSP(N)          !species of molecule n
!
      IF (ISPV(L) > 0) THEN
        DO KV=1,ISPV(L)                                           !vibrational mode
          NDVIB(NS,KV,L,IPVIB(KV,N))=1+NDVIB(NS,KV,L,IPVIB(KV,N)) !bin accumulation per species
        END DO
      END IF
!
      NDVIB(NS,0,L,0)=1+NDVIB(NS,0,L,0)   !total counter
    END DO
!
  END DO
END DO
!
RETURN
END SUBROUTINE SAMPLE_PDF
!
!*****************************************************************************
!
SUBROUTINE UPDATE_MP
!
!--author: isebasti
!--update macroscopic properties
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: NS,K,KK,L,N,NMCR,NBC,KCELLS(2,10)     !need to set NBC and KCELLS in some commom block
REAL(KIND=8) :: A,B,C,SMCR,DOF,AVW,UU,SVDF,VDOFM,TVIBM,EVIBM,DSUM(0:12),AA,BB,SN,UBMEAN(2)
REAL(KIND=8), DIMENSION(MSP) :: TVIB,VDOF,SDOF
REAL(KIND=8), DIMENSION(MMVM,MSP) :: TV,THCOL
REAL(KIND=8), DIMENSION(NCELLS,MMVM,MSP) :: DF
!
!--set cells close to boundaries and some inner cells
!
!IF (NBC < NCELLS) NBC=NCELLS
!DO N=1,NBC
!  KCELLS(1,N)=N
!  KCELLS(2,N)=NCELLS-(N-1)
!END DO
!
!----------------------------------------------------------------------------
!
VAR=0.D00
VARSP=0.
SMCR=0
NMCR=0
VDOFM=0.
!
!--can use the lines below to consider only cells close to boundaries
!NBC=1  !consider a maximum of NBC sampling cells from each boundary
!DO NS=1,2
! DO KK=1,NBC
!  N=KCELLS(NS,KK)
!
!--consider all NCELLS sampling cells
DO N=1,NCELLS  !loop is the same as in OUTPUT subroutine
!
  A=FNUM/(CELL(4,N)*NSAMP)
  IF (IVB == 1) A=A*((XB(2)-XB(1))/(XB(2)+VELOB*0.5D00*(FTIME+TISAMP)-XB(1)))**(IFX+1)
!--check the above for non-zero XB(1)
  DSUM=0.
  NMCR=NMCR+1
  DO L=1,MSP
    DSUM(0)=DSUM(0)+CS(0,N,L)
    DSUM(1)=DSUM(1)+CS(1,N,L)
    DSUM(2)=DSUM(2)+SP(5,L)*CS(1,N,L)
    DO K=1,3
      DSUM(K+2)=DSUM(K+2)+SP(5,L)*CS(K+1,N,L)
      IF (CS(1,N,L) > 0.1D00) THEN
        VARSP(K+1,N,L)=CS(K+4,N,L)/CS(1,N,L)  !--VARSP(2,3,4 are temporarily the mean of the squares of the velocities
        VARSP(K+8,N,L)=CS(K+1,N,L)/CS(1,N,L)  !--VARSP(9,10,11 are temporarily the mean of the velocities
      END IF
    END DO
    DSUM(6)=DSUM(6)+SP(5,L)*(CS(5,N,L)+CS(6,N,L)+CS(7,N,L))
    DSUM(10)=DSUM(10)+SP(5,L)*CS(5,N,L)
    DSUM(11)=DSUM(11)+SP(5,L)*CS(6,N,L)
    DSUM(12)=DSUM(12)+SP(5,L)*CS(7,N,L)
    IF (CS(1,N,L) > 0.5D00) THEN
      DSUM(7)=DSUM(7)+CS(5,N,L)+CS(6,N,L)+CS(7,N,L)
    END IF
    IF (ISPR(1,L) > 0) THEN
      DSUM(8)=DSUM(8)+CS(8,N,L)
      DSUM(9)=DSUM(9)+CS(1,N,L)*ISPR(1,L)
    END IF
  END DO
  AVW=0.
  DO L=1,MSP
    VARSP(0,N,L)=CS(1,N,L)
    VARSP(1,N,L)=0.D00
    VARSP(6,N,L)=0.
    VARSP(7,N,L)=0.
    VARSP(8,N,L)=0.
    IF (DSUM(1) > 0.1) THEN
      VARSP(1,N,L)=CS(1,N,L)/DSUM(1)  !isebasti: deleted 100* factor
      AVW=AVW+SP(3,L)*CS(1,N,L)/DSUM(1)
      IF ((ISPR(1,L) > 0).AND.(CS(1,N,L) > 0.5)) VARSP(6,N,L)=(2.D00/BOLTZ)*CS(8,N,L)/(DFLOAT(ISPR(1,L))*CS(1,N,L))
    END IF
    VARSP(5,N,L)=0.
    DO K=1,3
      VARSP(K+1,N,L)=(SP(5,L)/BOLTZ)*(VARSP(K+1,N,L)-VARSP(K+8,N,L)**2)
      VARSP(5,N,L)=VARSP(5,N,L)+VARSP(K+1,N,L)
    END DO
    VARSP(5,N,L)=VARSP(5,N,L)/3.D00
    VARSP(8,N,L)=(3.D00*VARSP(5,N,L)+DFLOAT(ISPR(1,L))*VARSP(6,N,L))/(3.D00+DFLOAT(ISPR(1,L))) !isebasti: included according to DSMC.f90
  END DO
!
  IF (IVB == 0) VAR(1,N)=CELL(1,N)
  IF (IVB == 1) THEN
    C=(XB(2)+VELOB*FTIME-XB(1))/DFLOAT(NDIV)   !--new DDIV
    VAR(1,N)=XB(1)+(DFLOAT(N-1)+0.5)*C
  END IF
  VAR(2,N)=DSUM(0)
  IF (DSUM(1) > 0.5) THEN
    VAR(3,N)=DSUM(1)*A               !--number density Eqn. (4.28)
    VAR(4,N)=VAR(3,N)*DSUM(2)/DSUM(1) !--density Eqn. (4.29)
    VAR(5,N)=DSUM(3)/DSUM(2)          !--u velocity component Eqn. (4.30)
    VAR(6,N)=DSUM(4)/DSUM(2)          !--v velocity component
    VAR(7,N)=DSUM(5)/DSUM(2)          !--w velocity component
    UU= VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2
    IF (DSUM(1) > 1) THEN
      VAR(8,N)=(ABS(DSUM(6)-DSUM(2)*UU))/(3.D00*BOLTZ*DSUM(1))  !--translational temperature Eqn. (4.39)
      VAR(19,N)=(ABS(DSUM(10)-DSUM(2)*VAR(5,N)**2))/(BOLTZ*DSUM(1))
      VAR(20,N)=(ABS(DSUM(11)-DSUM(2)*VAR(6,N)**2))/(BOLTZ*DSUM(1))
      VAR(21,N)=(ABS(DSUM(12)-DSUM(2)*VAR(7,N)**2))/(BOLTZ*DSUM(1))
    ELSE
      VAR(8,N)=1.
      VAR(19,N)=1.
      VAR(20,N)=1.
      VAR(21,N)=1.
    END IF
!--rotational temperature
    IF (DSUM(9) > 0.01D00) THEN
      VAR(9,N)=(2.D00/BOLTZ)*DSUM(8)/DSUM(9)    !Eqn. (4.36)
    ELSE
      VAR(9,N)=0.
    END IF
    DOF=(3.D00+DSUM(9)/DSUM(1))
!--vibration temperature default
    VAR(10,N)=FTMP(1)
!--overall temperature based on translation and rotation
    VAR(11,N)=(3.*VAR(8,N)+(DSUM(9)/DSUM(1))*VAR(9,N))/DOF
!--scalar pressure (now (from V3) based on the translational temperature)
    VAR(18,N)=VAR(3,N)*BOLTZ*VAR(8,N)
!
!--Tvib calculations according to DSMC.f90
    IF (MMVM > 0) THEN
      DO L=1,MSP
        VDOF(L)=0.
        IF (ISPV(L) > 0) THEN
          DO K=1,ISPV(L)
            IF (CS(K+8,N,L) > 0.) THEN
              A=CS(K+8,N,L)/CS(1,N,L)
              TV(K,L)=SPVM(1,K,L)/DLOG(1.d0+BOLTZ*SPVM(1,K,L)/A)  !--Eqn.(4.45) - assuming SHO
              DF(N,K,L)=2.d0*A/(BOLTZ*TV(K,L)) !--Eqn. (11.28) Bird94 - general definition
            ELSE
              TV(K,L)=0.
              DF(N,K,L)=0.
            END IF
            VDOF(L)=VDOF(L)+DF(N,K,L)  !--Eqn.(4.49)
          END DO
          TVIB(L)=0.
          DO K=1,ISPV(L)
            IF (VDOF(L) > 1.D-6) THEN
              TVIB(L)=TVIB(L)+TV(K,L)*DF(N,K,L)/VDOF(L)  !--Eqn.(4.50)
            ELSE
              TVIB(L)=FVTMP(1)
            END IF
          END DO
        ELSE
          TVIB(L)=0. !TREF  !--isebasti: TREF is not defined
          VDOF(L)=0.
        END IF
        VARSP(7,N,L)=TVIB(L)
      END DO
      VDOFM=0.
      TVIBM=0.
      A=1.D00 !--isebasti: instead of 0
      DO L=1,MSP
        IF (ISPV(L) > 0) A=A+CS(1,N,L)
      END DO
      DO L=1,MSP
        IF (ISPV(L) > 0) THEN
          VDOFM=VDOFM+VDOF(L)*CS(1,N,L)/A  !--Eqn.(4.51)
          TVIBM=TVIBM+TVIB(L)*CS(1,N,L)/A  !--Eqn.(4.52)
        END IF
      END DO
      VAR(10,N)=TVIBM
    END IF
!
!--convert the species velocity components to diffusion velocities
    DO L=1,MSP
      IF (VARSP(0,N,L) > 0.5) THEN
        DO K=1,3
          VARSP(K+8,N,L)=VARSP(K+8,N,L)-VAR(K+4,N)
        END DO
      ELSE
        DO K=1,3
          VARSP(K+8,N,L)=0.D00
        END DO
      END IF
    END DO
!
!--reset the overall temperature and degrees of freedom (now including vibrational modes)
    IF (MMVM > 0) THEN
      DO L=1,MSP
        SDOF(L)=3.D00+ISPR(1,L)+VDOF(L)
        VARSP(8,N,L)=(3.*VARSP(5,N,L)+ISPR(1,L)*VARSP(6,N,L)+VDOF(L)*VARSP(7,N,L))/SDOF(L)  !species overall T
      END DO
      A=0.D00
      B=0.D00
      DO L=1,MSP
        A=A+SDOF(L)*VARSP(8,N,L)*CS(1,N,L)
        B=B+SDOF(L)*CS(1,N,L)
      END DO
      VAR(11,N)=A/B !mixture overall T
      DOF=DOF+VDOFM !--isebasti: included
    END IF
!
!--Mach number
    VAR(17,N)=DSQRT(VAR(5,N)**2+VAR(6,N)**2+VAR(7,N)**2)
    VAR(12,N)=VAR(17,N)/SQRT((DOF+2.D00)*VAR(11,N)*(DSUM(1)*BOLTZ/DSUM(2))/DOF)
!--average number of molecules in (collision) cell
    VAR(13,N)=DSUM(0)/NSAMP/DFLOAT(NCIS)
    IF (COLLS(N) > 2.) THEN
!--mean collision time
      VAR(14,N)=0.5D00*(FTIME-TISAMP)*(DSUM(1)/NSAMP)/WCOLLS(N)
!--mean free path (based on r.m.s speed with correction factor based on equilibrium)
      VAR(15,N)=0.92132D00*DSQRT(DABS(DSUM(7)/DSUM(1)-UU))*VAR(14,N)
      VAR(16,N)=CLSEP(N)/(COLLS(N)*VAR(15,N))
    ELSE
!--m.f.p set by nominal values
      VAR(14,N)=1.D10
      VAR(15,N)=1.D10/VAR(3,N)
    END IF
  ELSE
    DO L=3,19
      VAR(L,N)=0.
    END DO
  END IF
! END DO  !loop is the same as in OUTPUT
END DO
!
RETURN
END SUBROUTINE UPDATE_MP
!
!*****************************************************************************
!
SUBROUTINE UPDATE_BC
!
!--author: isebasti
!--update inflow boundary conditions
!
USE MOLECS
USE CALC
USE GEOM
USE GAS
USE OUTPUT
!
IMPLICIT NONE
!
INTEGER :: NS,K,KK,L,N,NBC,KCELLS(2,10)      !need to set NBC and KCELLS in some commom block
REAL(KIND=8) :: A,B,AA,BB,SN,UBMEAN(2)
!
!--set cells close to boundaries
!
NBC=1  !consider a maximum of NBC sampling cells from each boundary
IF (NCELLS < NBC) NBC=NCELLS
DO N=1,NBC
  KCELLS(1,N)=N
  KCELLS(2,N)=NCELLS-(N-1)
END DO
!
!--take mean velocities weighted by density
!
UBMEAN=0.
!
DO NS=1,2
  A=0.d0
  DO KK=1,NBC
    N=KCELLS(NS,KK)
    UBMEAN(NS)=UBMEAN(NS)+VAR(4,N)*VAR(5,N)
    A=A+VAR(4,N)
  END DO
  UBMEAN(NS)=UBMEAN(NS)/A !=DSUM(rho*u)/DSUM(rho)
END DO
!
!--update inlet macroscopic properties velocities
!
IF (IUBC == 1) THEN !trying to ensure boundary velocities near to VFX (assumed UFND and UFTMP are equal to freestream values)
  A=UBMEAN(1)-VFX(1)
  B=UBMEAN(2)-VFX(2)
!
  IF (DABS(A) >= 50.d0 .AND. DABS(A) < 80.d0) UVFX(1)=UVFX(1)-0.2d0*A
  IF (DABS(B) >= 50.d0 .AND. DABS(B) < 80.d0) UVFX(2)=UVFX(2)-0.2d0*B
!
  IF (DABS(A) >= 10.d0 .AND. DABS(A) < 50.d0) UVFX(1)=UVFX(1)-0.1d0*A
  IF (DABS(B) >= 10.d0 .AND. DABS(B) < 50.d0) UVFX(2)=UVFX(2)-0.1d0*B
!
  IF (DABS(A) >= 00.d0 .AND. DABS(A) < 10.d0) UVFX(1)=UVFX(1)-0.05d0*A
  IF (DABS(B) >= 00.d0 .AND. DABS(B) < 10.d0) UVFX(2)=UVFX(2)-0.05d0*B
END IF
!
IF (IUBC == 2) THEN !pressure based boundary conditions (p and T are known/fixed at boundaries); see pg 66 in my Master's thesis
  DO K=1,2
    IF(K==1) N=1
    IF(K==2) N=NCELLS
!
    A=VAR(17,N)/VAR(12,N)     !interior sound speed
    UFND(K)=FND(K)            !inlet number density
    UFTMP(K)=FTMP(K)          !inlet temperature
    B=UFND(K)*BOLTZ*UFTMP(K)  !inlet pressure
!
    IF (K == 1) UVFX(K)=VAR(5,N)+(B-VAR(18,N))/(VAR(4,N)*A)   !upstream velocity
    IF (K == 2) UVFX(K)=VAR(5,N)-(B-VAR(18,N))/(VAR(4,N)*A)   !downstream velocity
  END DO
END IF
!
!--update inlet most probable speeds
!
DO K=1,2
  DO L=1,MSP
    UVMP(L,K)=SQRT(2.D0*BOLTZ*UFTMP(K)/SP(5,L)) !most probable speed
  END DO
END DO
!
!--update the entry quantities
!
DO K=1,2
  IF (ITYPE(K) == 0) THEN
    DO L=1,MSP
      IF (K == 1) SN=UVFX(1)/UVMP(L,1)
      IF (K == 2) SN=-UVFX(2)/UVMP(L,2)
      AA=SN
      A=1.D00+ERF(AA)
      BB=DEXP(-SN**2)
      ENTR(3,L,K)=SN
      ENTR(4,L,K)=SN+SQRT(SN**2+2.D00)
      ENTR(5,L,K)=0.5D00*(1.D00+SN*(2.D00*SN-ENTR(4,L,K)))
      ENTR(6,L,K)=3.D00*UVMP(L,K)
      B=BB+SPI*SN*A
      ENTR(1,L,K)=(UFND(K)*FSP(L,K)*UVMP(L,K))*B/(FNUM*2.D00*SPI)
      ENTR(2,L,K)=0.D00
    END DO
  END IF
END DO
!
RETURN
END SUBROUTINE UPDATE_BC
!
!*****************************************************************************
!
SUBROUTINE VIB_LEVEL(B,IV,KV,KS)
!
!--for a given energy EVIB and vibrational mode KV of species KS,
!--calculate the corresponding truncated vibrational level
!
USE GAS
USE CALC
USE MOLECS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: KS,IV,KV
REAL(KIND=8) :: EVIB,A(0:IVMODEL(KS,2)),B
!
EVIB=B*1.000001d0
!
IF (IVMODEL(KS,1) == 0) THEN
  IV=EVIB/(BOLTZ*SPVM(1,KV,KS)) !SHO model
ELSE
  A(:)=1.d0
  A(:)=EVIB-VIBEN(:,KS)
  !find the first, smallest positive value in A vector
  IV=MINLOC(A,DIM=1,MASK=A.GT.0.d0)-1 !subtract 1 because index starts at zero
  IF (IV < 0)  IV=0
  IF (IV > IVMODEL(KS,2)) IV=IVMODEL(KS,2)
END IF
!
RETURN
END SUBROUTINE VIB_LEVEL
!
!*****************************************************************************
!
SUBROUTINE VIB_ENERGY(EVIB,IV,KV,KS)
!
!--for a given vibrational level IV and mode KV of species KS,
!--calculate the corresponding truncated vibrational energy EVIB
!
USE GAS
USE CALC
USE MOLECS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: KS,IV,KV
REAL(KIND=8) :: EVIB
!
IF (IVMODEL(KS,1) == 0) THEN
  EVIB=IV*BOLTZ*SPVM(1,KV,KS) !SHO model
ELSE
  EVIB=VIBEN(IV,KS) !QCT for N2 and O2 and NO
END IF
!
RETURN
END SUBROUTINE VIB_ENERGY
!
!*****************************************************************************
!
SUBROUTINE ROOTF_SECANT(I,J,L,M,A,B,X2)
!
!--use the root finding secant method
!
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: NSTEP,I,J,L,M
REAL(KIND=8) :: A,B,X1,X2,ROOTF,F1,F2,XERROR,FERROR
!
X1=A                     !1st initial guess
X2=A*1.1d0               !2nd initial guess
F1=ROOTF(I,J,L,M,B,X1)
F2=ROOTF(I,J,L,M,B,X2)
XERROR=1.d0-X1/X2        !initializing error
FERROR=1.d0-F1/F2        !initializing error
NSTEP=0
!
DO WHILE (DABS(XERROR)>1.d-6 .AND. DABS(F2)/=0.d0 .AND. DABS(FERROR)/=0.d0)
  X1=X2                       !update 1st guess
  X2=X2*(1.d0-XERROR/FERROR)  !update 2nd guess; same as X3=X2-(X2-X1)*F2/(F2-F1)
  F1=F2
  F2=ROOTF(I,J,L,M,B,X2)
  XERROR=1.d0-X1/X2
  FERROR=1.d0-F1/F2
  NSTEP=NSTEP+1
  IF (NSTEP > 10000) THEN
    WRITE(*,*) 'ROOTF_SECANT is not converging:',DABS(XERROR),DABS(F2)
    XERROR=0.d0
    X2=A                 !same as initial guess
  END IF
  !WRITE(*,*) X1,X2,DABS(XERROR),DABS(F1),DABS(F2),DABS(FERROR)
END DO
!pause
!
RETURN
END SUBROUTINE ROOTF_SECANT
!
!*****************************************************************************
!
FUNCTION ROOTF(I,J,L,M,A,X)
!
!--evaluate the functions required in the root finding methods
!
USE CALC
USE GAS
USE OMP_LIB
!
IMPLICIT NONE
!
INTEGER :: I,J,K,L,M,MAXLEV,IVIB
REAL(KIND=8) :: ROOTF,A,X,QSUM,QVIB,EVIB,AEVIB,TEMP
REAL(KIND=8) :: VDOF
!
IF (I == 1) THEN !to calculate TVIB(AEVIB)
  K=1            !assuming single vib mode
  AEVIB=A        !average vib energy
  TEMP=X         !vib temperature
  QSUM=0.d0
  QVIB=0.d0
  MAXLEV=IVMODEL(J,2)
  DO IVIB=0,MAXLEV
    CALL VIB_ENERGY(EVIB,IVIB,K,J)
    QSUM=QSUM+EVIB*DEXP(-EVIB/BOLTZ/TEMP)
    QVIB=QVIB+DEXP(-EVIB/BOLTZ/TEMP)
  END DO
  ROOTF=AEVIB-QSUM/QVIB
END IF
!
IF (I == 2) THEN !to calculate Zcorrection according to Gimelshein PoF 14 (2012)
  K=1            !assuming single vib mode
  TEMP=X         !equilibrium temperature
  VDOF=2.d0*(SPVM(1,K,J)/X)/(DEXP(SPVM(1,K,J)/X)-1.d0) !vibrational dof
  ROOTF=A-(VDOF+(5.d0-2.d0*SPM(3,L,M)))*X
END IF
!
RETURN
END FUNCTION ROOTF
!
!*****************************************************************************
!
!--this is a collection of subroutine used for multi-dimensional interpolation
!
!     ******************************************************************
      subroutine interpolate(x1a,x2a,imax,jmax,ya,x1,x2,y)
!     ******************************************************************
!     Purpose: given a position (x1,x2) bounded by data in array
!              ya(imax,jmax) (see bellow), calculate the corresponding
!              interpolated value y
!     ------------------------------------------------------------------
!              |   x1a(1)       x1a(2)      ...     x1a(imax)
!     ------------------------------------------------------------------
!     x2a(1)   |   ya(1,1)     ya(2,1)      ...     ya(imax,1)
!     x2a(2)   |   ya(1,2)     ya(2,2)      ...     ya(imax,2)
!       .      |       .          .          .          .
!       .      |       .          .          .          .
!       .      |       .          .          .          .
!     x2a(jmax)|   ya(1,jmax)   ya(2,jmax)  ...    ya(imax,jmax)
!     ------------------------------------------------------------------
      implicit none
      integer:: i,j,imax,jmax,npi,npj,ik,jk,ikf,jkf
      real*8 :: x1a(imax),x2a(jmax),ya(imax,jmax)
      real*8 :: x1,x2,y,dy
!     ------------------------------------------------------------------
!     x1 - value in x1-coord (input)
!     x2 - value in x1-coord (input)
!     y  - interpolated value for the given pair of coordinates (x1,x2)
!     ------------------------------------------------------------------
      call locate(x1a,imax,x1,i)
      call locate(x2a,jmax,x2,j)
!
      npi=2 !interpolates with a polynomial of (npi-1) degree
      npj=2 !interpolates with a polynomial of (npj-1) degree
!
      ik=min(max(i-(npi-1)/2,1),imax+1-npi)
      jk=min(max(j-(npj-1)/2,1),jmax+1-npj)
      ikf=ik+npi-1
      jkf=jk+npj-1
!
      call polin2(x1a(ik:ikf),x2a(jk:jkf),ya(ik:ikf,jk:jkf),&
                  npi,npj,x1,x2,y,dy)
!     ------------------------------------------------------------------
      return
      end subroutine interpolate
!     ******************************************************************
!
!
!     ******************************************************************
      subroutine locate(xx,n,x,j)
!     ******************************************************************
!     See Numerical Recipes
!     ------------------------------------------------------------------
      implicit none
      integer:: j,n,jl,jm,ju
      real*8 :: x,xx(n)
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end subroutine locate
!     ******************************************************************
!
!
!     ******************************************************************
      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
!     ******************************************************************
!     2D interpolation of arbitrary polinomial order
!     Given arrays x1a(1:m) and x2a(1:n) of independent variables,and an
!     m by n array of function values ya(1:m,1:n) tabulated at the grid
!     points defined by x1a,x2a; and given values x1,x2 of the indepen-
!     dent variable, this routine returns an interpolated function value
!     y with error dy
!     ------------------------------------------------------------------
      implicit none
      integer,parameter:: nmax=10,mmax=10
      integer:: m,n,j,k
      real*8 :: dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      real*8 :: ymtmp(mmax),yntmp(nmax)
!     uses polint
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      end subroutine polin2
!     ******************************************************************
!
!
!     ******************************************************************
      subroutine polint(xa,ya,n,x,y,dy)
!     ******************************************************************
!     order n polynomial interpolation using Lagrange's formula as
!     described in Numerical Recipes:
!     Given arrays xa and ya each of length n, and given a value x, this
!     routine returns a value y and an error estimate dy. If P(x) is the
!     polynomial of degree N-1 such that P(xa_i)=ya_i,i=1,...,n, then
!     the returned value is y=P(x)
!     ------------------------------------------------------------------
      implicit none
      integer,parameter:: nmax=10
      integer:: n,i,m,ns
      real*8 :: dy,x,y,xa(n),ya(n)
      real*8 :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end subroutine polint
!     ******************************************************************
!
!*****************************************************************************
! include external functions
function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Allan Macleod
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ), parameter :: xlge = 5.10D+05
  real ( kind = 8 ), parameter :: xlgst = 1.0D+30
  real ( kind = 8 ) xvalue
  real ( kind = 8 ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end
function gamain ( x, p, ifault )

!*****************************************************************************80
!
!! GAMAIN computes the incomplete gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!    PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    G Bhattacharjee
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    G Bhattacharjee,
!    Algorithm AS 32:
!    The Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 19, Number 3, 1970, pages 285-287.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma routine.
!
!    Output, real ( kind = 8 ) GAMAIN, the value of the incomplete
!    gamma ratio.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: acu = 1.0D-08
  real ( kind = 8 ) alngam
  real ( kind = 8 ) an
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) dif
  real ( kind = 8 ) factor
  real ( kind = 8 ) g
  real ( kind = 8 ) gamain
  real ( kind = 8 ) gin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), parameter :: oflo = 1.0D+37
  real ( kind = 8 ) p
  real ( kind = 8 ) pn(6)
  real ( kind = 8 ) rn
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: uflo = 1.0D-37
  real ( kind = 8 ) x
!
!  Check the input.
!
  if ( p <= 0.0D+00 ) then
    ifault = 1
    gamain = 0.0D+00
    return
  end if

  if ( x < 0.0D+00 ) then
    ifault = 2
    gamain = 0.0D+00
    return
  end if

  if ( x == 0.0D+00 ) then
    ifault = 0
    gamain = 0.0D+00
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0D+00
    return
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0D+00
    return
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    gin = 1.0D+00
    term = 1.0D+00
    rn = p

    do

      rn = rn + 1.0D+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0D+00 - p
  b = a + x + 1.0D+00
  term = 0.0D+00

  pn(1) = 1.0D+00
  pn(2) = x
  pn(3) = x + 1.0D+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    a = a + 1.0D+00
    b = b + 2.0D+00
    term = term + 1.0D+00
    an = a * term
    do i = 1, 2
      pn(i+4) = b * pn(i+2) - an * pn(i)
    end do

    if ( pn(6) /= 0.0D+00 ) then

      rn = pn(5) / pn(6)
      dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
      if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
        if ( dif <= acu * rn ) then
          gamain = 1.0D+00 - factor * gin
          exit
        end if

      end if

      gin = rn

    end if

    do i = 1, 4
      pn(i) = pn(i+2)
    end do

    if ( oflo <= abs ( pn(5) ) ) then

      do i = 1, 4
        pn(i) = pn(i) / oflo
      end do

    end if

  end do

  return
end
subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function P(A,X) is defined as:
!
!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= PN(A,X) <= 1
!
!    and
!
!      PN(A,INFINITY) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
!
!  Modified:
!
!    20 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7382350532339351D+00, &
    0.9083579897300343D+00, &
    0.9886559833621947D+00, &
    0.3014646416966613D+00, &
    0.7793286380801532D+00, &
    0.9918490284064973D+00, &
    0.9516258196404043D-01, &
    0.6321205588285577D+00, &
    0.9932620530009145D+00, &
    0.7205974576054322D-01, &
    0.5891809618706485D+00, &
    0.9915368159845525D+00, &
    0.1018582711118352D-01, &
    0.4421745996289254D+00, &
    0.9927049442755639D+00, &
    0.4202103819530612D-01, &
    0.9796589705830716D+00, &
    0.9226039842296429D+00, &
    0.4470785799755852D+00, &
    0.7444549220718699D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

! subroutine to search index ii and ij in x meets
!       x(ii) < x0 < x(ij)
! if x0 equals to some value in x, ii=ij=index of the element
      subroutine binary_search(ii,ij,x0,x,lx,rx,ierror)
        implicit none
        integer :: ii, ij
        integer,intent(in) ::  rx,lx
        double precision,intent(in) :: x0,x(lx:rx)
        logical :: ierror
        integer :: L,R,M

        ierror = .false.
        L = lx
        R = rx
        M = (L+R)/2

        if ( (x0 .lt. x(L)) .or. (x0 .gt. x(R)))then
          write(6,*) 'No correct phase is found'
          ierror = .true.
          return
        elseif (x0 .eq. x(L)) then
          ii = L
          ij = L
          return
        elseif (x0 .eq. x(R)) then
          ii = R
          ij = R
          return
        endif

        do while ( (x(L) .le. x0) .and. (x(R) .ge. x0) .and. (R-L)>1 )
          if (x0 .eq. x(M)) then
            ii = M
            ij = M
            R = M
            L =M
            exit
          elseif (x0 .gt. x(M))then
            L = M
            M = (L+R)/2
          elseif (x0 .le. x(M)) then
            R = M
            M = (L+R)/2
          endif
        enddo
        if ( x(R) .eq. x0) then
          ii = R
          ij = R
        elseif ( x(L) .eq. x0) then
          ii = L
          ij = L
        else
          ii = L
          ij = R
        end if

      end subroutine binary_search

! the following one is very similar to the previous one
! the difference is that, this one also handles list like the folowing
! x=      [0,0,0,2,5,5,5,8]
! index = [0,1,2,3,4,5,6,7]
! x0=1 -> ii=2, ij=3
! x0=5 -> ii=ij=4
! noted that x0 must be larger than 0
      subroutine binary_search_int(ii,ij,x0,x,lx,rx,ierror)
        implicit none
        integer :: ii, ij
        integer,intent(in) ::  rx,lx
        integer,intent(in) :: x0,x(lx:rx)
        logical :: ierror
        integer :: L,R,M

        ierror = .false.
        L = lx
        R = rx
        M = (L+R)/2

        if ( (x0 .lt. x(L)) .or. (x0 .gt. x(R)))then
          write(6,*) 'No correct phase is found'
          ierror = .true.
          return
        elseif (x0 .eq. x(L)) then
          ii = L
          ij = L
          return
        elseif (x0 .eq. x(R)) then
          ii = R
          ij = R
          return
        endif

        do while ( (x(L) .le. x0) .and. (x(R) .ge. x0) .and. (R-L)>1 )
          if (x0 .eq. x(M)) then
            ii = M
            ij = M
            R = M
            L =M
            exit
          elseif (x0 .gt. x(M))then
            L = M
            M = (L+R)/2
          elseif (x0 .le. x(M)) then
            R = M
            M = (L+R)/2
          endif
        enddo
        if ( x(R) .eq. x0) then
          ii = R
          ij = R
        elseif ( x(L) .eq. x0) then
          ii = L
          ij = L
        else
          ii = L
          ij = R
        end if

      end subroutine binary_search_int

! -----------------------------------------------------------------------------------------
! Brent algorithm for solving equation
      double precision function brent(ax,bx,f,tol,fa0,fb0)
!xxx      double precision function zeroin(ax,bx,f,tol)
!
! This is a slightly modified version of the zeroin.f source code from netlib,
! which implements Brent's zero finding algorithm
! We have added two arguments so that if the function values at the end
! points are already known they won't be recalculated. If these are set
! to 0d0 initially they will be calculated within the routine.
! We have also determined the machine precision using a fortran 90
! intrinsic rather than the original use of d1mach, and we have
! reformatted for free format
!
      double precision, intent(in) :: ax, bx, fa0, fb0
      double precision :: tol
      double precision :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision, external :: f
!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
!xxx      double precision  dabs, d1mach
!xxx   10 eps = d1mach(4)

!     add by Han, initialize brent
      brent = -1.0d0
      eps=epsilon(1d0)
      tol1 = eps+1.0d0
!
      a=ax
      b=bx
      if(fa0.ne.0d0)then
        fa=fa0
      else
        fa=f(a)
      endif
      if(fb0.ne.0d0)then
        fb=fb0
      else
        fb=f(b)
      endif
!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(6,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' brent is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.  &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
!xxx  150 zeroin=b
  150 brent=b
      end function brent


! -----------------------------------------------------------------------------------------
! All the following is borrowed from
!      https://people.sc.fsu.edu/~jburkardt/f_src/chebyshev1_rule/chebyshev1_rule.html
! Solve the nodes ang weights for gauss-quad integration

! program main
    ! implicit none
    ! integer :: nt, kind
    ! real(8) :: alpha, beta,a,b
    ! real(8),allocatable :: t(:), wts(:)
    ! write(*,*) "NT = "
    ! read(*,*) nt
    ! write(*,*) "kind = "
    ! write(*,*) "1, Legendre,             (a,b)       1.0"
    ! write(*,*) "2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)"
    ! write(*,*) "3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha"
    ! write(*,*) "4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta"
    ! write(*,*) "5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))"
    ! write(*,*) "6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)"
    ! write(*,*) "7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha"
    ! write(*,*) "8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta"
    ! write(*,*) "9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)"
    ! read(*,*) kind
    ! write(*,*) "a b"
    ! read(*,*) a,b
    ! allocate(t(nt), wts(nt))
    ! call cgqf ( nt, kind, alpha, beta, a, b, t, wts )
    ! write(*,*) "Nodes:"
    ! write(*,*) t
    ! write(*,*) "weight:"
    ! write(*,*) wts
    ! deallocate(t, wts)
    !
!
! end
subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu

  call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end
subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

!*****************************************************************************80
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints, or
!    other parameters.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with
!  valid A and B.
!
  allocate ( mlt(1:nt) )

  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  do i = 1, nt
    ndx(i) = i
  end do

  call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end
subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) aba
  real ( kind = 8 ) abi
  real ( kind = 8 ) abj
  real ( kind = 8 ) abti
  real ( kind = 8 ) aj(m)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) apone
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) zemu

  temp = epsilon ( temp )

  call parchk ( kind, 2 * m - 1, alpha, beta )

  temp2 = 0.5D+00

  if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
    stop
  end if

  if ( kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 )**2 &
      / r8_gamma ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi**2
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 5 ) then

    zemu = r8_gamma ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 6 ) then

    zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti**2 ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 9 ) then

    zemu = pi / 2.0D+00
    aj(1:m) = 0.0D+00
    bj(1:m) = 0.5D+00

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine parchk ( kind, m, alpha, beta )

!*****************************************************************************80
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer ( kind = 4 ) M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) m
  real ( kind = 8 ) tmp

  if ( kind <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND <= 0.'
    stop
  end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
  if ( 3 <= kind .and. kind <= 8 .and. alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  3 <= KIND and ALPHA <= -1.'
    stop
  end if
!
!  Check BETA for Jacobi.
!
  if ( kind == 4 .and. kind <= 8 .and. beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND == 4 and BETA <= -1.0.'
    stop
  end if
!
!  Check ALPHA and BETA for rational.
!
  if ( kind == 8 ) then
    tmp = alpha + beta + m + 1.0D+00
    if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  KIND == 8 but condition on ALPHA and BETA fails.'
      stop
    end if
  end if

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) swts(nwts)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)

  temp = epsilon ( temp )

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0'
      stop
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0.'
      stop
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0.'
      stop
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end
subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************80
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu
!
!  Exit if the zero-th moment is not positive.
!
  if ( zemu <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGQF - Fatal error!'
    write ( *, '(a)' ) '  ZEMU <= 0.'
    stop
  end if
!
!  Set up vectors for IMTQLX.
!
  t(1:nt) = aj(1:nt)

  wts(1) = sqrt ( zemu )
  wts(2:nt) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

  return
end
