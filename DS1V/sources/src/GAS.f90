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
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: IQCTVT(:,:)

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
!--INONVHS(L,M): 0: VHS model between specie L and M
!                1: Special QCT based model for N2+O
!                2: Collision model based on exponential potential
!--IQCTVT(L,M): .false.  ME-QCT-VT model doesn't exist
!               .true.   ME-QCT-VT model exists
!--IVFP   0 set y velocity based on VFY, 1 set a Couette Flow type profile
END MODULE GAS
