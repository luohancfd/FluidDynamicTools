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
