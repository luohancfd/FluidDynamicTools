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
