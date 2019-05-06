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
OPEN (3,FILE='DS1VD.TXT')
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

READ(4,*) nonVHS
WRITE(3, *) ' noVHS = ',nonVHS
IF (nonVHS == 3 .and. IRM .ne. 250 .and. IRM .ne. 150 .and. IRM .ne. 160 .and. IRM .ne. 0) THEN
  WRITE(3,*) ' nonVHS = 3 should be used with IRM == 250/150/160'
  WRITE(*,*) ' nonVHS = 3 should be used with IRM == 250/150/160'
  STOP
END IF

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
IF ( IMF .EQ. 0 ) IMFS = 0

IF (nonVHS .ne. 0 .and. IMF == 0) THEN
  WRITE(*,*) ' nonVHS has only been tested with MF-DSMC model'
  WRITE(3,*) ' nonVHS has only been tested with MF-DSMC model'
  STOP
END IF

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

IF (QCTMODEL .ne. 2 .and. QCTMODEL .ne. 1 .AND. IMF .NE. 0)THEN
  WRITE(3,*) ' WARNING: MF MODEL CAN ONLY BE USED WIITH QCTMODEL=1/2'
  WRITE(*,*) ' WARNING: MF MODEL CAN ONLY BE USED WIITH QCTMODEL=1/2'
  stop
ENDIF
!==============================================================
!
CLOSE (3)
CLOSE (4)
!
RETURN
END SUBROUTINE READ_DATA
