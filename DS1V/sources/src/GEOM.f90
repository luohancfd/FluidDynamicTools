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
