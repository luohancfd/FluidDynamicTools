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

!-- Han's trick for nonVHS
IF (nonVHS .eq. 1) CCELL(4,:) = CCELL(4,:)*1.1D0
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
