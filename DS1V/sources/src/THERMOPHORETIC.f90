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
