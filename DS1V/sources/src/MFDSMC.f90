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
      OPEN(3,FILE="DS1VD.TXT", ACCESS="APPEND")
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
    ELSE IF (IMF == 3 .or. IMF == 4) THEN
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

  subroutine MF_INIT_SAMPLE
    implicit none
    NMFEV0 = 0.; NMFER0 = 0.; NMFET0 = 0.
    NMFEV = 0.; NMFER = 0.; NMFET = 0.
    NMFEVR = 0.; NMFERR = 0.; NMFETR = 0.
    NMFVT0 = 0.; NMFVT = 0.; NMFVTR = 0.
    !NCANGLE = 0; NCRANGLE = 0.
  end subroutine MF_INIT_SAMPLE

!
!--IMFS = 0, not sample anything related to MF model =1 sample
!--NMFANG = number of angle to sample for MF model, <8 for atom-diatom, >= 8 for diatom-diatom
!--NMFpair number of molecular pairs to use MF dissociation model
!--IMFpair store the index of different molecular pair
!--MFRMASS store reduced mass

END MODULE MFDSMC
