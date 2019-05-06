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
