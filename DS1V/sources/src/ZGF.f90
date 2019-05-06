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
