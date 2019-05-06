!
!*****************************************************************************
!
SUBROUTINE ISHR(ISHR3L)
!
!--author: isebasti
!--calculations required by Ziggurate Method; initialize different seeds
!
USE CALC
!
IMPLICIT NONE
!
INTEGER :: jsrl,jsr_inputl,ishr3l  !local variables
!
jsrl=ijsr
jsr_inputl = jsrl
jsrl = ieor(jsrl,ishft(jsrl,13))
jsrl = ieor(jsrl,ishft(jsrl,-17))
jsrl = ieor(jsrl,ishft(jsrl,  5))
ishr3l = jsr_inputl+jsrl
ijsr=jsrl
!
RETURN
END SUBROUTINE ISHR
