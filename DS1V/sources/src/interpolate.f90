!
!*****************************************************************************
!
!--this is a collection of subroutine used for multi-dimensional interpolation
!
!     ******************************************************************
      subroutine interpolate(x1a,x2a,imax,jmax,ya,x1,x2,y)
!     ******************************************************************
!     Purpose: given a position (x1,x2) bounded by data in array
!              ya(imax,jmax) (see bellow), calculate the corresponding
!              interpolated value y
!     ------------------------------------------------------------------
!              |   x1a(1)       x1a(2)      ...     x1a(imax)
!     ------------------------------------------------------------------
!     x2a(1)   |   ya(1,1)     ya(2,1)      ...     ya(imax,1)
!     x2a(2)   |   ya(1,2)     ya(2,2)      ...     ya(imax,2)
!       .      |       .          .          .          .
!       .      |       .          .          .          .
!       .      |       .          .          .          .
!     x2a(jmax)|   ya(1,jmax)   ya(2,jmax)  ...    ya(imax,jmax)
!     ------------------------------------------------------------------
      implicit none
      integer:: i,j,imax,jmax,npi,npj,ik,jk,ikf,jkf
      real*8 :: x1a(imax),x2a(jmax),ya(imax,jmax)
      real*8 :: x1,x2,y,dy
!     ------------------------------------------------------------------
!     x1 - value in x1-coord (input)
!     x2 - value in x1-coord (input)
!     y  - interpolated value for the given pair of coordinates (x1,x2)
!     ------------------------------------------------------------------
      call locate(x1a,imax,x1,i)
      call locate(x2a,jmax,x2,j)
!
      npi=2 !interpolates with a polynomial of (npi-1) degree
      npj=2 !interpolates with a polynomial of (npj-1) degree
!
      ik=min(max(i-(npi-1)/2,1),imax+1-npi)
      jk=min(max(j-(npj-1)/2,1),jmax+1-npj)
      ikf=ik+npi-1
      jkf=jk+npj-1
!
      call polin2(x1a(ik:ikf),x2a(jk:jkf),ya(ik:ikf,jk:jkf),&
                  npi,npj,x1,x2,y,dy)
!     ------------------------------------------------------------------
      return
      end subroutine interpolate
