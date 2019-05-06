!     ******************************************************************
!
!
!     ******************************************************************
      subroutine polint(xa,ya,n,x,y,dy)
!     ******************************************************************
!     order n polynomial interpolation using Lagrange's formula as
!     described in Numerical Recipes:
!     Given arrays xa and ya each of length n, and given a value x, this
!     routine returns a value y and an error estimate dy. If P(x) is the
!     polynomial of degree N-1 such that P(xa_i)=ya_i,i=1,...,n, then
!     the returned value is y=P(x)
!     ------------------------------------------------------------------
      implicit none
      integer,parameter:: nmax=10
      integer:: n,i,m,ns
      real*8 :: dy,x,y,xa(n),ya(n)
      real*8 :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end subroutine polint
