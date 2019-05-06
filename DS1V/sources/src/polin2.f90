!     ******************************************************************
!
!
!     ******************************************************************
      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
!     ******************************************************************
!     2D interpolation of arbitrary polinomial order
!     Given arrays x1a(1:m) and x2a(1:n) of independent variables,and an
!     m by n array of function values ya(1:m,1:n) tabulated at the grid
!     points defined by x1a,x2a; and given values x1,x2 of the indepen-
!     dent variable, this routine returns an interpolated function value
!     y with error dy
!     ------------------------------------------------------------------
      implicit none
      integer,parameter:: nmax=10,mmax=10
      integer:: m,n,j,k
      real*8 :: dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      real*8 :: ymtmp(mmax),yntmp(nmax)
!     uses polint
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      end subroutine polin2
