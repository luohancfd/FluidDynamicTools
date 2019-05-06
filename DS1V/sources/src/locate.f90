!     ******************************************************************
!
!
!     ******************************************************************
      subroutine locate(xx,n,x,j)
!     ******************************************************************
!     See Numerical Recipes
!     ------------------------------------------------------------------
      implicit none
      integer:: j,n,jl,jm,ju
      real*8 :: x,xx(n)
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end subroutine locate
