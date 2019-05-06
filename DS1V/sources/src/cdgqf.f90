
! -----------------------------------------------------------------------------------------
! Brent algorithm for solving equation
      double precision function brent(ax,bx,f,tol,fa0,fb0)
!xxx      double precision function zeroin(ax,bx,f,tol)
!
! This is a slightly modified version of the zeroin.f source code from netlib,
! which implements Brent's zero finding algorithm
! We have added two arguments so that if the function values at the end
! points are already known they won't be recalculated. If these are set
! to 0d0 initially they will be calculated within the routine.
! We have also determined the machine precision using a fortran 90
! intrinsic rather than the original use of d1mach, and we have
! reformatted for free format
!
      double precision, intent(in) :: ax, bx, fa0, fb0
      double precision :: tol
      double precision :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision, external :: f
!
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
!xxx      double precision  dabs, d1mach
!xxx   10 eps = d1mach(4)

!     add by Han, initialize brent
      brent = -1.0d0
      eps=epsilon(1d0)
      tol1 = eps+1.0d0
!
      a=ax
      b=bx
      if(fa0.ne.0d0)then
        fa=fa0
      else
        fa=f(a)
      endif
      if(fb0.ne.0d0)then
        fb=fb0
      else
        fb=f(b)
      endif
!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
         write(6,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,', &
                   ' brent is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.  &
      dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
!xxx  150 zeroin=b
  150 brent=b
      end function brent


! -----------------------------------------------------------------------------------------
! All the following is borrowed from
!      https://people.sc.fsu.edu/~jburkardt/f_src/chebyshev1_rule/chebyshev1_rule.html
! Solve the nodes ang weights for gauss-quad integration

! program main
    ! implicit none
    ! integer :: nt, kind
    ! real(8) :: alpha, beta,a,b
    ! real(8),allocatable :: t(:), wts(:)
    ! write(*,*) "NT = "
    ! read(*,*) nt
    ! write(*,*) "kind = "
    ! write(*,*) "1, Legendre,             (a,b)       1.0"
    ! write(*,*) "2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)"
    ! write(*,*) "3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha"
    ! write(*,*) "4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta"
    ! write(*,*) "5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))"
    ! write(*,*) "6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)"
    ! write(*,*) "7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha"
    ! write(*,*) "8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta"
    ! write(*,*) "9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)"
    ! read(*,*) kind
    ! write(*,*) "a b"
    ! read(*,*) a,b
    ! allocate(t(nt), wts(nt))
    ! call cgqf ( nt, kind, alpha, beta, a, b, t, wts )
    ! write(*,*) "Nodes:"
    ! write(*,*) t
    ! write(*,*) "weight:"
    ! write(*,*) wts
    ! deallocate(t, wts)
    !
!
! end
subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu

  call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end
