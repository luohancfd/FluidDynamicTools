subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
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
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) swts(nwts)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)

  temp = epsilon ( temp )

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0'
      stop
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0.'
      stop
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0.'
      stop
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end
