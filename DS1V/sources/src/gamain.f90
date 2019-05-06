function gamain ( x, p, ifault )

!*****************************************************************************80
!
!! GAMAIN computes the incomplete gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!    PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    G Bhattacharjee
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    G Bhattacharjee,
!    Algorithm AS 32:
!    The Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 19, Number 3, 1970, pages 285-287.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma routine.
!
!    Output, real ( kind = 8 ) GAMAIN, the value of the incomplete
!    gamma ratio.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: acu = 1.0D-08
  real ( kind = 8 ) alngam
  real ( kind = 8 ) an
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) dif
  real ( kind = 8 ) factor
  real ( kind = 8 ) g
  real ( kind = 8 ) gamain
  real ( kind = 8 ) gin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), parameter :: oflo = 1.0D+37
  real ( kind = 8 ) p
  real ( kind = 8 ) pn(6)
  real ( kind = 8 ) rn
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: uflo = 1.0D-37
  real ( kind = 8 ) x
!
!  Check the input.
!
  if ( p <= 0.0D+00 ) then
    ifault = 1
    gamain = 0.0D+00
    return
  end if

  if ( x < 0.0D+00 ) then
    ifault = 2
    gamain = 0.0D+00
    return
  end if

  if ( x == 0.0D+00 ) then
    ifault = 0
    gamain = 0.0D+00
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0D+00
    return
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0D+00
    return
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    gin = 1.0D+00
    term = 1.0D+00
    rn = p

    do

      rn = rn + 1.0D+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0D+00 - p
  b = a + x + 1.0D+00
  term = 0.0D+00

  pn(1) = 1.0D+00
  pn(2) = x
  pn(3) = x + 1.0D+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    a = a + 1.0D+00
    b = b + 2.0D+00
    term = term + 1.0D+00
    an = a * term
    do i = 1, 2
      pn(i+4) = b * pn(i+2) - an * pn(i)
    end do

    if ( pn(6) /= 0.0D+00 ) then

      rn = pn(5) / pn(6)
      dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
      if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
        if ( dif <= acu * rn ) then
          gamain = 1.0D+00 - factor * gin
          exit
        end if

      end if

      gin = rn

    end if

    do i = 1, 4
      pn(i) = pn(i+2)
    end do

    if ( oflo <= abs ( pn(5) ) ) then

      do i = 1, 4
        pn(i) = pn(i) / oflo
      end do

    end if

  end do

  return
end
