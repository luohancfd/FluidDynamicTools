subroutine gamma_inc_values ( n_data, a, x, fx )

!*****************************************************************************80
!
!! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
!
!  Discussion:
!
!    The (normalized) incomplete Gamma function P(A,X) is defined as:
!
!      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
!
!    With this definition, for all A and X,
!
!      0 <= PN(A,X) <= 1
!
!    and
!
!      PN(A,INFINITY) = 1.0
!
!    In Mathematica, the function can be evaluated by:
!
!      1 - GammaRegularized[A,X]
!
!  Modified:
!
!    20 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, the parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.10D+00, &
    0.10D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+01, &
    0.10D+01, &
    0.10D+01, &
    0.11D+01, &
    0.11D+01, &
    0.11D+01, &
    0.20D+01, &
    0.20D+01, &
    0.20D+01, &
    0.60D+01, &
    0.60D+01, &
    0.11D+02, &
    0.26D+02, &
    0.41D+02 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7382350532339351D+00, &
    0.9083579897300343D+00, &
    0.9886559833621947D+00, &
    0.3014646416966613D+00, &
    0.7793286380801532D+00, &
    0.9918490284064973D+00, &
    0.9516258196404043D-01, &
    0.6321205588285577D+00, &
    0.9932620530009145D+00, &
    0.7205974576054322D-01, &
    0.5891809618706485D+00, &
    0.9915368159845525D+00, &
    0.1018582711118352D-01, &
    0.4421745996289254D+00, &
    0.9927049442755639D+00, &
    0.4202103819530612D-01, &
    0.9796589705830716D+00, &
    0.9226039842296429D+00, &
    0.4470785799755852D+00, &
    0.7444549220718699D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.30D-01, &
    0.30D+00, &
    0.15D+01, &
    0.75D-01, &
    0.75D+00, &
    0.35D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.10D+00, &
    0.10D+01, &
    0.50D+01, &
    0.15D+00, &
    0.15D+01, &
    0.70D+01, &
    0.25D+01, &
    0.12D+02, &
    0.16D+02, &
    0.25D+02, &
    0.45D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
