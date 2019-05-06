subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
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
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) aba
  real ( kind = 8 ) abi
  real ( kind = 8 ) abj
  real ( kind = 8 ) abti
  real ( kind = 8 ) aj(m)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) apone
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) zemu

  temp = epsilon ( temp )

  call parchk ( kind, 2 * m - 1, alpha, beta )

  temp2 = 0.5D+00

  if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
    stop
  end if

  if ( kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 )**2 &
      / r8_gamma ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi**2
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 5 ) then

    zemu = r8_gamma ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 6 ) then

    zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti**2 ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 9 ) then

    zemu = pi / 2.0D+00
    aj(1:m) = 0.0D+00
    bj(1:m) = 0.5D+00

  end if

  return
end
