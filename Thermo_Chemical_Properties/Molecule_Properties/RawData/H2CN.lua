-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

H2CN = {}
H2CN.M = {
   value = 28.033280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
H2CN.atomic_constituents = {C=1,H=2,N=1}
H2CN.charge = 0
H2CN.gamma = {
   value = 1.2769e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
H2CN.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.851661000e00, 5.695233100e-03, 1.071140000e-06, -1.622612000e-09, -2.351108100e-13, 2.863782000e04, 8.992751100e00, }
   },
   { T_low  = 1000.0,
     T_high = 4000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.209703000e00, 2.969291100e-03, -2.855589100e-07, -1.635550000e-10, 3.043258900e-14, 2.767710900e04, -4.444478000e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
