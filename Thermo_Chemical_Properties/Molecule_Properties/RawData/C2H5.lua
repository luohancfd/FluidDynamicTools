-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H5 = {}
C2H5.M = {
   value = 29.061100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H5.atomic_constituents = {C=2,H=5}
C2H5.charge = 0
C2H5.gamma = {
   value = 1.1963e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H5.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.306465680e+00, -4.186588920e-03, 4.971428070e-05, -5.991266060e-08, 2.305090040e-11, 1.284162650e+04, 4.707209240e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 1.954656420e+00, 1.739727220e-02, -7.982066680e-06, 1.752176890e-09, -1.496415760e-13, 1.285752000e+04, 1.346243430e+01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
