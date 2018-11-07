-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH3OH = {}
CH3OH.M = {
   value = 32.041860e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH3OH.atomic_constituents = {C=1,H=4,O=1}
CH3OH.charge = 0
CH3OH.gamma = {
   value = 1.2320e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH3OH.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.715395820e00, -1.523091290e-02, 6.524411550e-05, -7.108068890e-08, 2.613526980e-11, -2.564276560e04, -1.504098230e00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 1.789707910e00, 1.409382920e-02, -6.365008350e-06, 1.381710850e-09, -1.170602200e-13, -2.537487470e04, 1.450236230e01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
