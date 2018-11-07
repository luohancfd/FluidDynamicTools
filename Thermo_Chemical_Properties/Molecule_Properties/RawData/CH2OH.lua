-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH2OH = {}
CH2OH.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH2OH.atomic_constituents = {C=1,H=3,O=1}
CH2OH.charge = 0
CH2OH.gamma = {
   value = 1.2070e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH2OH.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.863889180e00, 5.596723040e-03, 5.932717910e-06, -1.045320120e-08, 4.369672780e-12, -3.193913670e03, 5.473022430e00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.692665690e00, 8.645767970e-03, -3.751011200e-06, 7.872346360e-10, -6.485542010e-14, -3.242506270e03, 5.810432150e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
