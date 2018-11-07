-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH3O = {}
CH3O.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH3O.atomic_constituents = {C=1,H=3,O=1}
CH3O.charge = 0
CH3O.gamma = {
   value = 1.2802e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH3O.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.106204000e+00, 7.216595000e-03, 5.338472000e-06, -7.377636000e-09, 2.075610000e-12, 9.786011000e+02, 1.315217700e+01, }
   },
   { T_low  = 1000.0,
     T_high = 3000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.770799000e+00, 7.871497000e-03, -2.656384000e-06, 3.944431000e-10, -2.112616000e-14, 1.278325200e+02, 2.929575000e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
