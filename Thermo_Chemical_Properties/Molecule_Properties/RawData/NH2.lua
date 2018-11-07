-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

NH2 = {}
NH2.M = {
   value = 16.022580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
NH2.atomic_constituents = {N=1,H=2}
NH2.charge = 0
NH2.gamma = {
   value = 1.3254e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
NH2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.204002900e00, -2.106138500e-03, 7.106834800e-06, -5.611519700e-09, 1.644071700e-12, 2.188591000e04, -1.418424800e-01, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.834742100e00, 3.207308200e-03, -9.339080400e-07, 1.370295300e-10, -7.920614400e-15, 2.217195700e04, 6.520416300e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
