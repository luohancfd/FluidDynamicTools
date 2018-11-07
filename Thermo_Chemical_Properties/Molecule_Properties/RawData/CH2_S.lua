-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH2_S = {}
CH2_S.M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH2_S.atomic_constituents = {C=1,H=2}
CH2_S.charge = 0
CH2_S.gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH2_S.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.198604110e00, -2.366614190e-03, 8.232962200e-06, -6.688159810e-09, 1.943147370e-12, 5.049681630e04, -7.691189670e-01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.292038420e00, 4.655886370e-03, -2.011919470e-06, 4.179060000e-10, -3.397163650e-14, 5.092599970e04, 8.626501690e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
