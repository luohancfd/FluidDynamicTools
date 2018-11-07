-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH2CHO = {}
CH2CHO.M = {
   value = 43.044620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH2CHO.atomic_constituents = {C=2,H=3,O=1}
CH2CHO.charge = 0
CH2CHO.gamma = {
   value = 1.1776e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH2CHO.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.409062000e00, 1.073857400e-02, 1.891492000e-06, -7.158583000e-09, 2.867385000e-12, 1.521476600e03, 9.558290000e00, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.975670000e00, 8.130591000e-03, -2.743624000e-06, 4.070304000e-10, -2.176017000e-14, 4.903218000e02, -5.045251000e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
