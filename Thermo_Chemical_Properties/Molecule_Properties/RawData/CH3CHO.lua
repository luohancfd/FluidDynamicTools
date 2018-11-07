-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH3CHO = {}
CH3CHO.M = {
   value = 44.052560e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH3CHO.atomic_constituents = {C=2,H=4,O=1}
CH3CHO.charge = 0
CH3CHO.gamma = {
   value = 1.1762e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH3CHO.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.729459500e00, -3.193285800e-03, 4.753492100e-05, -5.745861100e-08, 2.193111200e-11, -2.157287800e04, 4.103015900e00, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.404110800e00, 1.172305900e-02, -4.226313700e-06, 6.837245100e-10, -4.098486300e-14, -2.259312200e04, -3.480791700e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
