-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HCCO = {}
HCCO.M = {
   value = 41.028740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HCCO.atomic_constituents = {C=2,H=1,O=1}
HCCO.charge = 0
HCCO.gamma = {
   value = 1.2067e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HCCO.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.251721400e00, 1.765502100e-02, -2.372910100e-05, 1.727575900e-08, -5.066481100e-12, 2.005944900e04, 1.249041700e01, }
   },
   { T_low  = 1000.0,
     T_high = 4000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.628205800e00, 4.085340100e-03, -1.593454700e-06, 2.862605200e-10, -1.940783200e-14, 1.932721500e04, -3.930259500e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
