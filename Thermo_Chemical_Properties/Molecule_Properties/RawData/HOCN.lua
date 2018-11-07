-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HOCN = {}
HOCN.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HOCN.atomic_constituents = {C=1,H=1,N=1,O=1}
HOCN.charge = 0
HOCN.gamma = {
   value = 1.2185e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HOCN.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1368.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.786049520e00, 6.886679220e-03, -3.214878640e-06, 5.171957670e-10, 1.193607880e-14, -2.826984000e03, 5.632921620e00, }
   },
   { T_low  = 1368.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.897848850e00, 3.167893930e-03, -1.118010640e-06, 1.772431440e-10, -1.043391770e-14, -3.706533310e03, -6.181678250e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
