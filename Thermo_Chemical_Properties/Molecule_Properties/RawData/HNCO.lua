-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HNCO = {}
HNCO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HNCO.atomic_constituents = {C=1,H=1,N=1,O=1}
HNCO.charge = 0
HNCO.gamma = {
   value = 1.2173e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HNCO.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1478.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.630963170e00, 7.302823570e-03, -2.280500030e-06, -6.612712980e-10, 3.622357520e-13, -1.558736360e04, 6.194577270e00, }
   },
   { T_low  = 1478.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 6.223951340e00, 3.178640040e-03, -1.093787550e-06, 1.707351630e-10, -9.950219550e-15, -1.665993440e04, -8.382247410e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
