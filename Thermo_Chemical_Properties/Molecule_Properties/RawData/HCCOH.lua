-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HCCOH = {}
HCCOH.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HCCOH.atomic_constituents = {C=2,H=2,O=1}
HCCOH.charge = 0
HCCOH.gamma = {
   value = 1.1656e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HCCOH.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 1.242373300e00, 3.107220100e-02, -5.086686400e-05, 4.313713100e-08, -1.401459400e-11, 8.031614300e03, 1.387431900e01, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.923829100e00, 6.792360000e-03, -2.565856400e-06, 4.498784100e-10, -2.994010100e-14, 7.264626000e03, -7.601774200e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
