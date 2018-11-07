-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HCNO = {}
HCNO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HCNO.atomic_constituents = {C=1,H=1,N=1,O=1}
HCNO.charge = 0
HCNO.gamma = {
   value = 1.2154e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HCNO.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1382.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.647279890e00, 1.275053420e-02, -1.047942360e-05, 4.414328360e-09, -7.575214660e-13, 1.929902520e04, 1.073329720e01, }
   },
   { T_low  = 1382.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 6.598604560e00, 3.027786260e-03, -1.077043460e-06, 1.716665280e-10, -1.014393910e-14, 1.796613390e04, -1.033065990e01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
