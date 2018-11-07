-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C3H7 = {}
C3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C3H7.atomic_constituents = {C=3,H=7}
C3H7.charge = 0
C3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C3H7.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 1.051551800e00, 2.599198000e-02, 2.380054000e-06, -1.960956900e-08, 9.373247000e-12, 1.063186300e04, 2.112255900e01, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 7.702698700e00, 1.604420300e-02, -5.283322000e-06, 7.629859000e-10, -3.939228400e-14, 8.298433600e03, -1.548018000e01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
