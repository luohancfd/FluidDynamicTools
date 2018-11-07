-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CHO = {}
CHO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CHO.atomic_constituents = {C=1,H=1,O=1}
CHO.charge = 0
CHO.gamma = {
   value = 1.3161e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CHO.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.221185840e+00, -3.243925320e-03, 1.377994460e-05, -1.331440930e-08, 4.337688650e-12, 3.839564960e+03, 3.394372430e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.772174380e+00, 4.956955260e-03, -2.484456130e-06, 5.891617780e-10, -5.335087110e-14, 4.011918150e+03, 9.798344920e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
