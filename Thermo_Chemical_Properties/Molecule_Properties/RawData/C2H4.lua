-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H4 = {}
C2H4.M = {
   value = 28.053160e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H4.atomic_constituents = {C=2,H=4}
C2H4.charge = 0
C2H4.gamma = {
   value = 1.2393e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H4.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.959201480e+00, -7.570522470e-03, 5.709902920e-05, -6.915887530e-08, 2.698843730e-11, 5.089775930e+03, 4.097330960e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.036111160e+00, 1.464541510e-02, -6.710779150e-06, 1.472229230e-09, -1.257060610e-13, 4.939886140e+03, 1.030536930e+01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
C2H4.T_c = {
   value = 282.34,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C2H4.p_c = {
   value = 50.41e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
