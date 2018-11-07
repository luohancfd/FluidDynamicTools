-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H6 = {}
C2H6.M = {
   value = 30.069040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H6.atomic_constituents = {C=2,H=6}
C2H6.charge = 0
C2H6.gamma = {
   value = 1.1872e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H6.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.291424920e00, -5.501542700e-03, 5.994382880e-05, -7.084662850e-08, 2.686857710e-11, -1.152220550e+04, 2.666823160e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 1.071881500e+00, 2.168526770e-02, -1.002560670e-05, 2.214120010e-09, -1.900028900e-13, -1.142639320e+04, 1.511561070e+01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
C2H6.T_c = {
   value = 305.32,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C2H6.p_c = {
   value = 48.72e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
