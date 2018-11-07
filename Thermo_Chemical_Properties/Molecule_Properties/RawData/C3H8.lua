-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C3H8 = {}
C3H8.M = {
   value = 44.095620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C3H8.atomic_constituents = {C=3,H=8}
C3H8.charge = 0
C3H8.gamma = {
   value = 1.1267e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C3H8.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 9.335538100e-01, 2.642457900e-02, 6.105972700e-06, -2.197749900e-08, 9.514925300e-12, -1.395852000e04, 1.920169100e01, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 7.534136800e00, 1.887223900e-02, -6.271849100e-06, 9.147564900e-10, -4.783806900e-14, -1.646751600e04, -1.789234900e01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
C3H8.T_c = {
   value = 369.83,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C3H8.p_c = {
   value = 42.48e05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
