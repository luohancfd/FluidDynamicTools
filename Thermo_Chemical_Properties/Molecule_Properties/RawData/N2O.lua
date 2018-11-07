-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

N2O = {}
N2O.M = {
   value = 44.012800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
N2O.atomic_constituents = {N=2,O=1}
N2O.charge = 0
N2O.gamma = {
   value = 1.2735e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
N2O.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.257150200e00, 1.130472800e-02, -1.367131900e-05, 9.681980600e-09, -2.930718200e-12, 8.741774400e03, 1.075799200e01, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.823072900e00, 2.627025100e-03, -9.585087400e-07, 1.600071200e-10, -9.775230300e-15, 8.073404800e03, -2.201720700e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
N2O.T_c = {
   value = 309.60,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
N2O.p_c = {
   value = 72.55e05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
