-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CH2CO = {}
CH2CO.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CH2CO.atomic_constituents = {C=2,H=2,O=1}
CH2CO.charge = 0
CH2CO.gamma = {
   value = 1.1908e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
CH2CO.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.135836300e00, 1.811887210e-02, -1.739474740e-05, 9.343975680e-09, -2.014576150e-12, -7.042918040e03, 1.221564800e01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.511297320e00, 9.003597450e-03, -4.169396350e-06, 9.233458820e-10, -7.948382010e-14, -7.551053110e03, 6.322472050e-01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
