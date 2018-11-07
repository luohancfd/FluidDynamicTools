-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H3 = {}
C2H3.M = {
   value = 27.045220e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H3.atomic_constituents = {C=2,H=3}
C2H3.charge = 0
C2H3.gamma = {
   value = 1.2408e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H3.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.212466450e00, 1.514791620e-03, 2.592094120e-05, -3.576578470e-08, 1.471508730e-11, 3.485984680e04, 8.510540250e00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.016724000e00, 1.033022920e-02, -4.680823490e-06, 1.017632880e-09, -8.626070410e-14, 3.461287390e04, 7.787323780e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
