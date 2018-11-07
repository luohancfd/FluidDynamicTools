-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

NH3 = {}
NH3.M = {
   value = 17.030520e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
NH3.atomic_constituents = {N=1,H=3}
NH3.charge = 0
NH3.gamma = {
   value = 1.3036e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
NH3.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.286027400e00, -4.660523000e-03, 2.171851300e-05, -2.280888700e-08, 8.263804600e-12, -6.741728500e03, -6.253727700e-01, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.634452100e00, 5.666256000e-03, -1.727867600e-06, 2.386716100e-10, -1.257878600e-14, -6.544695800e03, 6.566292800e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
