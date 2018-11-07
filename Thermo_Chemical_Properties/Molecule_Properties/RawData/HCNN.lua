-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HCNN = {}
HCNN.M = {
   value = 41.032040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
HCNN.atomic_constituents = {C=1,H=1,N=2}
HCNN.charge = 0
HCNN.gamma = {
   value = 1.2032e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
HCNN.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.524319400e00, 1.596061900e-02, -1.881635400e-05, 1.212554000e-08, -3.235737800e-12, 5.426198400e04, 1.167587000e01, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 5.894636200e00, 3.989595900e-03, -1.598238000e-06, 2.924939500e-10, -2.009468600e-14, 5.345294100e04, -5.103050200e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
