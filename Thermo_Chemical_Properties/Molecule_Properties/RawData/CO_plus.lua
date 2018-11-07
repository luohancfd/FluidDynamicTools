-- Collater: Daniel F. Potter
-- Date: 07 Aug 2009

CO_plus = {}
CO_plus.M = {
   value = 28.0095514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'cea2::thermo.inp'
}
CO_plus.atomic_constituents = {C=1,O=1}
CO_plus.charge = 1
CO_plus.gamma = {
   value = 1.4,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'from theory assuming negligible electronic and vibrational contributions'
}
CO_plus.CEA_coeffs = {
   { T_low  = 298.150,
     T_high = 1000.0,
     coeffs = { -2.178786658e+04,  1.288857032e+02,  3.769057550e+00, 
                -3.431730130e-03,  8.193945750e-06, -6.463814690e-09,
                 1.803727574e-12,  1.482345898e+05,  3.990547070e+00 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {  2.316847506e+05, -1.057646148e+03,  4.554257780e+00,
     	         4.495520320e-04, -2.489507047e-07,  5.267566420e-11,
                -3.289510270e-15,  1.555050724e+05, -3.873462640e+00 }
   },
   { T_low  = 6000.0,
     T_high = 20000.0,
     coeffs = { -3.035604054e+08,  2.393118392e+05, -7.034999240e+01,
     	         1.139551440e-02, -8.315173100e-07,  2.863705515e-11,
                -3.803269410e-16, -1.688617704e+06,  6.291980420e+02 }
   },
   ref='Gurvich (1991) from cea2::thermo.inp'
}

-- Nonequilibrium data

CO_plus.species_type = "polar diatomic"
CO_plus.s_0 = {
   value = 2998.98,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CO_plus.h_f = {
   value = 44548703.63,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_plus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CO_plus.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_plus.electronic_levels = {
   n_levels = 4,
   ref = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/',
   -- NOTE: dzero from Reddy and Viswanath J. Astrophys. Astr. (1990) 11, 67-72
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1151,  2,  67186.06,  2214.240,  15.1640, -0.0007,     0.000E+00,  1.97720,  1.896e-02,  6.350E-06,  0.000e+00,  0.000E+00,  0,  2 },
   ilev_1  = {   20733.3,  1.2438,  4,  46452.76,  1562.060,  13.5320,  0.0131,     0.000E+00,  1.58940,  1.942e-02,  6.600e-06,  0.000e+00,  0.000e+00,  1,  2 },
   ilev_2  = {   45876.7,  1.1688,  2,  21309.36,  1734.180,  27.9270,  0.3283,     0.000E+00,  1.79920,  3.025e-02,  7.750e-06,  2.200e-07,  0.000e+00,  0,  2 },
   ilev_3  = {   63012.0,  1.3460,  4,   4174.06,  1144.000,  33.3000,  0.0000,     0.000E+00,  1.35700,  2.400e-02,  0.000e+00,  0.000e+00,  0.000e+00,  2,  2 }
   -- ===========================================================================================================================================================
}
