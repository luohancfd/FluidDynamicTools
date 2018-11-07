-- Collater: Daniel F Potter
-- Date: 21-Jun-2009

CN_plus = {}
CN_plus.M = {
   value = 26.0168514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CN_plus.atomic_constituents = {C=1,N=1}
CN_plus.charge = 1
CN_plus.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CN_plus.CEA_coeffs = {
  { T_low  = 298.15,
    T_high = 1000.0,
    coeffs = { -8.302909570e+05,  8.775687500e+03, -2.977443560e+01,
    	        4.976897060e-02, -1.302225951e-05, -2.058325353e-08,
    	        1.126843895e-11,  1.703860539e+05,  2.039918818e+02
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -7.153463080e+06,  1.857250421e+04, -1.084534159e+01,
    	        6.106681430e-03, -1.191208566e-06,  1.184848778e-10,
    	       -4.799838730e-15,  9.242644960e+04,  1.135340573e+02
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = { -2.354919695e+08,  1.433776703e+05, -2.975360271e+01,
    	        4.280545600e-03, -2.707260413e-07,  8.178340660e-12,
    	       -9.629506200e-17, -9.229047140e+05,  2.964624987e+02
    }
  },

  ref="from CEA2::thermo.inp"
}
CN_plus.viscosity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.62526577e+00, B = -0.31779652e+02, C = -0.16407983e+04, D = 0.17454992e+01 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.87395209e+00, B = 0.56152222e+03, C = -0.17394809e+06, D = -0.39335958e+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.88503551e+00, B = 0.90902171e+03, C = -0.73129061e+06, D = -0.53503838e+00 }
  },
  ref = "from CEA2::trans.inp data for CO"
}
CN_plus.thermal_conductivity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.85439436E+00, B = 0.10573224E+03, C = -0.12347848E+05, D = 0.47793128E+00 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.88407146E+00, B = 0.13357293E+03, C = -0.11429640E+05, D = 0.24417019E+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.24175411E+01, B = 0.80462671E+04, C =  0.31090740E+07, D = -0.14516932E+02 }
  },
  ref = "from CEA2::trans.inp data for CO"
}

-- Nonequilibrium data

CN_plus.species_type = "polar diatomic"
CN_plus.s_0 = {
   value = 8203.91,
   units = 'J/kg-K',
   description = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/',
   reference = 'NA'
}
CN_plus.h_f = {
   value = 69143297.79,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN_plus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CN_plus.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN_plus.electronic_levels = {
   n_levels = 5,
   ref = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/',
   -- NOTE: dzero estimated from cea2 h_f values for CN+ -> C+ + N
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1729,  1,  40395.14,  2033.050,  16.1400,  0.000e+00,  0.000e+00,  1.89640,  1.880e-02,  7.000e-06,  0.000e+00,  0.000e+00,  0,  1 },
   ilev_1  = {   8313.60,  1.2473,  2,  32081.54,  1688.350,  15.1200,  0.000e+00,  0.000e+00,  1.67670,  1.910e-02,  6.840e-06,  0.000e+00,  0.000e+00,  1,  1 },
   ilev_2  = {  31771.00,  1.3640,  1,   8624.14,  1265.000,  11.0000,  0.000e+00,  0.000e+00,  1.40300,  0.200e-02,  1.300e-05,  0.000e+00,  0.000e+00,  0,  1 },
   ilev_3  = {  45533.60,  1.1710,  1,  40000.00,  2670.500,  46.9000,  0.000e+00,  0.000e+00,  1.90300,  3.200e-02,  4.700e-06,  0.000e+00,  0.000e+00,  0,  1 },
   ilev_4  = {  46253.00,  1.2762,  2,  20000.00,   890.760,   0.0000,  0.000e+00,  0.000e+00,  1.60180,  0.000e+00,  1.950e-05,  0.000e+00,  0.000e+00,  1,  1 }
   -- ===========================================================================================================================================================
}
