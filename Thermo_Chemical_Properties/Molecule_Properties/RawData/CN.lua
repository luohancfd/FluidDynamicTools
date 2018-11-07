-- Collater: Rowan J. Gollan
-- Date: 21-Jun-2009

CN = {}
CN.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CN.atomic_constituents = {C=1,N=1}
CN.charge = 0
CN.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CN.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  3.949148570e+03, -1.391590572e+02,  4.930835320e+00,
               -6.304670510e-03,  1.256836472e-05, -9.878300500e-09,
                2.843137221e-12,  5.228455380e+04, -2.763115585e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -2.228006270e+06,  5.040733390e+03, -2.121897722e-01,
                1.354901134e-03,  1.325929798e-07, -6.937006370e-11,
                5.494952270e-15,  1.784496132e+04,  3.282563919e+01
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = { -1.794798118e+08,  1.054346069e+05, -1.729624170e+01,
                2.194895530e-03, -8.508938030e-08,  9.318692990e-13,
                6.358139930e-18, -7.962594120e+05,  1.913139639e+02
    }
  },

  ref="from CEA2::thermo.inp"
}
-- CN.viscosity = {
--    model = "Blottner",
--    parameters = { 
--       A_mu = -0.0250, B_mu = 0.6810, C_mu = -12.4327,
--       ref = "Blottner (1971)"
--    }
-- }
CN.viscosity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.62526577e+00, B = -0.31779652e+02, C = -0.16407983e+04, D = 0.17454992e+01 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.87395209e+00, B = 0.56152222e+03, C = -0.17394809e+06, D = -0.39335958e+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.88503551e+00, B = 0.90902171e+03, C = -0.73129061e+06, D = -0.53503838e+00 }
  },
  ref = "from CEA2::trans.inp data for CO"
}
CN.thermal_conductivity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.85439436E+00, B = 0.10573224E+03, C = -0.12347848E+05, D = 0.47793128E+00 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.88407146E+00, B = 0.13357293E+03, C = -0.11429640E+05, D = 0.24417019E+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.24175411E+01, B = 0.80462671E+04, C =  0.31090740E+07, D = -0.14516932E+02 }
  },
  ref = "from CEA2::trans.inp data for CO"
}

-- Nonequilibrium data

CN.species_type = "polar diatomic"
CN.eps0 = {
   value = 1.0354935e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
CN.sigma = {
   value = 3.856e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
CN.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NA'
}
CN.h_f = {
   value = 16861160.30,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN.I = {
   value = 52280679.54,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'Derived from CEA formation energies'
}
CN.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN.electronic_levels = {
   -- n_levels = 7,
   n_levels = 5,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1718,  2,  61631.00,  2068.640,  13.0988, -1.130E-02,  3.200E-04,  1.89978,  1.736E-02,  6.400E-06,  1.200E-08,  0.000E+00,  0,  2 },
   ilev_1  = {   9245.28,  1.2333,  4,  52580.00,  1812.560,  12.6090, -1.180E-02,  0.000E+00,  1.71510,  1.708E-02,  5.930E-06,  4.200E-08, -5.264E+01,  1,  2 },
   ilev_2  = {  25752.80,  1.1511,  2,  55130.00,  2161.420,  18.1200, -4.309E-01,  0.000E+00,  1.96882,  2.005E-02,  6.600E-06,  0.000E+00,  0.000E+00,  0,  2 },
   ilev_3  = {  32400.00,  1.2300,  4,  40030.00,  1000.000,   0.0000,  0.000E+00,  0.000E+00,  1.00000,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  4 },
   ilev_4  = {  54486.30,  1.4980,  4,  26970.00,  1004.710,   8.7800,  0.000E+00,  0.000E+00,  1.16200,  1.300E-02,  7.000E-06,  0.000E+00, -3.030E+00,  1,  2 },
   ilev_5  = {  59151.18,  1.3245,  2,  21970.00,  1681.430,   3.6000, -1.020E+00,  0.000E+00,  1.48710,  6.430E-03,  5.000E-06,  0.000E+00,  0.000E+00,  0,  2 },
   ilev_6  = {  60095.64,  1.3732,  4,  21250.00,  1239.500,  12.7500,  0.000E+00,  0.000E+00,  1.38340,  1.870E-02,  7.000E-06,  0.000E+00,  2.786E+01,  2,  2 },
   -- ===========================================================================================================================================================
}
