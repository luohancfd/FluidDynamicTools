-- Collater: Rowan J. Gollan
-- Date: 21-Jun-2009

NH = {}
NH.M = {
   value = 15.0146400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
NH.atomic_constituents = {N=1,H=1}
NH.charge = 0
NH.gamma = {
   value = 1.398,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
NH.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  1.359651320e+04, -1.900296604e+02,  4.518496790e+00,
               -2.432776899e-03,  2.377587464e-06, -2.592797084e-10,
               -2.659680792e-13,  4.280972190e+04, -3.886561616e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.958141991e+06, -5.782861300e+03,  9.335742020e+00,
               -2.292910311e-03,  6.076092480e-07, -6.647942750e-11,
                2.384234783e-15,  7.898912340e+04, -4.116970400e+01
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = {  9.524636790e+07, -8.585826910e+04,  2.980445181e+01,
               -2.979563697e-03,  1.656334158e-07, -4.744791840e-12,
                5.570148290e-17,  6.961434270e+05, -2.229027419e+02
    }
  },
  ref="from CEA2::thermo.inp"
}
NH.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.83724737e+00, B=0.43997150e+03, C=-0.17450753e+06, D=0.10365689e+00},
      {T_low=5000.0, T_high=15000.0, A=0.89986588e+00, B=0.14112801e+04, C=-0.18200478e+07, D=-0.55811716e+00},
      ref = 'from CEA2::trans.inp data for N'
   }
}
NH.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.83771661e+00, B=0.44243270e+03, C=-0.17578446e+06, D=0.89942915e+00},
      {T_low=5000.0, T_high=15000.0, A=0.90001710e+00, B=0.14141175e+04, C=-0.18262403e+07, D=0.24048513e+00},
      ref = 'from CEA2::trans.inp data for N'
   }
}

-- Nonequilibrium data

NH.species_type = "polar diatomic"
NH.s_0 = {
   value = 12071.55,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
NH.h_f = {
   value = 23778925.17,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
NH.I = {
   value = 84181637.35,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
NH.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
NH.electronic_levels = {
   n_levels = 6,
   ref = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.0362,  3,  26133.46,  3282.200,  78.300,   0.000E+00,  0.000E+00, 16.69930,  6.490E-01, 17.097E-04,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_1  = {  12566.00,  1.0340,  2,  13567.46,  3188.000,  68.000,   0.000E+00,  0.000E+00, 16.43900,  6.600E-01, 16.200E-04,  0.000E+00,  0.000E+00,  2,  1 },
   ilev_2  = {  21202.00,  1.0360,  1,   4931.46,  3352.400,  74.200,   7.000E-01,  0.000E+00, 16.70500,  5.910E-01, 16.000E-04,  0.000E+00,  0.000E+00,  0,  1 },
   ilev_3  = {  29807.00,  1.0369,  6,  26000.00,  3231.200,  98.600,   0.000E+00,  0.000E+00, 16.67450,  7.454E-01, 17.800E-04,  0.000E+00,  0.000E+00,  1,  3 },
   ilev_4  = {  43744.00,  1.1100,  2,  15000.00,  2122.600,  00.000,   0.000E+00,  0.000E+00, 14.53700,  5.930E-01, 22.000E-04,  0.000E+00,  0.000E+00,  1,  1 },
   ilev_5  = {  83160.00,  1.1163,  1,  20000.00,  2672.600,  71.200,   0.000E+00,  0.000E+00, 14.39000,  6.210E-01, 16.000E-04,  0.000E+00,  0.000E+00,  0,  1 }
   -- ===========================================================================================================================================================
}
