-- Collater: Rowan J. Gollan
-- Date: 20-Jun-2009

CH = {}
CH.M = {
   value = 13.0186400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH.atomic_constituents = {C=1,H=1}
CH.charge = 0
CH.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  2.220590133e+04, -3.405411530e+02,  5.531452290e+00,
	       -5.794964260e-03,  7.969554880e-06, -4.465911590e-09,
	        9.596338320e-13,  7.240783270e+04, -9.107673050e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.060763440e+06, -5.396206660e+03,  7.856293850e+00, 
               -7.965907450e-04,  1.764308305e-07, -1.976386267e-11,
	        5.030429510e-16,  1.062236592e+05, -3.154757439e+01
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = { -8.068368690e+08,  4.575450540e+05, -9.843975080e+01,
	        1.235244098e-02, -8.485608570e-07,  3.040410624e-11,
               -4.400315170e-16, -3.595851590e+06,  8.953477440e+02
    }
  },

  ref="from CEA2::thermo.inp"
}
CH.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
CH.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}

-- Nonequilibrium data

CH.species_type = "polar diatomic"
CH.s_0 = {
   value = 14059.84,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CH.h_f = {
   value = 45885791.76,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CH.I = {
   value = 78856621.05,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CH.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CH.electronic_levels = {
   n_levels = 5,
   ref = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1199,  4,  23422.80,  2858.500,  63.020,   0.000E+00,  0.000E+00, 14.45700,  5.340E-01, 1.4500E-03,  0.000E+00,  0.000E+00,  1,  2 },
   ilev_1  = {   5844.00,  1.0850,  4,  17578.80,  3145.000,  72.000,   0.000E+00,  0.000E+00, 15.40000,  5.500E-01, 0.0000E+00,  0.000E+00,  0.000E+00,  0,  4 },
   ilev_2  = {  23189.80,  1.1019,  4,  20000.00,  2930.700,  96.650,   0.000E+00,  0.000E+00, 14.93400,  6.970E-01, 1.5400E-03,  4.000E-05,  0.000E+00,  2,  2 },
   ilev_3  = {  26044.00,  1.1975,  2,  20000.00,  1794.900,   0.000,   0.000E+00,  0.000E+00, 12.64500,  0.000E+00, 2.2200E-03,  0.000E+00,  0.000E+00,  0,  2 },
   ilev_4  = {  31801.50,  1.1143,  2,  20000.00,  2840.200, 125.960,   1.355E+01,  0.000E+00, 14.60300,  7.185E-01, 1.5550E-03,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
