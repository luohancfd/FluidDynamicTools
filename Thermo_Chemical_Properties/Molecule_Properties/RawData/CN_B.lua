-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of cyanogen, B2Sigma+

CN_B = {}
CN_B.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CN_B.atomic_constituents = {C=1,N=1}
CN_B.charge = 0
CN_B.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

CN_B.species_type = "polar diatomic"
CN_B.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
CN_B.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
CN_B.h_f = {
   value = 28702164.21,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN_B.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN_B.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  0.0,       1.1511,  2,  55130.00,  2161.420,  18.1200, -4.309E-01,  0.000E+00,  1.96882,  2.005E-02,  6.600E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
