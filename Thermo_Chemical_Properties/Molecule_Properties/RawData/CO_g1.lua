-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- excited states 1 (a3)to 4 (e3) of carbon monoxide

CO_g1 = {}
CO_g1.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CO_g1.atomic_constituents = {C=1,O=1}
CO_g1.charge = 0
CO_g1.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

CO_g1.species_type = "polar diatomic"
CO_g1.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
CO_g1.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
CO_g1.h_f = {
   value = 16847048.60,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_g1.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_g1.electronic_levels = {
   n_levels = 4,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  0.0,       1.2057,  6,  41020.00,  1743.410,  14.3600, -4.500E-02,  0.000E+00,  1.69124,  1.904E-02,  6.360E-06,  4.000E-08,  4.153E+01,  1,  3 },
   ilev_1  = {   7138.79,  1.3523,  3,  34140.00,  1228.600,  10.4680,  9.100E-03,  0.000E+00,  1.34460,  1.892E-02,  6.410E-06,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_2  = {  12433.40,  1.3696,  6,  28870.00,  1171.940,  10.6350,  7.850E-02,  0.000E+00,  1.31080,  1.782E-02,  6.590E-06,  0.000E+00, -1.600E+01,  1,  3 },
   ilev_3  = {  15543.54,  1.3840,  3,  25790.00,  1117.720,  10.6860,  1.174E-01,  0.000E+00,  1.28360,  1.753E-02,  6.770E-06,  0.000E+00,  0.000E+00,  0,  3 },
   -- ===========================================================================================================================================================
}
