-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of dicarbon, X1Sigmag+

C2_X = {}
C2_X.M = { 
   value = 24.0214000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
C2_X.atomic_constituents = {C=2}
C2_X.charge = 0
C2_X.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

C2_X.species_type = "nonpolar diatomic"
C2_X.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
C2_X.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
C2_X.h_f = {
   value = 34571562.11,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2_X.I = {
   value = 45829859.25,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C2_X.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2_X.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.2425,  1,  50104.00,  1854.710,  13.3400, -1.720E-01,  0.000E+00,  1.81984,  1.765E-02,  6.920E-06,  8.100E-08,  0.000E+00,  0,  1 }
   -- ===========================================================================================================================================================
}
