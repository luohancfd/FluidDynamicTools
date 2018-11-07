-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 5th excited state of dicarbon, d3Pig

C2_d3 = {}
C2_d3.M = { 
   value = 24.0214000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
C2_d3.atomic_constituents = {C=2}
C2_d3.charge = 0
C2_d3.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

C2_d3.species_type = "nonpolar diatomic"
C2_d3.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
C2_d3.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
C2_d3.h_f = {
   value = 44542774.12,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2_d3.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2_d3.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  0.0,       1.2661,  6,  30120.00,  1788.220,  16.4400, -5.067E-01,  0.000E+00,  1.75270,  1.608E-02,  6.740E-06,  1.030E-07, -1.690E+01,  1,  3 }
   -- ===========================================================================================================================================================
}
