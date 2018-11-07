-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of nitrogen, A3Sigmau+

N2_A = {}
N2_A.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_A.atomic_constituents = {N=2}
N2_A.charge = 0
N2_A.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

N2_A.species_type = "nonpolar diatomic"
N2_A.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
N2_A.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
N2_A.h_f = {
   value = 21438654.26,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_A.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_A.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  0.0,       1.2864,  3,  28980.00,  1460.941,  13.9800,  2.400E-02, -2.560E-03,  1.45390,  1.750E-02,  5.780E-06,  0.000E+00,  0.000E+00,  0,  3 }
   -- ===========================================================================================================================================================
}
   
