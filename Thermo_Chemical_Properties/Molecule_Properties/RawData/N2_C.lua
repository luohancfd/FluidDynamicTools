-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 10th excited state of nitrogen, C3Piu

N2_C = {}
N2_C.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_C.atomic_constituents = {N=2}
N2_C.charge = 0
N2_C.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

N2_C.species_type = "nonpolar diatomic"
N2_C.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
N2_C.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
N2_C.h_f = {
   value = 38064450.92,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_C.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_C.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0 = {  0.0,       1.1487,  6,   8831.00,  2047.780,  28.9488,  0.000E+00,  0.000E+00,  1.82473,  1.868E-02,  0.000E+00,  0.000E+00,  3.920E+01,  1,  3 }
   -- ===========================================================================================================================================================
}
   
