-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of carbon monoxide, a 3Pir

CO_a3 = {}
CO_a3.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CO_a3.atomic_constituents = {C=1,O=1}
CO_a3.charge = 0
CO_a3.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

CO_a3.species_type = "polar diatomic"
CO_a3.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
CO_a3.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
CO_a3.h_f = {
   value = 16847048.60,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_a3.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_a3.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  0.0,       1.2057,  6,  41020.00,  1743.410,  14.3600, -4.500E-02,  0.000E+00,  1.69124,  1.904E-02,  6.360E-06,  4.000E-08,  4.153E+01,  1,  3 }
   -- ===========================================================================================================================================================
}
