-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of nitrogen, X1Sigmag+

N2_X = {}
N2_X.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_X.atomic_constituents = {N=2}
N2_X.charge = 0
N2_X.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

N2_X.species_type = "nonpolar diatomic"
N2_X.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
N2_X.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
N2_X.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_X.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_X.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.0977,  1,  78740.00,  2358.569,  14.3244, -2.258E-03, -2.400E-04,  1.99824,  1.732E-02,  5.760E-06,  0.000E+00,  0.000E+00,  0,  1 }
   -- ===========================================================================================================================================================
}
   
