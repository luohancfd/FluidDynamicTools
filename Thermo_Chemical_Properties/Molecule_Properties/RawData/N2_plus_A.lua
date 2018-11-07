-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Diatomic nitrogen cation excited state 'A' (A2 Pi u)

N2_plus_A = {}
N2_plus_A.M = { 
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_A.atomic_constituents = {N=2}
N2_plus_A.charge = 1
N2_plus_A.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- thermal nonequilibrium data

N2_plus_A.species_type = "nonpolar diatomic"
N2_plus_A.h_f = {
   value = 5.3886282e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_A.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_plus_A.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat (slightly different to NIST data eg Te = 9166.9)',
   -- ===========================================================================================================================================================
   --   n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  9167.46,  1.1749,  4,  61280.00,  1903.700,  15.1110,  1.120E-02, -2.700E-04,  1.74450,  1.883E-02,  5.600E-06,  1.800E-07, -7.462E+01,  1,  2 }
   -- ===========================================================================================================================================================
}
