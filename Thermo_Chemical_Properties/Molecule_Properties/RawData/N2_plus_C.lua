-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Diatomic nitrogen cation excited state 'C' (C2 Sigma u+)

N2_plus_C = {}
N2_plus_C.M = { 
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_C.atomic_constituents = {N=2}
N2_plus_C.charge = 1
N2_plus_C.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- thermal nonequilibrium data

N2_plus_C.species_type = "nonpolar diatomic"
N2_plus_C.s_0 = {
   value = 7055.69,
   units = 'J/kg-k',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
N2_plus_C.h_f = {
   value = 5.3886282e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_C.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
N2_plus_C.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_plus_C.electronic_levels = {
   n_levels = 5,
   ref = 'Spradian07::diatom.dat (confirmed with NIST)',
   -- ===========================================================================================================================================================
   --   n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = { 64609.03,  1.2628,  2,  24990.00,  2069.400,   8.3000, -6.300E-01,  1.300E-02,  1.50980, -1.000E-03,  4.000E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
