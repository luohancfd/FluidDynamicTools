-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Diatomic nitrogen cation excited state 'B' (B2 Sigma u+)

N2_plus_B = {}
N2_plus_B.M = { 
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_B.atomic_constituents = {N=2}
N2_plus_B.charge = 1
N2_plus_B.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- thermal nonequilibrium data

N2_plus_B.species_type = "nonpolar diatomic"
N2_plus_B.s_0 = {
   value = 7055.69,
   units = 'J/kg-k',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
N2_plus_B.h_f = {
   value = 5.3886282e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_B.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
N2_plus_B.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_plus_B.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat (confirmed with NIST, slightly different data for upper levels)',
   -- ===========================================================================================================================================================
   --   n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = { 25461.11,  1.0742,  2,  44730.00,  2421.140,  24.0700, -3.000E-01, -6.670E-02,  2.08507,  2.120E-02,  6.170E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
