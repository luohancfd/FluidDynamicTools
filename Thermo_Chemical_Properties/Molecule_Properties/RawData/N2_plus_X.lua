-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Diatomic nitrogen cation ground state (X2 Sigma g+)

N2_plus_X = {}
N2_plus_X.M = { 
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_X.atomic_constituents = {N=2}
N2_plus_X.charge = 1
N2_plus_X.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- thermal nonequilibrium data

N2_plus_X.species_type = "nonpolar diatomic"
N2_plus_X.s_0 = {
   value = 7055.69,
   units = 'J/kg-k',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
N2_plus_X.h_f = {
   value = 5.3886282e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_plus_X.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
N2_plus_X.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_plus_X.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat (confirmed with NIST)',
   -- ===========================================================================================================================================================
   --   n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {     0.00,  1.1164,  2,  70300.00,  2207.220,  16.2260,  4.000E-03, -6.100E-03,  1.93171,  1.882E-02,  6.100E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
