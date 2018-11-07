-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of atomic Nitrogen, 2s22p3 2P¡

N_2Po = {}
N_2Po.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
N_2Po.atomic_constituents = {N=1}
N_2Po.charge = 0
N_2Po.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

N_2Po.species_type = "monatomic"
N_2Po.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
N_2Po.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
N_2Po.h_f = {
   value = 58377313.00,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N_2Po.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N_2Po.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,    0.0,         6,   -1,    1,    1,    1 }
   -- ===========================================================
}
