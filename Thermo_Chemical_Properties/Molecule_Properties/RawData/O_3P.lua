-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of atomic Oxygen, 2s22p4 3P

O_3P = {}
O_3P.M = {
   value = 15.9994e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermO_3P.inp'
}
O_3P.atomic_constituents = {O=1}
O_3P.charge = 0
O_3P.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

O_3P.species_type = "monatomic"
O_3P.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
O_3P.h_f = {
   value = 15632317.75,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermO_3P.inp'
}
O_3P.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
O_3P.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
O_3P.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,       0.0,     9,   -1,    1,    1,    2 }
   -- ===========================================================
}

