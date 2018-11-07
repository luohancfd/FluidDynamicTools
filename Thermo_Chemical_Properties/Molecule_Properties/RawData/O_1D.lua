-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of atomic Nitrogen, 2s22p4 1D

O_1D = {}
O_1D.M = {
   value = 15.9994e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermO_1D.inp'
}
O_1D.atomic_constituents = {O=1}
O_1D.charge = 0
O_1D.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

O_1D.species_type = "monatomic"
O_1D.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
O_1D.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
O_1D.h_f = {
   value = 27438333.86,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermO_1D.inp'
}
O_1D.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
O_1D.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,     0.0,         5,   -1,    2,    0,    2 }
   -- ===========================================================
}

