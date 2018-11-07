-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of atomic Nitrogen, 2s22p3 2D¡

N_2Do = {}
N_2Do.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
N_2Do.atomic_constituents = {N=1}
N_2Do.charge = 0
N_2Do.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
   },
  ref="none"
}

-- Nonequilibrium data

N_2Do.species_type = "monatomic"
N_2Do.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'none'
}
N_2Do.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ionization energy',
   reference = 'none'
}
N_2Do.h_f = {
   value = 50168674.94,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N_2Do.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N_2Do.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,    0.0,         10,  -1,    2,    1,    1 }
   -- ===========================================================
}
