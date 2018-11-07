-- Collater: Daniel F. Potter
-- Date: 23-Feb-2010

H_plus = {}
H_plus.M = {
   value = 1.0073914e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
H_plus.atomic_constituents = {H=1}
H_plus.charge = 1
H_plus.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
H_plus.CEA_coeffs = {
  { T_low  = 298.150,
    T_high = 1000.0,
    coeffs = { 0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
    	       0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
               0.000000000e+00,  1.840214877e+05, -1.140646644e+00 },
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { 0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
    	       0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
               0.000000000e+00,  1.840214877e+05, -1.140646644e+00 },
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = { 0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
    	       0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
               0.000000000e+00,  1.840214877e+05, -1.140646644e+00 },
  },
  ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
H_plus.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.74226149e+00, B=-0.40132865e+03, C=0.18554165e+06, D=0.46741844e-01},
      {T_low=5000.0, T_high=15000.0, A=0.87486623e+00, B=-0.25022902e+04, C=0.70955048e+07, D=-0.93888455e+00},
      ref = 'from CEA2::trans.inp data for H)'
   }
}
H_plus.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.74166119e+00, B=-0.40487203e+03, C=0.18775642e+06, D=0.34843121e+01},
      {T_low=5000.0, T_high=15000.0, A=0.87447639e+00, B=-0.25089452e+04, C=0.71081294e+07, D=0.24970991e+01},
      ref = 'from CEA2::trans.inp data for H'
   }
}

-- Thermal nonequilibrium data

H_plus.species_type = "monatomic"
H_plus.s_0 = {
   value = 108150.62,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
H_plus.h_f = {
   value = 1524974233.45,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
H_plus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
H_plus.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
H_plus.electronic_levels = {
   -- NOTE: H+ has no electrons - we model this is a zero state with a degeneracy of one
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.     n      E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0  = {  1,        0.00,    1,   -1,   -1,   -1,   -1 },
   -- ===========================================================
}

