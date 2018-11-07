-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Free electron

e_minus = {}
e_minus.M = {
   value = 0.000548579903e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
e_minus.atomic_constituents = {}
e_minus.charge = -1
e_minus.gamma = {
   value = 1.667, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
e_minus.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=2000.0, T_high=5000.0, A=0.59319174e+01, B=0.56594215e+04, C=-0.22576125e+07, D=-0.53458874e+02},
      ref = 'from CEA2 data for O'
   }
}
e_minus.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=2000.0, T_high=5000.0, A=0.59320964e+01, B=0.56601476e+04, C=-0.22577332e+07, D=-0.42512600e+02},
      ref = 'from CEA2 data for O'
   }
}
e_minus.CEA_coeffs = {
  { T_low = 298.150,
    T_high = 1000.0,
    coeffs = {  0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
                0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
                0.000000000e+00, -7.453750000e+02, -1.172081224e+01
    }
  },
  { T_low = 1000.0,
    T_high = 6000.0,
    coeffs = {  0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
                0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
                0.000000000e+00, -7.453750000e+02, -1.172081224e+01
    }
  },
  { T_low = 6000.0,
    T_high = 20000.0,
    coeffs = {  0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
                0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
                0.000000000e+00, -7.453750000e+02, -1.172081224e+01
    }
  }
}

-- thermal nonequilibrium data

e_minus.species_type = "free electron"
e_minus.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NA'
}
e_minus.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::therm.inp'
}
e_minus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
e_minus.Z = {
   value = -1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
e_minus.electronic_levels = {
   n_levels = 1,
   ref = 'Park (1991) <Nonequilibrium hypersonic aerothermodynamics>',
   -- ===================================
   --  Level   n      E(cm-1)        g
   -- ===================================
   ilev_0  = { 1,       0.00,        2 },
   -- ===================================
}
  
