-- Author: Rowan J. Gollan
-- Date: 30-Oct-2008

-- hydrogen iodine

HI = {}
HI.M = { 
   value = 127.9124100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
HI.atomic_constituents = {H=1,I=1}
HI.charge = 0
HI.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

HI.e_zero = {
   value = 206070.7,
   units = 'J/kg',
   description = 'reference energy, from CEA (enthalpy of formation)'
}

HI.CEA_coeffs = {
   { T_low = 200.000,
     T_high = 1000.000,
     coeffs = { 1.872881730e+04, -3.431788840e+02,  5.956712430e+00,
               -8.543439600e-03,  1.454780274e-05, -1.049104164e-08,
                2.839734003e-12,  3.682950720e+03, -8.149756090e+00 }
   },
   { T_low = 1000.000,
     T_high = 6000.000,
     coeffs = { 4.724921450e+05, -1.923465741e+03,  5.758048970e+00,
               -4.066266380e-04,  9.474332050e-08, -1.033534431e-11,
                4.611614790e-16,  1.394857037e+04, -1.182487652e+01 }
   }
}

-- Presently these values are oxygen values.
HI.d = { 
   value = 3.433e-10,
   units = 'm',
   description = 'equivalent hard sphere diameter, based on L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
HI.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.919e-05, T_ref = 273.0, S = 139.0,
      ref = "Table 1-2, White (2006)"
   }
}
HI.thermal_conductivity = { 
   model = "Sutherland",
   parameters = { 
      k_ref = 0.0244, T_ref = 273.0, S = 240.0,
      ref = "Table 1-3, White (2006)"
   }
}
