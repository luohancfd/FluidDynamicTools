-- Author: Rowan J. Gollan
-- Date: 30-Oct-2008

-- diatomic iodine

I2 = {}
I2.M = { 
   value = 253.8089400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
I2.atomic_constituents = {I=2}
I2.charge = 0
I2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

I2.e_zero = {
   value = 245933.0, 
   units = 'J/kg',
   description = 'reference energy, from CEA enthalpy of formation'
}

I2.CEA_coeffs = {
   { T_low = 200.0,
     T_high = 1000.0,
     coeffs = { -5.087968770e+03, -1.249585210e+01,  4.504219090e+00,
		 1.370962533e-04, -1.390523014e-07,  1.174813853e-10,
		-2.337541043e-14,  6.213469810e+03,  5.583836940e+00 }
   },
   { T_low = 1000.000,
     T_high = 6000.000,
     coeffs = { -5.632594160e+06,  1.793961560e+04, -1.723055169e+01,
                 1.244214080e-02, -3.332768580e-06,  4.125477940e-10,
		-1.960461713e-14, -1.068505292e+05,  1.600531883e+02 }
   }
}

-- Presently these values are oxygen values.
I2.d = { 
   value = 3.433e-10,
   units = 'm',
   description = 'equivalent hard sphere diameter, based on L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
I2.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.919e-05, T_ref = 273.0, S = 139.0,
      ref = "Table 1-2, White (2006)"
   }
}
I2.thermal_conductivity = { 
   model = "Sutherland",
   parameters = { 
      k_ref = 0.0244, T_ref = 273.0, S = 240.0,
      ref = "Table 1-3, White (2006)"
   }
}
