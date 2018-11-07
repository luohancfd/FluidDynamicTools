-- Author: Rowan J. Gollan
-- Date: 18-February-2009

air13 = {}
air13.M = { 
   value = 28.964e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864 & 867'
}
air13.atomic_constituents = {air13=1}
air13.charge = 0
air13.gamma = {
   value = 1.3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at moderate temperatures',
   reference = ''
}
air13.d = {
   value = 3.617e-10,
   units = 'm',
   description = 'equivalent hard-sphere diameter, sigma from L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
air13.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.716e-05, T_ref = 273.0, S = 111.0,
      ref = "Table 1-2, White (2006)"
   }
}
air13.thermal_conductivity = { 
   model = "Sutherland",
   parameters = {
      k_ref = 0.0241, T_ref = 273.0, S = 194.0,
      ref = "Table 1-3, White (2006)"
   }
}

