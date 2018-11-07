-- Author: Rowan J. Gollan
-- Date: 28-July-2008

air = {}
air.M = { 
   value = 28.964e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864 & 867'
}
air.atomic_constituents = {air=1}
air.charge = 0
air.gamma = {
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Anderson Jr. (1989), p. 441'
}
air.d = {
   value = 3.617e-10,
   units = 'm',
   description = 'equivalent hard-sphere diameter, sigma from L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
air.p_a = {
   value = 101325.0,
   units = 'Pa',
   description = 'reference pressure for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
air.rho_a = {
   value = 10.0,
   units = 'kg/m^3',
   description = 'reference density for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
air.k_s = {
   value = 1.0,
   units = '-',
   description = 'reference k_s for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
air.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.716e-05, T_ref = 273.0, S = 111.0,
      ref = "Table 1-2, White (2006)"
   }
}
air.thermal_conductivity = { 
   model = "Sutherland",
   parameters = {
      k_ref = 0.0241, T_ref = 273.0, S = 194.0,
      ref = "Table 1-3, White (2006)"
   }
}
air.CEA_coeffs = {
  { T_low = 200.0,
    T_high = 1000.0,
    coeffs = {  1.009950160e+04, -1.968275610e+02,  5.009155110e+00,
               -5.761013730e-03,  1.066859930e-05, -7.940297970e-09,
                2.185231910e-12, -1.767967310e+02, -3.921504225e+00 }
  },
  { T_low = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.415214430e+05, -1.257874600e+03,  5.144558670e+00,
               -2.138541790e-04,  7.065227840e-08, -1.071483490e-11,
                6.577800150e-16,  6.462263190e+03, -8.147411905e+00 }
  }
}

