-- Data collected by Rowan J. Gollan
-- Date: 29-Jan-2013

Cl2 = {}
Cl2.M = {
   value = 70.9060e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
Cl2.atomic_constituents = {Cl=2}
Cl2.charge = 0
Cl2.species_type = "nonpolar diatomic"
Cl2.gamma = {
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
Cl2.d = {
   value = 4.217e-10,
   units = 'm',
   description = 'hard sphere diameter',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
Cl2.sigma = {
   value = 4.217e-10,
   note = 'see above for d'
}
Cl2.r0 = {
   value = 4.217e-10,
   note = 'see above for d'
}
Cl2.eps0 = {
   value = 4.36285558e-21,
   units = 'J',
   description = 'Depth of intermolecular potential minimum',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
Cl2.r_eq = {
   value = 1.987e-10,
   units = 'm',
   description = 'equilibrium intermolecular distance',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
Cl2.f_m = {
   value = 1.0,
   units = 'non-dimensional',
   description = 'mass factor'
}
Cl2.mu = {
   value = 5.88711e-26,
   units = 'kg/particle',
   description = 'reduced mass of consituent atoms'
}
Cl2.alpha = {
   value = 1.09,
   units = 'Angstrom^3',
   description = 'polarizability',
   reference = 'value for N2'
}
Cl2.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'heat of formation',
   reference = 'from CEA2::thermo.inp'
}
Cl2.s_0 = {
   value = 3146.15,
   units = 'J/(kg.K)',
   description = 'standard state entropy',
   reference = 'NIST Chemistry Webbook: 223.081/70.9060e-3'
}
Cl2.CEA_coeffs = {
   { T_low = 200.0,
     T_high = 1000.0,
     coeffs = { 3.462815170e+04, -5.547126520e+02,  6.207589370e+00,
		-2.989632078e-03, 3.173027290e-06, -1.793629562e-09,
		4.260043590e-13,  1.534069331e+03, -9.438331107e+00 }
   },
   { T_low = 1000.0,
     T_high = 6000.0,
     coeffs = { 6.092569420e+06, -1.949627662e+04,  2.854535795e+01,
		-1.449968764e-02,  4.463890770e-06, -6.358525860e-10,
		3.327360290e-14,  1.212117724e+05, -1.690778824e+02 }
   }
}
Cl2.viscosity = {
   model = "CEA",
   parameters = {
      {T_low = 300.0, T_high = 1000.0, A = 0.53516134e+00, B = -0.23624735e+03, C= 0.13738454e+05, D = 0.24970463e+01 },
      {T_low = 1000.0, T_high = 5000.0, A = 0.63348430e+00, B = -0.38786240e+02, C = -0.35830615e+05, D = 0.16699633e+01}
   },
   reference = 'from CEA2::trans.inp which cites Svehla (1994)'
}
Cl2.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low = 300.0, T_high = 1000.0, A = 0.34156262e+00, B = -0.46059166e+03, C = 0.34712872e+05, D = 0.37412367e+01},
      {T_low = 1000.0, T_high = 5000.0, A = 0.87392526e+00, B = 0.19876120e+03, C = -0.28784264e+05, D = -0.53204988e+00 }
   },
   reference = 'from CEA2::trans.inp which cites Svehla (1994)'
}
Cl2.theta_v = {
   value = 805.3,
   units = 'K',
   description = 'characteristic vibrational temperature'
}

