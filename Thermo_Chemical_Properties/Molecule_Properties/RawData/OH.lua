-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

OH = {}
OH.M = {
   value = 17.00734e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
OH.atomic_constituents = {O=1,H=1}
OH.charge = 0
OH.gamma = {
   value = 1.386,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
OH.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.59711536e+00, B=-0.46100678e+03, C=0.37606286e+05, D=0.24041761e+01},
      {T_low=5000.0, T_high=15000.0, A=0.64287721e+00, B=-0.18173747e+03, C=-0.88543767e+05, D=0.19636057e+01},
      ref = 'from CEA2::trans.inp which cites Svehla (1994)'
   }
}
OH.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.68627561e+00, B=-0.74033274e+03, C=0.27559033e+05, D=0.28308741e+01},
      {T_low=5000.0, T_high=15000.0, A=-0.47918112e+00, B=-0.93769908e+04, C=0.70509952e+07, D=0.14203688e+02},
      ref = 'from CEA2::trans.inp which cites Svehla (1994)'
   }
}
OH.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -1.998858990e+03,  9.300136160e+01,  3.050854229e+00,
	        1.529529288e-03, -3.157890998e-06,  3.315446180e-09,
	       -1.138762683e-12,  2.991214235e+03,  4.674110790e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.017393379e+06, -2.509957276e+03,  5.116547860e+00,
		1.305299930e-04, -8.284322260e-08,  2.006475941e-11,
               -1.556993656e-15,  2.019640206e+04, -1.101282337e+01
    }
  },
  { T_low = 6000.0,
    T_high = 20000.0,
    coeffs = {  2.847234193e+08, -1.859532612e+05,  5.008240900e+01,
               -5.142374980e-03,  2.875536589e-07, -8.228817960e-12,
                9.567229020e-17,  1.468393908e+06, -4.023555580e+02
    }
 },
  ref="from CEA2::thermo.inp"
}

-- Nonequilibrium data

OH.species_type = "polar diatomic"
OH.eps0 = {
   value = 1.0,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Find me!'
}
OH.sigma = {
   value = 1.0,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Find me!'
}
OH.r0 = {
   value = 1.0,
   units = 'm',
   description = 'Zero of the intermolecular potential',
   reference = 'Find me!'
}
OH.r_eq = {
   value = 1.0,
   units = 'm',
   description = 'Equilibrium intermolecular distance',
   reference = 'Find me!'
}
OH.f_m = {
   value = 1.0,
   units = 'ND',
   description = 'Mass factor = ( M ( Ma^2 + Mb^2 ) / ( 2 Ma Mb ( Ma + Mb ) )',
   reference = 'Calculate me!'
}
OH.mu = {
   value = 1.0,
   units = 'kg/particle',
   description = 'Reduced mass of constituent atoms',
   reference = 'Calculate me!'
}
OH.alpha = {
   value = 1.0,
   units = 'Angstrom^3',
   description = 'Polarizability',
   reference = 'Find me!'
}
OH.mu_B = {
   value = 1.0,
   units = 'Debye',
   description = 'Dipole moment',
   reference = 'Find me!'
}
OH.s_0 = {
   value = 10801.81,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
OH.h_f = {
   value = 2191889.27,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
OH.I = {
   value = 73847505.07,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
OH.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
OH.electronic_levels = {
   n_levels = 2,
   ref = 'Huber & Herzberg (1979), Luque & Crosley (1998)',
   -- =====================================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye            weze          be         alphae       de         betae        spn-orb     l   s  
   -- =====================================================================================================================================================================
   ilev_0  = {      0.00,  0.96966, 1,  37275.0,  3737.7941,  84.91456,  0.558406E+00, -2.59739E-02,  18.89638,  0.7242E+00,  19.38E-04, 4.320E-05,  -139.21E+00,  1,  2 },
   ilev_1  = {  32683.97,  1.0121,  1,  20406.2,  3178.3554,  92.68141, -1.773050E+00,  3.07923E-01,  17.38922,  0.7868E+00,  20.39E-04, 0.000E-00,    0.000E+00,  0,  2 },
   -- =====================================================================================================================================================================
}
