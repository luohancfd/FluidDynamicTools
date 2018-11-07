-- Author: Chris James
-- Date: 08-October-2012
--  08-October-2012: starting building this from the He file.

Ne = {}
Ne.M = {
   value = 20.1797000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'cea2::thermo.inp'
}
Ne.atomic_constituents = {Ne=1}
Ne.charge = 0
Ne.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
Ne.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 6000.0,
     coeffs = { 0.000000000e+00, 0.000000000e+00, 2.500000000e+00,
		0.000000000e+00, 0.000000000e+00, 0.000000000e+00,
		0.000000000e+00, -7.453750000e+02, 3.355322720e+00,
	      }
   },
   { T_low  = 6000.0,
     T_high = 20000.0,
     coeffs = { -1.238252746e+07, 6.958579580e+03, 1.016709287e+00,
		1.424664555e-04, -4.803933930e-09,-1.170213183e-13,
		8.415153652e-18, -5.663933630e+04, 1.648438697e+01
	      }
   },
   ref='cea2::thermo.inp'
   -- NOTE: Data for first 2 ranges has been combined as the coefficients are the same - CJ
}
Ne.d = {
   value = 2.789e-10,
   units = 'm',
   description = 'equivalent hard-sphere diameter, sigma from L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2007), p. 864'
}
Ne.viscosity = {
   model = "CEA",
   parameters = {
	{T_low=200.0, T_high=1000.0, A=0.68398511E+00, B=0.18732366E+02, C=-0.23663189E+04, D=0.18284755E+01},
	{T_low=1000.0, T_high=5000.0, A=0.72333495E+00, B=0.10420872E+03, C=-0.25429545E+05, D=0.14942434E+01},
	{T_low=5000.0, T_high=15000.0, A=0.77549350E+00, B=0.59414850E+03, C=-0.69670786E+06, D=0.97885712E+00},
      ref = "from CEA2::trans.inp which cites Bich et al. (1990)"
   }
}
Ne.thermal_conductivity = {
   model = "CEA",
   parameters = {
	{T_low=200.0, T_high=1000.0, A=0.68509965E+00, B=0.19794924E+02, C=-0.24525539E+04, D=0.22586136E+01},
	{T_low=1000.0, T_high=5000.0, A=0.72278122E+00, B=0.10528290E+03, C=-0.26355706E+05, D=0.19367337E+01},  
	{T_low=5000.0,  T_high=15000.0, A=0.77589413E+00, B=0.61283778E+03, C=-0.74015705E+06, D=0.14114011E+01},

      ref = "from CEA2::trans.inp which cites Bich et al. (1990)"
   }
}
Ne.T_c = {
   value = 44.40,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
Ne.p_c = {
   value = 2.76e+06,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}

-- Thermal nonequilibrium data

Ne.species_type = "monatomic"
Ne.eps0 = {
   value = 4.52852544e-22,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
Ne.sigma = {
   value = 2.820e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
Ne.s_0 = {
   value = 7251.34,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
Ne.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
Ne.I = {
   value = 103106674.71,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
Ne.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
