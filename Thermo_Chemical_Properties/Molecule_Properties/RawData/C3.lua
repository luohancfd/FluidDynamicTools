-- Author: Elise J. Fahy
-- 05 Feb 2013

C3 = {}
C3.M = {
   value = 36.0321e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C3.atomic_constituents = {C=3}
C3.charge = 0
C3.gamma = { 
   value = 1.2469,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K (=Cp/(Cp-R))',
   reference = 'evaluated using Cp @ 300.0K from Capitelli (2005) ESA STR-246'
}
C3.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {-4.354614480e+04,  6.660183220e+02,  1.451033157e+00,
                7.434513120e-03, -3.810152990e-06, -2.336961396e-11, 
                4.407054530e-13,  9.635170200e+04,  2.025173297e+01 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = { 4.508098930e+06, -1.461033761e+04,  2.281974644e+01,
	       -8.544340610e-03,  2.146069341e-06, -2.103867761e-10,
	        6.351589060e-15,  1.911976065e+05, -1.271869723e+02 }
   },
   { T_low  = 6000.0,
     T_high = 20000.0,
     coeffs = { 1.539589859e+08, -2.089057498e+05,  7.681111210e+01,
	       -8.939056190e-03,  5.594033240e-07, -1.743774353e-11,
	        2.181541208e-16,  1.650801763e+06, -6.081693320e+02 }
   },
   ref='Gurvich (1991) from cea2::thermo.inp'
}

-- C3.viscosity = {
--    model = "Blottner",
--    parameters = { 
--       A_mu = -0.01470, B_mu = 0.88110, C_mu = -13.5051,
--       ref = "Blottner (1971)"
--    }
-- }

C3.T_c = {
   value = 0.0,
   units = 'K',
   description = 'critical temperature',
   reference = 'none'
}

C3.p_c = { 
   value = 0.0,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'none'
}

C3.rho_c = { 
   value = 0.0,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'none',
}

-- Nonequilibrium data

C3.species_type = "linear nonpolar polyatomic"
C3.eps0 = {
   value = 7.3906e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = "Park (2001) Chemical-Kinetic Parameters of Hyperbolic Earth Entry",
}
C3.sigma = {
   value = 3.245e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = "Park (2001) Chemical-Kinetic Parameters of Hyperbolic Earth Entry",
}
C3.h_f = {
   value = 23311121.08,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C3.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'none'
}
C3.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
C3.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C3.theta_v = {
   value = 90.95,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'we[1] from electronic level table below'
}
C3.electronic_levels = {
   n_levels = 10,
   ref = 'Capitelli (2005) ESA STR-246',
   -- NOTE: Have put the we[1] vibrational temperature first so it is used for VT exchange as done by Park.
   --       Previously we[0] was the primary vibrational temperature. 
   --       Table has original notation in top row, and alternate notation underneath, since sources such as 	
   --       Capitelli use the alternate notation
-- =====================================================================================================================================
   --   n          Te    re       g    dzero       A0     B0        C0       sigma   sigma_rot      we[1]   we[0]     we[2]   we[3]
   --                           (p_i) (E_diff)                                                      nu[1]   nu[0]     nu[2]   nu[3]
--   =====================================================================================================================================
   ilev_0  = {    0.0,    1.277,    1,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_1  = {14000.0,    1.277,    6,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_2  = {21500.0,    1.277,    6,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_3  = {23800.0,    1.277,    3,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_4  = {24675.5,    1.277,    2,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_5  = {29100.0,    1.277,    6,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_6  = {32800.0,    1.277,    3,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_7  = {32900.0,    1.277,    1,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_8  = {33700.0,    1.277,    2,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   ilev_9  = {40500.0,    1.277,    2,    61570.75,    0.0,    0.4305,    0.0,          2,           0,      63.1,    1225.0,    63.1,    2040.0},
   -- ==================================================================================================================================
}

-- Real gas data  - not required at present for C3 because we are using thermally perfect cases only.

