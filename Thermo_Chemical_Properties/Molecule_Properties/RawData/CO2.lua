-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

CO2 = {}
CO2.M = {
   value = 44.009500e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
CO2.atomic_constituents = {C=1,O=2}
CO2.charge = 0
CO2.gamma = {
   value = 1.2877e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
--[[
CO2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 0.000000000e+00,  0.000000000e+00, 2.356773520e+00,
                8.984596770e-03, -7.123562690e-06, 2.459190220e-09,
               -1.436995480e-13, -4.837196970e+04, 9.901052220e+00 }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = { 0.000000000e+00,  0.000000000e+00, 3.857460290e+00,
                4.414370260e-03, -2.214814040e-06, 5.234901880e-10,
               -4.720841640e-14, -4.875916600e+04, 2.271638060e+00 }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
--]]
CO2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 4.943650540e+04, -6.264116010e+02,  5.301725240e+00,
                2.503813816e-03, -2.127308728e-07, -7.689988780e-10,
                2.849677801e-13, -4.528198460e+04, -7.048279440e+00 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = { 1.176962419e+05, -1.788791477e+03,  8.291523190e+00,
               -9.223156780e-05,  4.863676880e-09, -1.891053312e-12,
                6.330036590e-16, -3.908350590e+04, -2.652669281e+01 }
   },
   { T_low  = 6000.0,
     T_high = 20000.0,
     coeffs = {-1.544423287e+09,  1.016847056e+06, -2.561405230e+02,
                3.369401080e-02, -2.181184337e-06,  6.991420840e-11,
               -8.842351500e-16, -8.043214510e+06,  2.254177493e+03 }
   },
   ref='Gurvich (1991) from cea2::thermo.inp'
}
-- CO2.viscosity = {
--    model = "Blottner",
--    parameters = { 
--       A_mu = -0.026654, B_mu = 1.107305, C_mu = -14.291274,
--       ref = "Table 4, ESA Radiation Test Case 8 (2014)"
--    }
-- }
CO2.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.51137258e+00, B=-0.22951321e+03, C=0.13710678e+05, D=0.27075538e+01},
      {T_low=1000.0, T_high=5000.0, A=0.63978285e+00, B=-0.42637076e+02, C=-0.15522605e+05, D=0.16628843e+01},
      {T_low=5000.0, T_high=10000.0, A=0.72150912e+00, B=0.75012895e+03, C=-0.11825507e+07, D=0.85493645e+00},
      ref = 'from CEA2::trans.inp which cites Boushehri et al (1987) and Svehla (1994)'
   }
}
CO2.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.48056568e+00, B=-0.50786720e+03, C=0.35088811e+05, D=0.36747794e+01},
      {T_low=1000.0, T_high=5000.0, A=0.69857277e+00, B=-0.11830477e+03, C=-0.50688859e+05, D=0.18650551e+01},
      {T_low=5000.0, T_high=10000.0, A=0.105183581e+01, B=-0.42555944e+04, C=0.14288688e+08, D=-0.88950473e+00},
      ref = 'from CEA2::trans.inp which cites Boushehri et al (1987) and Svehla (1994)'
   }
}
CO2.T_c = {
   value = 304.12,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
CO2.p_c = {
   value = 73.74e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
CO2.rho_c = {
   value = 466.5007, -- 10.60 mol/L
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'Ely, JF, JW Magee and WM Haynes (1987). Thermophysical properties for special high CO2 content mixtures',
}

-- Nonequilibrium data

CO2.species_type = "linear nonpolar polyatomic"
CO2.eps0 = {
   value = 2.695044416e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
CO2.sigma = {
   value = 3.941e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
CO2.s_0 = {
   value = 4857.70,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CO2.h_f = {
   value = -8941478.54,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO2.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CO2.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO2.theta_v = {
   value = 960.1,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Bending mode from Capitelli (2005) ESA STR-246'
}
CO2.electronic_levels = {
   -- n_levels = 22,
   n_levels = 5,
   ref = 'Capitelli (2005) ESA STR-246',
   -- NOTE: Have put the we[1] vibrational temperature first so it is used for VT exchange as done by Park.
   --       Previously we[0] was the primary vibrational temperature. 
   -- ==============================================================================================================================================
   --   n       Te         re       g       dzero      A0        B0          C0       sigma   sigma_rot   we[1]      we[0]      we[2]      we[3]
   -- ==============================================================================================================================================
   ilev_0  = {      0.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_1  = {  30000.00,  1.1621,  3,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_2  = {  33000.00,  1.1621,  6,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_3  = {  36000.00,  1.1621,  3,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_4  = {  46000.00,  1.1621,  2,      43984.06,  5.3,      0.42600,    1.246,   2,      2,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_5  = {  72480.00,  1.1621,  2,      43984.06,  5.3,      0.42600,    1.246,   2,      2,            667.3,  1225.00,   667.3,     2349.3 },
   ilev_6  = {  73100.00,  1.1621,  1,      43984.06,  5.3,      0.42600,    1.246,   2,      2,            667.3,  1225.00,   667.3,     2349.3 },
   ilev_7  = {  85160.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_8  = {  85840.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_9  = {  88535.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1458.00,   667.3,     2349.3 },
   ilev_10 = {  89111.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_11 = {  91830.00,  1.1621,  1,      43984.06,  0.0,      0.39150,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_12 = {  92360.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_13 = {  96600.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_14 = {  99331.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_15 = { 100570.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_16 = { 100650.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_17 = { 100940.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_18 = { 127065.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_19 = { 127443.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_20 = { 130774.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   ilev_21 = { 132977.00,  1.1621,  1,      43984.06,  0.0,      0.39021,    0.000,   2,      0,            667.3,  1384.86,   667.3,     2349.3 },
   -- ==============================================================================================================================================
}

-- Real gas data

CO2.reference_state = {
   description = 'Reference temperature, internal energy and entropy',
   reference = 'Reynolds, WC (1979). Thermodynamic Properties in SI.',
   units = 'K, J/kg, J/kg.K',
   T0 = 216.54, -- Triple point temperature
   u0 = 3.2174105e+5,
   s0 = 2.1396056e+3,
}
CO2.Cp0_coeffs = {
   description = 'Coefficients for Cp0 polynomial of form T^(-2) to T^4.',
   reference = 'Reynolds, WC (1979). Thermodynamic Properties in SI.',
   note = 'Cv0 coefficients converted to Cp0 model using Cp0 = Cv0 + R',
   units = 'J/kg.K', -- Resulting unit when using these coefficients in Cp0 equation.
   T_low  = 50.0, -- Range of validity.
   T_high = 1500.0,
   G = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
        8.726361e+03,
        3.729294e+02, -- = 1.840040e+02 + 188.92535
        1.914025e+00,
       -1.667825e-03,
        7.305950e-07,
       -1.255290e-10}
}
CO2.Bender_EOS_coeffs = {
   description = 'Coefficients for Bender equation of state',
   reference = 'Reynolds, WC (1979). Thermodynamic Properties in SI. Equation P-3',
   units = 'Pa', -- Resulting unit when using these coefficients in p-rho-T equation.
   A = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
        2.2488558e-01,  -1.3717965e+02,  -1.4430214e+04,
       -2.9630491e+06,  -2.0606039e+08,   4.5554393e-05,
        7.7042840e-02,   4.0602371e+01,   4.0029509e-07,
       -3.9436077e-04,   1.2115286e-10,   1.0783386e-07,
        4.3962336e-11,  -3.6505545e+04,   1.9490511e+07,
       -2.9186718e+09,   2.4358627e-02,  -3.7546530e+01,
        1.1898141e+04,
        5.0e-6} -- Final element is "gamma", roughly 1/(rho_c^2) as used in the MBWR EOS.
}
CO2.MBWR_EOS_coeffs = {
   description = 'Coefficients for MBWR equation of state',
   reference = 'Ely, JF, JW Magee and WM Haynes (1987). Thermophysical properties for special high CO2 content mixtures',
   units = 'bar', -- Resulting unit when using these coefficients in p-rho-T equation.
   b = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
       -0.981851065838e-02,   0.995062267309e+00,  -0.228380160313e+02,
        0.281827634529e+04,  -0.347001262699e+06,   0.394706709102e-03,
       -0.325550000110e+00,   0.484320083063e+01,  -0.352181542995e+06,
       -0.324053603343e-04,   0.468596684665e-01,  -0.754547012075e+01,
       -0.381894354016e-04,  -0.442192933859e-01,   0.516925168095e+02,
        0.212450985237e-02,  -0.261009474785e-04,  -0.888533388977e-01,
        0.155226179403e-02,   0.415091004940e+06,  -0.110173967489e+08,
        0.291990583344e+04,   0.143254606508e+08,   0.108574207533e+02,
       -0.247799657039e+03,   0.199293590763e-01,   0.102749908059e+03,
        0.377618865158e-04,  -0.332276512346e-02,   0.179196707121e-07,
        0.945076627807e-05,  -0.123400943061e-02}
}
