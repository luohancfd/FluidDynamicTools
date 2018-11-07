-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

NCO = {}
NCO.M = {
   value = 42.016800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
NCO.atomic_constituents = {C=1,N=1,O=1}
NCO.charge = 0
NCO.gamma = {
   value = 1.2609e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
-- NCO.CEA_coeffs = {
--    { T_low  = 200.0,
--      T_high = 1000.0,
--      coeffs = {0.000000000e00, 0.000000000e00, 2.826930800e00, 8.805168800e-03, -8.386613400e-06, 4.801696400e-09, -1.331359500e-12, 1.468247700e04, 9.550464600e00, }
--    },
--    { T_low  = 1000.0,
--      T_high = 6000.0,
--      coeffs = {0.000000000e00, 0.000000000e00, 5.152184500e00, 2.305176100e-03, -8.803315300e-07, 1.478909800e-10, -9.097799600e-15, 1.400412300e04, -2.544266000e00, }
--    },
--    ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
-- }
NCO.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 1.136503036e+04, -2.444613367e+02,  4.671376100e+00,
     	        2.309387548e-03,  2.798649599e-06, -4.546357380e-09,
                1.692880931e-12,  1.577649188e+04, -2.171476903e-01 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = { 1.089445289e+05, -1.735459316e+03,  8.655610330e+00,
     	       -4.053229260e-04,  7.599716410e-08, -7.253804150e-12,
                3.244872410e-16,  2.365792776e+04, -2.619532970e+01 }
   },
   ref='East (1993) and Jacox (1998) from cea2::thermo.inp'
}

-- Nonequilibrium data

NCO.species_type = "linear polar polyatomic"
NCO.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NA'
}
NCO.h_f = {
   value = 3137964.84,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
NCO.I = {
   value = 270027882.67,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
NCO.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
NCO.theta_v = {
   value = 766.4,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Bending mode from Capitelli (2005) ESA STR-246'
}
NCO.electronic_levels = {
   n_levels = 4,
   ref = 'Capitelli (2005) ESA STR-246',
   -- ==========================================================================================================================================
   --   n       Te         re       g       dzero      A0        B0          C0       sigma   sigma_rot   we[0]      we[1]      we[1]      we[2]
   -- ==========================================================================================================================================
   ilev_0  = {      0.00,  2.4080,  2,      33875.29,  0.000,    0.38940,    0.000,   1,      0,          1275.0,    532.7,     532.7,     1922.0 },
   ilev_1  = {     95.20,  2.4080,  2,      33875.29,  0.000,    0.38940,    0.000,   1,      0,          1275.0,    532.7,     532.7,     1922.0 },
   ilev_2  = {  22802.00,  2.3690,  2,      33875.29,  0.000,    0.40210,    0.000,   1,      0,          2338.0,    680.8,     680.8,     1289.3 },
   ilev_3  = {  31811.00,  2.4500,  4,      33875.29,  0.000,    0.37650,    0.000,   1,      0,          2303.0,    532.7,     532.7,     1047.0 },
}
