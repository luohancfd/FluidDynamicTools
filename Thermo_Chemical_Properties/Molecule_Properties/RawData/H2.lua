-- Author: Rowan J. Gollan
-- Date: 30-Oct-2008

-- diatomic hydrogen

H2 = {}
H2.M = { 
   value = 2.0158800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
H2.atomic_constituents = {H=2}
H2.charge = 0
H2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
-- Presently these values are oxygen values.
-- H2.d = { 
--    value = 3.433e-10,
--    units = 'm',
--    description = 'equivalent hard sphere diameter, based on L-J parameters',
--    reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
-- }
-- H2.viscosity = { 
--    model = "Sutherland",
--    parameters = { 
--       mu_ref = 1.919e-05, T_ref = 273.0, S = 139.0,
--       ref = "Table 1-2, White (2006)"
--    }
-- }
-- H2.thermal_conductivity = { 
--    model = "Sutherland",
--    parameters = { 
--       k_ref = 0.0244, T_ref = 273.0, S = 240.0,
--       ref = "Table 1-3, White (2006)"
--    }
-- }
-- End oxygen values? DFP
H2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  4.078323210e+04, -8.009186040e+02,  8.214702010e+00,
               -1.269714457e-02,  1.753605076e-05, -1.202860270e-08,
                3.368093490e-12,  2.682484665e+03, -3.043788844e+01
	    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  5.608128010e+05, -8.371504740e+02,  2.975364532e+00,
                1.252249124e-03, -3.740716190e-07,  5.936625200e-11,
               -3.606994100e-15,  5.339824410e+03, -2.202774769e+00
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = {  4.966884120e+08, -3.147547149e+05,  7.984121880e+01,
               -8.414789210e-03,  4.753248350e-07, -1.371873492e-11,
                1.605461756e-16,  2.488433516e+06, -6.695728110e+02
    }
  },
  ref="from CEA2::thermo.inp"
}
H2.viscosity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.74553182e+00, B=0.43555109e+02, C=-0.32579340e+04, D=0.13556243e+00 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.96730605e+00, B=0.67931897e+03, C=-0.21025179e+06, D=-0.18251697e+01 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.10126129e+01, B=0.14973739E+04, C=-0.14428484e+07, D=-0.23254928e+01 }
  },
  ref = "from CEA2::trans.inp which cites ASSAEL ET AL (1986)  SVEHLA (1994)"
}
H2.thermal_conductivity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.10059461E+01, B = 0.27951262E+03, C = -0.29792018E+05, D = 0.11996252E+01 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.10582450E+01, B = 0.24875372E+03, C = 0.11736907E+05, D = 0.82758695E+00 },
     {T_low = 5000.0, T_high = 15000.0, A = -0.22364420E+00, B = -0.69650442e+04, C = -0.77771313e+05, D = 0.13189369e+02 }
  },
  ref = "from CEA2::trans.inp which cites ASSAEL ET AL (1986)  SVEHLA (1994)"
}
H2.T_c = {
   value = 32.98,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
H2.p_c = {
   value = 12.93e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}

-- Nonequilibrium data

H2.species_type = "nonpolar diatomic"
H2.eps0 = {
   value = 8.24252826e-22,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
H2.sigma = {
   value = 2.827e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
H2.r0 = {
   value = 2.827e-10,
   units = 'm',
   description = 'Zero of the intermolecular potential',
   reference = 'Take the value for sigma, suggested by Thivet et al.'
}
H2.r_eq = {
   value = 0.7414e-10,
   units = 'm',
   description = 'Equilibrium intermolecular distance',
   reference = 'See ilev_0 data below'
}
H2.f_m = {
   value = 1.0,
   units = 'ND',
   description = 'Mass factor = ( M ( Ma^2 + Mb^2 ) / ( 2 Ma Mb ( Ma + Mb ) )',
   reference = 'Thivet et al (1991) Phys. Fluids A 3 (11)'
}
H2.mu = {
   value = 1.67372395859e-27,
   units = 'kg/particle',
   description = 'Reduced mass of constituent atoms',
   reference = 'See molecular weight for H'
}
H2.alpha = {
   value = 0.787,
   units = 'Angstrom^3',
   description = 'Polarizability',
   reference = 'http://cccbdb.nist.gov cites Olney et al. (1997)'
}
H2.s_0 = {
   value = 64825.29,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
H2.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
H2.I = {
   value = 738326998.13,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
H2.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
H2.electronic_levels = {
   n_levels = 12,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  0.7414,  1,  36117.40,  4401.213, 121.3360,  0.000E+00,  0.000E+00, 60.85300,  3.062E+00,  4.710E-02, -2.740E-03,  0.000E+00,  0,  1 },
   ilev_1  = {  91700.80,  1.2928,  1,  28200.00,  1358.090,  20.8880,  0.000E+00,  0.000E+00, 20.01540,  1.184E+00,  1.625E-02, -2.165E-03,  0.000E+00,  0,  1 },
   ilev_2  = {  95838.50,  1.0376,  6,  21319.60,  2466.890,  63.5100,  0.000E+00,  0.000E+00, 31.07000,  1.425E+00,  1.950E-02,  0.000E+00, -1.249E-01,  1,  3 },
   ilev_3  = {  95936.10,  0.9888,  3,  21125.00,  2664.830,  71.6500,  0.000E+00,  0.000E+00, 34.21600,  1.671E+00,  2.160E-02,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_4  = { 100057.80,  1.0328,  2,  19260.00,  2443.770,  69.5240,  0.000E+00,  0.000E+00, 31.36290,  1.665E+00,  2.230E-02, -7.400E-04,  0.000E+00,  1,  1 },
   ilev_5  = { 100082.30,  2.0000,  1,  17694.00,  1199.000,   0.0000,  0.000E+00,  0.000E+00,  6.24000,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  1 },
   ilev_6  = { 107774.70,  1.1070,  3,   9519.30,  2196.130,  65.8000,  0.000E+00,  0.000E+00, 27.30000,  1.515E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_7  = { 111642.80,  1.1192,  1,   7873.00,  2039.520,  83.4060,  0.000E+00,  0.000E+00, 26.70500,  2.781E+00,  1.200E-02,  0.000E+00,  0.000E+00,  0,  1 },
   ilev_8  = { 112700.30,  1.0496,  6,  19739.00,  2371.580,  66.2700,  0.000E+00,  0.000E+00, 30.36400,  1.545E+00,  1.910E-02,  0.000E+00,  2.910E-02,  1,  3 },
   ilev_9  = { 112854.40,  1.0500,  3,  19635.00,  2290.860, 105.4300,  0.000E+00,  0.000E+00, 30.00000,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_10 = { 112913.00,  1.0450,  3,   7717.00,  2268.730,   0.0000,  0.000E+00,  0.000E+00, 30.62000,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_11 = { 113888.70,  1.0508,  2,  18557.20,  2359.910,  68.8160,  0.000E+00,  0.000E+00, 30.29600,  1.420E+00,  2.010E-02,  0.000E+00,  0.000E+00,  1,  1 }
   -- ===========================================================================================================================================================
}

