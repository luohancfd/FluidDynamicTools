-- Collater: Rowan J. Gollan
-- Date: 20-Jun-2009

C2 = {}
C2.M = { 
   value = 24.0214000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
C2.atomic_constituents = {C=2}
C2.charge = 0
C2.gamma = { 
   value = 1.236,
   units = 'non-dimensional',
    description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
C2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  5.559634510e+05, -9.980126440e+03,  6.681620370e+01,
               -1.743432724e-01,  2.448523051e-04, -1.703467580e-07,
                4.684527730e-11,  1.445869634e+05, -3.448229700e+02
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -9.689267930e+05,  3.561092990e+03, -5.064138930e-01,
                2.945154879e-03, -7.139441190e-07,  8.670657250e-11,
               -4.076906810e-15,  7.681796830e+04,  3.339985240e+01
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = {  6.315145920e+06,  1.365420661e+04, -3.996903670e+00,
                1.937561376e-03, -1.584446580e-07,  5.520861660e-12,
               -7.253735340e-17,  9.387024990e+03,  6.614329920e+01
    }
  },

  ref="from CEA2::thermo.inp"
}
--C2.viscosity = {
--    model = "Blottner",
--    parameters = { 
--       A_mu = -0.0031, B_mu = 0.6920, C_mu = -12.6127,
--       ref = "Blottner (1971)"
--    }
-- }

-- Nonequilibrium data

C2.species_type = "nonpolar diatomic"
C2.eps0 = {
   value = 2.693663758e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
C2.sigma = {
   value = 3.913e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
C2.s_0 = {
   value = 8300.51,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C2.h_f = {
   value = 34571562.11,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2.I = {
   value = 45829859.25,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C2.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2.electronic_levels = {
   -- n_levels = 11,
   n_levels = 5,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.2425,  1,  50104.00,  1854.710,  13.3400, -1.720E-01,  0.000E+00,  1.81984,  1.765E-02,  6.920E-06,  8.100E-08,  0.000E+00,  0,  1 },
   ilev_1  = {    716.24,  1.3119,  6,  49500.00,  1641.350,  11.6700,  0.000E+00,  0.000E+00,  1.63246,  1.661E-02,  6.440E-06,  0.000E+00, -1.525E+01,  1,  3 },
   ilev_2  = {   6434.30,  1.3693,  3,  43860.00,  1470.450,  11.1900,  2.800E-02,  0.000E+00,  1.49852,  1.634E-02,  6.220E-06,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_3  = {   8391.00,  1.3184,  2,  41840.00,  1608.350,  12.0780, -1.000E-02,  0.000E+00,  1.61634,  1.686E-02,  6.440E-06,  3.600E-08,  0.000E+00,  1,  1 },
   ilev_4  = {  13312.10,  1.2300,  3,  36740.00,  1961.600,  13.7000,  0.000E+00,  0.000E+00,  1.87000,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0,  3 },
   ilev_5  = {  20022.50,  1.2661,  6,  30120.00,  1788.220,  16.4400, -5.067E-01,  0.000E+00,  1.75270,  1.608E-02,  6.740E-06,  1.030E-07, -1.690E+01,  1,  3 },
   ilev_6  = {  34261.30,  1.2552,  2,  15870.00,  1809.100,  15.8100,  0.000E+00,  0.000E+00,  1.78340,  1.800E-02,  6.800E-06,  0.000E+00,  0.000E+00,  1,  1 },
   ilev_7  = {  40796.65,  1.5351,  6,  19880.00,  1106.560,  39.2600,  0.000E+00,  0.000E+00,  1.19220,  2.420E-02,  6.300E-06,  2.900E-07,  0.000E+00,  1,  3 },
   ilev_8  = {  43239.41,  1.2380,  1,  20000.00,  1829.569,  13.9400,  0.000E+00,  0.000E+00,  1.83320,  1.960E-02,  7.362E-06,  0.000E+00,  0.000E+00,  0,  1 },
   ilev_9  = {  55034.66,  1.2529,  1,  20000.00,  1671.499,  40.0199,  2.480E-01,  0.000E+00,  1.78970,  3.870E-02,  8.207E-06,  0.000E+00,  0.000E+00,  0,  1 },
   ilev_10 = {  71045.75,  1.3929,  3,  20000.00,  1360.499,  14.8000,  0.000E+00,  0.000E+00,  1.44800,  4.000E-02,  6.561E-06,  0.000E+00,  0.000E+00,  0,  3 },
   -- ===========================================================================================================================================================
}
