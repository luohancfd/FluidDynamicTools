-- Author: Daniel F. Potter
-- Date: 24-Sept-2009

-- Diatomic nitrogen cation

N2_plus = {}
N2_plus.M = { 
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
N2_plus.atomic_constituents = {N=2}
N2_plus.charge = 1
N2_plus.gamma = {
   value = 1.399, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
N2_plus.viscosity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.62526577e+00, B=-0.31779652e+02, C=-0.16407983e+04, D=0.17454992e+01 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.87395209e+00, B=0.56152222e+03, C=-0.17394809e+06, D=-0.39335958e+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.88503551e+00, B=0.90902171e+03, C=-0.73129061e+06, D=-0.53503838e+00 }
  },
  ref="cea2 data for N2"
}
N2_plus.thermal_conductivity = {
  model = "CEA",
  parameters = {
     {T_low = 200.0, T_high = 1000.0, A = 0.85439436e+00, B=0.10573224e+03, C=-0.12347848e+05, D=0.47793128e+00 },
     {T_low = 1000.0, T_high = 5000.0, A = 0.88407146e+00, B=0.13357293e+03, C=-0.11429640e+05, D=0.24417019e+00 },
     {T_low = 5000.0, T_high = 15000.0, A = 0.24176185e+01, B = 0.80477749e+04, C = 0.31055802e+07, D=-0.14517761e+02 }
  },
  ref="cea2 data for N2"
}
N2_plus.CEA_coeffs = {
  { T_low = 298.150,
    T_high = 1000.0,
    coeffs = { -3.474047470e+04,  2.696222703e+02,  3.164916370e+00,
               -2.132239781e-03,  6.730476400e-06, -5.637304970e-09, 
                1.621756000e-12,  1.790004424e+05,  6.832974166e+00
    }
  },
  { T_low = 1000.0,
    T_high = 6000.0,
    coeffs = { -2.845599002e+06,  7.058893030e+03, -2.884886385e+00,
                3.068677059e-03, -4.361652310e-07,  2.102514545e-11, 
                5.411996470e-16,  1.340388483e+05,  5.090897022e+01
    }
  },
  { T_low = 6000.0,
    T_high = 20000.0,
    coeffs = { -3.712829770e+08,  3.139287234e+05, -9.603518050e+01, 
                1.571193286e-02, -1.175065525e-06,  4.144441230e-11, 
               -5.621893090e-16, -2.217361867e+06,  8.436270947e+02
    }
  }
}

-- thermal nonequilibrium data

N2_plus.species_type = "nonpolar diatomic"
N2_plus.s_0 = {
   value = 7055.69,
   units = 'J/kg-k',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
N2_plus.h_f = {
   value = 5.3886282e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_plus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
N2_plus.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_plus.electronic_levels = {
   n_levels = 5,
   ref = 'Spradian07::diatom.dat (confirmed with NIST, slightly different data for upper levels)',
   -- ===========================================================================================================================================================
   --    n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {     0.00,  1.1164,  2,  70300.00,  2207.220,  16.2260,  4.000E-03, -6.100E-03,  1.93171,  1.882E-02,  6.100E-06,  0.000E+00,  0.000E+00,  0,  2 },
   ilev_1  = {  9167.46,  1.1749,  4,  61280.00,  1903.700,  15.1110,  1.120E-02, -2.700E-04,  1.74450,  1.883E-02,  5.600E-06,  1.800E-07, -7.462E+01,  1,  2 },
   ilev_2  = { 25461.11,  1.0742,  2,  44730.00,  2421.140,  24.0700, -3.000E-01, -6.670E-02,  2.08507,  2.120E-02,  6.170E-06,  0.000E+00,  0.000E+00,  0,  2 },
   ilev_3  = { 51663.20,  1.4710,  4,  18630.00,   907.710,  11.9100,  1.600E-02,  0.000E+00,  1.11300,  2.000E-02,  5.000E-06,  0.000E+00, -1.650E+01,  1,  2 },
   ilev_4  = { 64609.03,  1.2628,  2,  24990.00,  2069.400,   8.3000, -6.300E-01,  1.300E-02,  1.50980,  2.000E-02,  4.000E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
