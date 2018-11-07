-- Collater: Daniel F Potter
-- Date: 29-Apr-2010

C2_plus = {}
C2_plus.M = { 
   value = 24.0208514,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
C2_plus.atomic_constituents = {C=2}
C2_plus.charge = 1
C2_plus.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -9.913423840e+04,  1.347170609e+03, -3.476753160e+00,
    	        1.676429424e-02, -1.865908025e-05,  1.091134647e-08,
               -2.434913818e-12,  2.335454800e+05,  4.406644620e+01  }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { 3.836292810e+06, -6.242062450e+03,  2.779245639e+00,
    	       6.065865860e-03, -2.452799858e-06,  3.882942500e-10,
    	      -2.190639912e-14,  2.857447553e+05,  7.297383490e-01
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = { 4.992689800e+07, -2.017309121e+04,  6.342227540e+00,
    	       6.369922600e-04, -1.036760828e-07,  4.943272920e-12,
    	      -7.988804260e-17,  4.120857610e+05, -2.112967169e+01
    }
  },

  ref="from CEA2::thermo.inp"
}

-- Nonequilibrium data

C2_plus.species_type = "nonpolar diatomic"
C2_plus.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C2_plus.h_f = {
   value = 83459803.26,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2_plus.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C2_plus.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2_plus.electronic_levels = {
   n_levels = 2,
   ref = 'Capitelli (2005)',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.301,   4,  42908.35,  1350.0,    0.0000,   0.000E+00,  0.000E+00,  1.65900,   0.000E+00,  1.000E-05,  0.000E+00,  0.000E+00,  0,  0 },
   ilev_1  = {  40143.00,  1.306,   2,   2765.35,  1340.0,    0.0000,   0.000E+00,  0.000E+00,  1.64800,  0.000E+00,  1.000E-05,  0.000E+00,  0.000E+00,  0,  0 }
   -- ===========================================================================================================================================================
}
