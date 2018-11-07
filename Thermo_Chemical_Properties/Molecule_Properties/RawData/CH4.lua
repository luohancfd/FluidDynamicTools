-- Collater: Rowan J. Gollan
-- Date: 17-Apr-2009

CH4 = {}
CH4.M = {
   value = 16.0424600e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH4.atomic_constituents = {C=1,H=4}
CH4.charge = 0
CH4.gamma = {
   value = 1.303,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH4.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   }
}
CH4.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   }
}
CH4.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -1.766850998e+05,  2.786181020e+03, -1.202577850e+01,
	        3.917619290e-02, -3.619054430e-05,  2.026853043e-08,
               -4.976705490e-12, -2.331314360e+04,  8.904322750e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  3.730042760e+06, -1.383501485e+04,  2.049107091e+01,
	       -1.961974759e-03,  4.727313040e-07, -3.728814690e-11,
                1.623737207e-15,  7.532066910e+04, -1.219124889e+02
    }
  },
  ref="from CEA2::thermo.inp"
}
-- CH4.CEA_coeffs = {
--    { T_low  = 200.0,
--      T_high = 1000.0,
--      coeffs = {0.00000000e+00, 0.00000000e+00, 5.14987613e+00, -1.36709788e-02, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -1.02466476e+04, -4.64130376e+00},
--    },
--    { T_low  = 1000.0,
--      T_high = 3500.0,
--      coeffs = {0.00000000e+00, 0.00000000e+00, 7.48514950e-02, 1.33909467e-02, -5.73285809e-06, 1.22292535e-09, -1.01815230e-13, -9.46834459e+03, 1.84373180e+01},
-- },
--    ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
-- }
CH4.T_c = {
   value = 190.56,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
CH4.p_c = {
   value = 45.99e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}

-- noneq data

CH4.species_type = "nonlinear nonpolar polyatomic"
CH4.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CH4.h_f = {
   value = -4650159.64,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CH4.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CH4.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CH4.theta_v = {
   value = 1879.0,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Conversion from NIST WebBook'
}
CH4.electronic_levels = {
   n_levels = 1,
   ref = 'various',
   -- NOTES ON VALUES:
   -- All values from NIST CCCBDB unless otherwise stated.
   -- re is the average of the bond lengths: H-C.
   -- g from Ochkin (2009), "Appendix A: Statistical Weights and Statistical Sums".
   -- dzero from NIST CCCBDB, converted to kJ/kg.
   -- Sigma from Ochkin (2009),"Appendix A: Statistical Weights and Statistical Sums", point group from NIST CCCBDB
   -- A0, C0 and sigma_rot have been set to zero because at present, they are not required for any calculations
   -- (cf. Capitelli, 2005).
   --
   -- ====================================================================================================================================
   --   n          Te         re      g       dzero        A0      B0      C0   sigma  sigma_rot  we[0]     we[1]     we[2]      we[2] 
   --                               (p_i)    (E_diss)                                              
   -- ====================================================================================================================================
   ilev_0  = {     0.00,    1.087,    1,    36720.9522294,     5.24120,   5.24120,   5.24120,   12,      0,    2917.0,    1534.0,   3019.0,    1306.0 },
   -- ====================================================================================================================================
}
