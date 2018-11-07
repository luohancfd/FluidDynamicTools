-- Collater: Rowan J. Gollan
-- Date: 17-Apr-2009

CH3 = {}
CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH3.atomic_constituents = {C=1,H=3}
CH3.charge = 0
CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH3.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -2.876188806e+04,  5.093268660e+02,  2.002143949e-01,
	        1.363605829e-02, -1.433989346e-05,  1.013556725e-08,
               -3.027331936e-12,  1.408271825e+04,  2.022772791e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.760802663e+06, -9.336531170e+03,  1.487729606e+01,
               -1.439429774e-03,  2.444477951e-07, -2.224555778e-11,
                8.395065760e-16,  7.481809480e+04, -7.919682400e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
-- CH3.CEA_coeffs = {
--    { T_low  = 200.0,
--      T_high = 1000.0,
--      coeffs = {-2.87618881e+04, 5.09326866e+02, 2.00214395e-01, 1.36360583e-02, -1.43398935e-05, 1.01355673e-08, -3.02733194e-12, 1.40827182e+04, 2.02277279e+01}
--    },
--    { T_low  = 1000.0,
--      T_high = 6000.0,
--      coeffs = {2.76080266e+06, -9.33653117e+03, 1.48772961e+01, -1.43942977e-03, 2.44447795e-07, -2.22455578e-11, 8.39506576e-16, 7.48180948e+04, -7.91968240e+01}
--    },
--    ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
-- }
CH3.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
CH3.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}

-- noneq data

CH3.species_type = "nonlinear nonpolar polyatomic"
CH3.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CH3.h_f = {
   value = 9754753.73,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CH3.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CH3.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CH3.theta_v = {
   value = 872.6,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Conversion from NIST WebBook'
}
CH3.electronic_levels = {
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
   -- ==============================================================================================================================================
   --   n          Te      re      g       dzero        A0       B0           C0     sigma  sigma_rot     we[0]      we[1]      we[2]       we[2] 
   --                            (p_i)    (E_diss)                                              
   -- ==============================================================================================================================================
   ilev_0  = {  0.00,     1.079,    2,   38613.1357855,    9.57789,      9.57789,       4.74202,     6,      0,        3004.4,     606.5,    3160.8,     1396.0 },
   ilev_1  = {  46239.0,  1.079,    2,   38613.1357855,    9.57789,      9.57789,       4.74202,     6,      0,        2040.0,     606.5,    3160.8,     1396.0 },
   ilev_2  = {  59972.0,  1.079,    2,   38613.1357855,    9.57789,      9.57789,       4.74202,     6,      0,        2931.0,     1323.0,   3087.0,     1428.0 },
   -- ==============================================================================================================================================
}

