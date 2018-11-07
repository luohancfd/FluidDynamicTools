-- Collater: Rowan J. Gollan
-- Date: 20-Jun-2009

CH2 = {}
CH2.M = {
   value = 14.0265800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH2.atomic_constituents = {C=1,H=2}
CH2.charge = 0
CH2.gamma = {
   value = 1.311,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  3.218921730e+04, -2.877601815e+02,  4.203583820e+00,
		3.455405960e-03, -6.746193340e-06,  7.654571640e-09,
	       -2.870328419e-12,  4.733624710e+04, -2.143628603e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.550418031e+06, -7.971625390e+03,  1.228924487e+01,
	       -1.699122922e-03,  2.991728605e-07, -2.767007492e-11,
                1.051341740e-15,  9.642216890e+04, -6.094739910e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
CH2.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
CH2.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}

-- noneq data

CH2.species_type = "nonlinear nonpolar polyatomic"
CH2.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
CH2.h_f = {
   value = 27830341.89,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CH2.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
CH2.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CH2.theta_v = {
   value = 1946.0,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Conversion from NIST WebBook'
}
CH2.electronic_levels = {
   n_levels = 2,
   ref = 'various',
   -- NOTES ON VALUES:
   -- All values from NIST CCCBDB unless otherwise stated.
   -- re is the average of the bond lengths: H-C.
   -- g from Ochkin (2009), "Appendix A: Statistical Weights and Statistical Sums".
   -- dzero from NIST CCCBDB
   -- Sigma from Ochkin (2009),"Appendix A: Statistical Weights and Statistical Sums", point group from NIST CCCBDB
   -- A0, C0 and sigma_rot have been set to zero because at present, they are not required for any calculations
   -- (cf. Capitelli, 2005).
   --
   -- ========================================================================================================================
   --   n          Te      re      g       dzero        A0          B0         C0     sigma  sigma_rot  we[0]     we[1]     we[2] 
   --                            (p_i)    (E_diss)                                              
   -- ========================================================================================================================
   ilev_0 =  {     0.0,  1.085,    3,     35430.350,  73.05775,    8.42,    7.21927,    2,    0,     2806.0,     963.1,    3190.0},
   ilev_1 =  {  3147.0,  1.085,    1,     35430.350,  19.80540,   11.25,    7.23860,    2,    0,     2806.0,    1352.6,    2865.0},
   -- ========================================================================================================================
}
