-- Collater: Rowan J. Gollan
-- Date: 21-Jun-2009
-- Updated by Elise J. Fahy
-- Date: 26 Feb 2013

HCN = {}
HCN.M = {
   value = 27.0253400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HCN.atomic_constituents = {C=1,H=1,N=1}
HCN.charge = 0
HCN.gamma = {
   value = 1.301,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HCN.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  9.098286930e+04, -1.238657512e+03,  8.721307870e+00,
               -6.528242940e-03,  8.872700830e-06, -4.808886670e-09,
                9.317898500e-13,  2.098915450e+04, -2.746678076e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.236889278e+06, -4.446732410e+03,  9.738874850e+00,
               -5.855182640e-04,  1.072791440e-07, -1.013313244e-11,
                3.348247980e-16,  4.221513770e+04, -4.005774072e+01
    }
  },
  ref="from CEA2::thermo.inp"
}

HCN.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.94863717e+00, B=-0.14891490e+03, C=0.15258721e+05, D=-0.72592817e+00},
      {T_low=1000.0, T_high=5000.0, A=0.57370725e+00, B=-0.85239973e+03, C=0.17953641e+06, D=0.24032031e+01},
      ref = 'from CEA2::trans.inp which cites Zeleznik & Svehla (1970), Svehla (1994)'
   }
}
HCN.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.11749061e+01, B=-0.19100307e+03, C=0.15714065e+05, D=-0.13488014e+01},
      {T_low=1000.0, T_high=5000.0, A=0.50543688e+00, B=-0.13891056e+04, C=0.28003144e+06, D=0.42095130e+01},
      ref = 'from CEA2::trans.inp which cites Zeleznik & Svehla (1970), Svehla (1994)'
   }
}
HCN.T_c = {
   value = 0.0,
   units = 'K',
   description = 'critical temperature',
   reference = 'none'
}
HCN.p_c = {
   value = 0.0,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'none'
}
HCN.rho_c = {  
   value = 0.0,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'none',
}

-- Nonequilibrium data

HCN.species_type = "linear nonpolar polyatomic"
HCN.eps0 = {
   value = 7.7542e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = "Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5",
}
HCN.sigma = {
   value = 3.630e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = "Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5",
}
HCN.h_f = {
   value = 4924358.398,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
HCN.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'none'
}
HCN.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
HCN.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
HCN.theta_v = {
   value = 1009.1,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'we[1] from electronic level table below'
}
HCN.electronic_levels = {
   n_levels = 13,
   ref = "Nayak et al (2005), Theoretical study on the excited states of HCN",
   -- NOTE: Have put the we[1] vibrational temperature first so it is used for VT exchange as done for CO2.
   --       Table has original notation in top row, and alternate notation underneath, since sources such as 	
   --       Capitelli use the alternate notation.
   -- NOTES ON SELECTED VALUES:
   -- re is the average of the two bond lengths: H-C and C-N.
   -- g from Ochkin (2009), "Appendix A: Statistical Weights and Statistical Sums".
   -- All states have three singly degenerate modes. Ground state is linear, excited states are bent (non-linear).
   -- Sigma from Ochkin (2009),"Appendix A: Statistical Weights and Statistical Sums".
   -- A0, C0 and sigma_rot have been set to zero because at present, they are not required for any calculations
   -- (cf. Capitelli, 2005). However, there is reason to believe that sigma_rot = 1 (cf. Fernandez-Ramos et al, 
   -- 2007) and Ic > Ib > Ia (cf. Lovas - NIST Microwave Spectral Tables for Triatomic Molecules). 
   -- Further investigation may be required for these values, if they are needed in the code.
   --
   -- ==========================================================================================================================
   --   n          Te         re      g       dzero        A0      B0      C0   sigma  sigma_rot  we[0]     we[1]     we[2] 
   --                               (p_i)    (E_diss)                                              
   -- ==========================================================================================================================
   ilev_0 =  {     0.0,     1.110,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3351.0,    2101.0,    702.0},
   ilev_1 =  {35971.80,     1.219,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2818.0,    1525.0,    981.0},
   ilev_2 =  {43714.61,      1.22,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2605.0,    1413.0,    948.0},
   ilev_3 =  {50650.88,    1.2145,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3121.0,    1452.0,    875.0},
   ilev_4 =  {52263.97,    1.2185,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2606.0,    1413.0,    948.0},
   ilev_5 =  {54602.94,    1.2365,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     1987.0,    1433.0,    499.0},
   ilev_6 =  {57103.22,     1.218,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3072.0,    1439.0,    966.0},
   ilev_7 =  {57345.18,     1.205,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2764.0,    1814.0,   1032.0},
   ilev_8 =  {61377.90,    1.1615,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3356.0,    1759.0,    893.0},
   ilev_9 =  {61781.17,    1.2185,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3095.0,    1414.0,    935.0},
   ilev_10 = {65652.57,     1.208,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2524.0,    1595.0,    760.0},
   ilev_11 = {68798.09,      1.34,    3,     45327.70,    0.0,    1.48,    0.0,    1,    0,     2782.0,    1462.0,   1003.0},
   ilev_12 = {71620.99,     1.151,    1,     45327.70,    0.0,    1.48,    0.0,    1,    0,     3602.0,    1762.0,   1011.0},
   -- ==========================================================================================================================
}

-- Real gas data  - not required at present for HCN because we are using thermally perfect cases only.

