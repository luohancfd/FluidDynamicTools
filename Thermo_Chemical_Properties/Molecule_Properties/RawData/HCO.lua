-- Author: Elise J. Fahy
-- Date: 28 Feb 2013

HCO = {}
HCO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
HCO.atomic_constituents = {C=1,H=1,O=1}
HCO.charge = 0
HCO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
HCO.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {-1.189851887e+04,    2.151536111e+02,     2.730224028e+00,
               1.806516108e-03,    4.984300570e-06,    -5.814567920e-09,
               1.869689894e-12,    2.905755640e+03,     1.136772540e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { 6.949606120e+05,   -3.656223380e+03,    9.604731170e+00,
              -1.117129278e-03,    2.875328019e-07,   -3.626247740e-11,
               1.808329595e-15,    2.543704440e+04,   -3.582473720e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
HCO.T_c = {
   value = 0.0,
   units = 'K',
   description = 'critical temperature',
   reference = 'none'
}
HCO.p_c = {
   value = 0.0,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'none'
}
HCO.rho_c = {  
   value = 0.0,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'none',
}

-- Nonequilibrium data

HCO.species_type = "nonlinear polar polyatomic"
HCO.eps0 = {
   value = 6.8724e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = "Chemkin Transport User Manual (2001)",
}
HCO.sigma = {
   value = 3.59e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = "Chemkin Transport User Manual (2001)",
}
HCO.h_f = {
   value = 1461085.931,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
HCO.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'none'
}
HCO.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
HCO.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
HCO.theta_v = {
   value = 1557.6,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'we[1] from electronic level table below'
}
HCO.electronic_levels = {
   n_levels = 3,
   ref = "Adamson (1994), The Spectroscopy of the Formyl Radical",
   --       Table has original notation in top row, and alternate notation underneath, since sources such as 	
   --       Capitelli use the alternate notation.
   -- NOTES ON SELECTED VALUES:
   -- re is the average of the two bond lengths: H-C and C-O.
   -- g from Ochkin (2009), "Appendix A: Statistical Weights and Statistical Sums".
   -- All states have three singly degenerate modes. Ground state is linear, excited states are bent (non-linear).
   -- Sigma from Ochkin (2009),"Appendix A: Statistical Weights and Statistical Sums".
   -- sigma_rot have been set to zero because at present, it is not required for any calculations
   -- (cf. Capitelli, 2005). However, there is reason to believe that sigma_rot = 1 (cf. Fernandez-Ramos et al, 2007) 
   -- For n=0 and 2, A0 and C0 are included because the molecule is bent, but for n=1, they are not required
   -- because the molecule is linear. 
   --
   -- ============================================================================================================================
   --   n          Te       re      g      dzero       A0       B0          C0    sigma  sigma_rot   we[0]     we[1]      we[2] 
   --                             (p_i)   (E_diss)                                              
   -- ============================================================================================================================
   ilev_0 =  {    0.0,    1.150,    2,    5058.0,    24.29,    1.4937,    1.3988,    1,     0,     2488.0,    1080.7,    1861.0},
   ilev_1 =  { 9297.0,    1.125,    2,    5058.0,      0.0,    1.3389,       0.0,    1,     0,     3319.0,     805.4,    1812.2},
   ilev_2 =  {38695.0,    1.260,    2,    5058.0,    16.04,     1.192,     1.108,    1,     0,     2597.0,    1382.0,    1065.0},
   -- ============================================================================================================================
}

-- Real gas data  - not required at present for HCN because we are using thermally perfect cases only.

