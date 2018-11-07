-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009
-- Updated by Elise J. Fahy
-- Date: 20 Feb 2013

C2H = {}
C2H.M = {
   value = 25.029340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H.atomic_constituents = {C=2,H=1}
C2H.charge = 0
C2H.gamma = {
   value = 1.2465e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 1.343669487e+04,   -5.067970720e+02,    7.772107410e+00,
               -6.512339820e-03,    1.030117855e-05,   -5.880147670e-09,
                1.226901861e-12,    6.892269990e+04,   -1.871881626e+01 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = { 3.922334570e+06,    -1.204751703e+04,     1.756172920e+01,
               -3.655442940e-03,     6.987685430e-07,    -6.825162010e-11,
                2.719262793e-15,     1.433266627e+05,    -9.561634380e+01 }
   },
   ref="Ervin (1990), Jacox (1998), Peric (1990), and Kanamori (1988) from cea2::thermo.inp"
}
C2H.T_c = {
   value = 0.0,
   units = 'K',
   description = 'critical temperature',
   reference = 'none',
}
C2H.p_c = { 
   value = 0.0,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'none',
}
C2H.rho_c = {  
   value = 0.0,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'none',
}

-- Nonequilibrium data

C2H.species_type = "linear nonpolar polyatomic"
C2H.eps0 = {
   value = 7.4279e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = "Park (2001) Chemical-Kinetic Parameters of Hyperbolic Earth Entry",
}
C2H.sigma = {
   value = 3.243e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = "Park (2001) Chemical-Kinetic Parameters of Hyperbolic Earth Entry",
}
C2H.h_f = {
   value = 22621470.72,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2H.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'none'
}
C2H.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
C2H.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2H.theta_v = {
   value = 559.7,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'Bending mode from Cui & Morokuma (2011) ESA STR-246'
}
C2H.electronic_levels = {
   n_levels = 2,
   ref = "Cui & Morokuma (2011) Studies on Excited States of C2H",
   -- NOTE: Have put the we[1] vibrational temperature first so it is used for VT exchange as done for CO2.
   --       Table has original notation in top row, and alternate notation underneath, since sources such as 	
   --       Capitelli use the alternate notation.
   -- NOTES ON SELECTED VALUES:
   -- Only two energy levels have been included at present because the temperatures in the boundary layer
   -- should not be high enough to warrant the inclusion of higher levels. 
   -- Sigma has been inferred from Ochkin (2009) and Capitelli (2005). 
   -- B0 was selected from Tarroni & Carter (2003) based on Te and the electronic state configuration.
   -- A0, C0 and sigma_rot have been set to zero because at present, they are not required for any calculations
   -- (cf. Capitelli, 2005). However, there is reason to believe that sigma_rot = 1 (cf. Fernandez-Ramos et al, 
   -- 2007) and Ic > Ib > Ia (cf. Lovas - NIST Microwave Spectral Tables for Triatomic Molecules). 
   -- Further investigation may be required for these values, if they are needed in the code.
-- =======================================================================================================================================================
   --   n       Te      re       g        dzero        A0        B0          C0    sigma   sigma_rot     we[1]      we[0]     we[2]      we[3]
   --                          (p_i)    (E_diff)                                                         nu[1]      nu[0]     nu[2]      nu[3]
   -- ====================================================================================================================================================
   ilev_0 = {   0.0,  1.217,     2,      39388.0,     0.0,     1.4117,      0.0,     1,      0,          389.0,     3612.0,   389.0,     1848.0},
   ilev_1 = {3692.6,  1.288,     4,      39388.0,     0.0,     1.4117,      0.0,     1,      0,          571.0,     3424.0,   571.0,     1795.0}
   -- ====================================================================================================================================================
}

-- Real gas data  - not required at present for C2H because we are using thermally perfect cases only.

