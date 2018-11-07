-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009
-- Updated by Elise J. Fahy
-- Date: 26 Feb 2013

C2H2 = {}
C2H2.M = {
   value = 26.037280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H2.atomic_constituents = {C=2,H=2}
C2H2.charge = 0
C2H2.gamma = {
   value = 1.2321e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = { 1.598112089e+05,     -2.216644118e+03,     1.265707813e+01,
               -7.979651080e-03,      8.054992750e-06,    -2.433307673e-09,
               -7.529233180e-14,      3.712619060e+04,    -5.244338900e+01 }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = { 1.713847410e+06,     -5.929106660e+03,     1.236127943e+01,
                1.314186993e-04,     -1.362764431e-07,     2.712655786e-11,
               -1.302066204e-15,      6.266578970e+04,    -5.818960590e+01 }
   },
   ref='Gurvich (1991) from cea2::thermo.inp'
}
C2H2.T_c = {
   value = 308.30,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C2H2.p_c = {
   value = 61.14e05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C2H2.rho_c = {  
   value = 231.0,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'Compressed Gas Association (1999). Handbook of Compressed Gases, 4th Ed.',
}

-- Nonequilibrium data

C2H2.species_type = "linear nonpolar polyatomic"
C2H2.eps0 = {
   value = 2.8842e-21,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = "Chemkin Transport User Manual (2001)",
}
C2H2.sigma = {
   value = 4.10e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = "Chemkin Transport User Manual (2001)",
}
C2H2.h_f = {
   value = 8764356.338,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C2H2.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'none'
}
C2H2.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'NA'
}
C2H2.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C2H2.theta_v = {
   value = 875.8,
   units = 'K',
   description = 'Characteristic vibrational temperature',
   reference = 'we[1] from electronic level table below'
}
C2H2.electronic_levels = {
   n_levels = 2,
   ref = "various - see notes below",
   -- NOTE: Have put the we[1] vibrational temperature first so it is used for VT exchange as done for CO2.
   --       Table has original notation in top row, and alternate notation underneath, since sources such as 	
   --       Capitelli use the alternate notation.
   -- NOTES ON SELECTED VALUES:
   -- Te for first excited state from Ganot et al (2003), "Non-adiabatic dissociation of rovibrationally excited acetylene".
   -- re is the average of the three bond lengths: C-C and two C-H bonds
   -- re for ground state from Herman et al (2003), "NIST vibration spectral database on acetylene".
   -- re for first excited state from Tobiason et al (1993), "Normal modes analysis of state acetylene...".
   -- g from Ochkin (2009), "Appendix A: Statistical Weights and Statistical Sums".
   -- dzero from Partridge & Bauschlicher Jr (1995), "The dissociation energies of CH4 and C2H2 revisited".
   -- B0 from Herman et al (2003), "NIST vibration spectral database on acetylene".
   -- Ground state frequencies from Zou & Bowman (2003), "A new ab initio potential energy surface describing acetylene/vinylidene isomerisation".
   -- Excited state frequencies from Tobiason et al (1993), "Normal modes analysis of state acetylene...".
   -- Ground state is linear and has two doubly degenerate vibrational modes, whereas excited state is bent
   -- (i.e. non-linear) and has six singly degenerate modes only.
   -- Only two energy levels have been included at present because the temperatures in the boundary layer
   -- should not be high enough to warrant the inclusion of higher levels. 
   -- Sigma has been inferred from Ochkin (2009) and Capitelli (2005). 
   -- A0, C0 and sigma_rot have been set to zero because at present, they are not required for any calculations
   -- (cf. Capitelli, 2005). However, there is reason to believe that sigma_rot = 1 (cf. Fernandez-Ramos et al, 
   -- 2007) and Ic > Ib > Ia (cf. Lovas - NIST Microwave Spectral Tables for Triatomic Molecules). 
   -- Further investigation may be required for these values, if they are needed in the code.
   --
   -- =============================================================================================================================================================
   --   n          Te      re       g      dzero      A0       B0      C0   sigma   sigma_rot   we[0]    we[1]     we[2]     we[3]      we[4]    we[5]    we[6]
   --                             (p_i)   (E_diss)                                              
   -- =============================================================================================================================================================
   ilev_0 = {     0.0,    1.108,    2,    46061.1,    0.0,    1.18,    0.0,    2,      0,      3397.12,  608.73,   608.73,    3316.86,  1981.80,  729.08,  729.08},
   ilev_1 = {42197.56,    1.282,    1,    46061.1,    0.0,    1.18,    0.0,    2,      0,      1420.0,   3004.0,   1064.0,    765.0,    914.0,    785.0}
   -- =============================================================================================================================================================
}

-- Real gas data  - not required at present for C2H2 because we are using thermally perfect cases only.

