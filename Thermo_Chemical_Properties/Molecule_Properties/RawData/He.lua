-- Author: Rowan J. Gollan
-- Date: 30-July-2008
--  26-Aug-2009: added CEA coeffs (BTO)
--  23-Feb-2010: more accurate mol weight, high temperature entropy curve fit

He = {}
He.M = {
   value = 4.0026020e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'cea2::thermo.inp'
}
He.atomic_constituents = {He=1}
He.charge = 0
He.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
He.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 6000.0,
     coeffs = { 0.000000000e+00,  0.000000000e+00,  2.500000000e+00, 
                0.000000000e+00,  0.000000000e+00,  0.000000000e+00, 
	        0.000000000e+00, -7.453750000e+02,  9.287239740e-01,
	      }
   },
   { T_low  = 6000.0,
     T_high = 20000.0,
     coeffs = { 3.396845420e+06, -2.194037652e+03,  3.080231878e+00,
     	       -8.068957550e-05,  6.252784910e-09, -2.574990067e-13,
     	        4.429960218e-18,  1.650518960e+04, -4.048814390e+00
	      }
   },
   ref='cea2::thermo.inp'
   -- NOTE: Data for first 2 ranges has been combined and is same as Chemkin data
}
He.d = {
   value = 2.576e-10,
   units = 'm',
   description = 'equivalent hard-sphere diameter, sigma from L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
-- He.viscosity = {
--    model = "Sutherland",
--    parameters = {
--       mu_ref = 1.865e-05, T_ref = 273.1, S = 72.9,
--       ref = [[
-- mu_ref, T_ref: Table 11, Chapman and Cowling (1970)
-- S: Table 15, Chapman and Cowling (1970)
--       ]]
--    }
-- }
He.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.75015944e+00, B=0.35763243e+02, C=-0.22121291e+04, D=0.92126352e+00},
      {T_low=1000.0, T_high=5000.0, A=0.83394166e+00, B=0.22082656e+03, C=-0.52852591e+05, D=0.20809361e+00},
      {T_low=5000.0, T_high=15000.0, A=0.86316349e+00, B=0.96205176e+03, C=-0.12498705e+07, D=-0.14115714e+00},
      ref = "from CEA2::trans.inp which cites Bich et al. (1990)"
   }
}
He.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.75007833e+00, B=0.36577987e+02, C=-0.23636600e+04, D=0.29766475e+01},
      {T_low=1000.0, T_high=5000.0, A=0.83319259e+00, B=0.22157417e+03, C=-0.53304530e+05, D=0.22684592e+01},
      {T_low=5000.0, T_high=15000.0, A=0.85920953e+00, B=0.89873206e+03, C=-0.11069262e+07, D=0.19535742e+01},
      ref = "from CEA2::trans.inp which cites Bich et al. (1990)"
   }
}
He.T_c = {
   value = 5.19,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
He.p_c = {
   value = 2.27e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}

-- Thermal nonequilibrium data

He.species_type = "monatomic"
He.eps0 = {
   value = 1.411032476e-22,
   units = 'J',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
He.sigma = {
   value = 2.551e-10,
   units = 'm',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
He.s_0 = {
   value = 31517.75,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
He.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
He.I = {
   value = 592696726.91,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
He.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
He.electronic_levels = {
   -- n_levels = 103,
   n_levels = 2,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.    n       E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0   =  { 1,        0.00,     1,   -1,    0,    0,    2 },
   ilev_1   =  { 1,   159855.97,     3,    0,    0,    1,    2 },
   ilev_2   =  { 1,   166277.44,     1,    0,    0,    0,    2 },
   ilev_3   =  { 1,   169086.91,     9,    1,    1,    1,    1 },
   ilev_4   =  { 1,   171134.90,     3,    1,    1,    0,    1 },
   ilev_5   =  { 1,   183236.79,     3,    0,    0,    1,    2 },
   ilev_6   =  { 1,   184864.83,     1,    0,    0,    0,    2 },
   ilev_7   =  { 1,   185564.60,     9,    1,    1,    1,    1 },
   ilev_8   =  { 1,   186101.55,    15,    2,    2,    1,    2 },
   ilev_9   =  { 1,   186104.96,     5,    2,    2,    0,    2 },
   ilev_10  =  { 1,   186209.36,     3,    1,    1,    0,    1 },
   ilev_11  =  { 1,   190298.11,     3,    0,    0,    1,    2 },
   ilev_12  =  { 1,   190940.23,     1,    0,    0,    0,    2 },
   ilev_13  =  { 1,   191217.06,     9,    1,    1,    1,    1 },
   ilev_14  =  { 1,   191444.48,    15,    2,    2,    1,    2 },
   ilev_15  =  { 1,   191446.45,     5,    2,    2,    0,    2 },
   ilev_16  =  { 1,   191451.88,    21,    3,    3,    1,    1 },
   ilev_17  =  { 1,   191451.90,     7,    3,    3,    0,    1 },
   ilev_18  =  { 1,   191492.71,     3,    1,    1,    0,    1 },
   ilev_19  =  { 1,   193346.99,     3,    0,    0,    1,    2 },
   ilev_20  =  { 1,   193663.51,     1,    0,    0,    0,    2 },
   ilev_21  =  { 1,   193800.71,     9,    1,    1,    1,    1 },
   ilev_22  =  { 1,   193917.15,    15,    2,    2,    1,    2 },
   ilev_23  =  { 1,   193918.29,     5,    2,    2,    0,    2 },
   ilev_24  =  { 1,   193921.12,    21,    3,    3,    1,    1 },
   ilev_25  =  { 1,   193921.13,     7,    3,    3,    0,    1 },
   ilev_26  =  { 1,   193921.62,    27,   -1,   -1,    1,    2 },
   ilev_27  =  { 1,   193921.62,     9,   -1,   -1,    0,    2 },
   ilev_28  =  { 1,   193942.46,     3,    1,    1,    0,    1 },
   ilev_29  =  { 1,   194936.12,     3,    0,    0,    1,    2 },
   ilev_30  =  { 1,   195114.87,     1,    0,    0,    0,    2 },
   ilev_31  =  { 1,   195192.75,     9,    1,    1,    1,    1 },
   ilev_32  =  { 1,   195260.07,    15,    2,    2,    1,    2 },
   ilev_33  =  { 1,   195260.77,     5,    2,    2,    0,    2 },
   ilev_34  =  { 1,   195262.42,    21,    3,    3,    1,    1 },
   ilev_35  =  { 1,   195262.43,     7,    3,    3,    0,    1 },
   ilev_36  =  { 1,   195262.72,    27,   -1,   -1,    1,    2 },
   ilev_37  =  { 1,   195262.73,     9,   -1,   -1,    0,    2 },
   ilev_38  =  { 1,   195262.79,    33,   -1,   -1,    1,    1 },
   ilev_39  =  { 1,   195262.79,    11,   -1,   -1,    0,    1 },
   ilev_40  =  { 1,   195274.91,     3,    1,    1,    0,    1 },
   ilev_41  =  { 1,   195868.24,     3,    0,    0,    1,    2 },
   ilev_42  =  { 1,   195978.89,     1,    0,    0,    0,    2 },
   ilev_43  =  { 1,   196027.32,     9,    1,    1,    1,    1 },
   ilev_44  =  { 1,   196069.67,    15,    2,    2,    1,    2 },
   ilev_45  =  { 1,   196070.13,     5,    2,    2,    0,    2 },
   ilev_46  =  { 1,   196071.18,    21,    3,    3,    1,    1 },
   ilev_47  =  { 1,   196071.18,     7,    3,    3,    0,    1 },
   ilev_48  =  { 1,   196071.37,    27,   -1,   -1,    1,    2 },
   ilev_49  =  { 1,   196071.37,     9,   -1,   -1,    0,    2 },
   ilev_50  =  { 1,   196071.41,    33,   -1,   -1,    1,    1 },
   ilev_51  =  { 1,   196071.41,    11,   -1,   -1,    0,    1 },
   ilev_52  =  { 1,   196071.43,    39,   -1,   -1,    1,    2 },
   ilev_53  =  { 1,   196071.43,    13,   -1,   -1,    0,    2 },
   ilev_54  =  { 1,   196079.09,     3,    1,    1,    0,    1 },
   ilev_55  =  { 1,   196461.36,     3,    0,    0,    1,    2 },
   ilev_56  =  { 1,   196534.56,     1,    0,    0,    0,    2 },
   ilev_57  =  { 1,   196566.71,     9,    1,    1,    1,    1 },
   ilev_58  =  { 1,   196595.06,    15,    2,    2,    1,    2 },
   ilev_59  =  { 1,   196595.37,     5,    2,    2,    0,    2 },
   ilev_60  =  { 1,   196596.08,    21,    3,    3,    1,    1 },
   ilev_61  =  { 1,   196596.08,     7,    3,    3,    0,    1 },
   ilev_62  =  { 1,   196596.21,    27,   -1,   -1,    1,    2 },
   ilev_63  =  { 1,   196596.21,     9,   -1,   -1,    0,    2 },
   ilev_64  =  { 1,   196596.24,    33,   -1,   -1,    1,    1 },
   ilev_65  =  { 1,   196596.24,    11,   -1,   -1,    0,    1 },
   ilev_66  =  { 1,   196596.25,    39,   -1,   -1,    1,    2 },
   ilev_67  =  { 1,   196596.25,    13,   -1,   -1,    0,    2 },
   ilev_68  =  { 1,   196596.25,    45,   -1,   -1,    1,    1 },
   ilev_69  =  { 1,   196596.25,    15,   -1,   -1,    0,    1 },
   ilev_70  =  { 1,   196601.40,     3,    1,    1,    0,    1 },
   ilev_71  =  { 1,   196861.99,     3,    0,    0,    1,    2 },
   ilev_72  =  { 1,   196912.90,     1,    0,    0,    0,    2 },
   ilev_73  =  { 1,   196935.33,     9,    1,    1,    1,    1 },
   ilev_74  =  { 1,   196955.23,    15,    2,    2,    1,    2 },
   ilev_75  =  { 1,   196955.45,     5,    2,    2,    0,    2 },
   ilev_76  =  { 1,   196955.94,    21,    3,    3,    1,    1 },
   ilev_77  =  { 1,   196955.95,     7,    3,    3,    0,    1 },
   ilev_78  =  { 1,   196956.04,    27,   -1,   -1,    1,    2 },
   ilev_79  =  { 1,   196956.04,     9,   -1,   -1,    0,    2 },
   ilev_80  =  { 1,   196956.06,    33,   -1,   -1,    1,    1 },
   ilev_81  =  { 1,   196956.06,    11,   -1,   -1,    0,    1 },
   ilev_82  =  { 1,   196956.07,    39,   -1,   -1,    1,    2 },
   ilev_83  =  { 1,   196956.07,    13,   -1,   -1,    0,    2 },
   ilev_84  =  { 1,   196956.07,    45,   -1,   -1,    1,    1 },
   ilev_85  =  { 1,   196956.07,    15,   -1,   -1,    0,    1 },
   ilev_86  =  { 1,   196959.69,     3,    1,    1,    0,    1 },
   ilev_87  =  { 1,   197145.23,     3,    0,    0,    1,    2 },
   ilev_88  =  { 1,   197182.06,     1,    0,    0,    0,    2 },
   ilev_89  =  { 1,   197198.33,     9,    1,    1,    1,    1 },
   ilev_90  =  { 1,   197212.82,    15,    2,    2,    1,    2 },
   ilev_91  =  { 1,   197212.99,     5,    2,    2,    0,    2 },
   ilev_92  =  { 1,   197213.35,    21,    3,    3,    1,    1 },
   ilev_93  =  { 1,   197213.35,     7,    3,    3,    0,    1 },
   ilev_94  =  { 1,   197213.42,    27,   -1,   -1,    1,    2 },
   ilev_95  =  { 1,   197213.42,     9,   -1,   -1,    0,    2 },
   ilev_96  =  { 1,   197213.44,    33,   -1,   -1,    1,    1 },
   ilev_97  =  { 1,   197213.44,    11,   -1,   -1,    0,    1 },
   ilev_98  =  { 1,   197213.44,    39,   -1,   -1,    1,    2 },
   ilev_99  =  { 1,   197213.44,    13,   -1,   -1,    0,    2 },
   ilev_100  =  { 1,   197213.44,    45,   -1,   -1,    1,    1 },
   ilev_101  =  { 1,   197213.44,    15,   -1,   -1,    0,    1 },
   ilev_102  =  { 1,   197216.09,     3,    1,    1,    0,    1 },
   -- ===========================================================
}

