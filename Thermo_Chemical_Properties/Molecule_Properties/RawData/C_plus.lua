-- Collater: Daniel F. Potter
-- Date: 01-Feb-2010

C_plus = {}
C_plus.M = { 
   value = 12.0101514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
C_plus.atomic_constituents = {C=1}
C_plus.charge = 1
C_plus.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
C_plus.CEA_coeffs = {
  { T_low  = 298.150,
    T_high = 1000.0,
    coeffs = {  2.258535929e+03, -1.574575687e+00,  2.503637730e+00,
    	       -5.202878370e-06,  4.516908390e-09, -2.181431053e-12,
    	        4.495047033e-16,  2.168951913e+05,  4.345699505e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.255112551e+04, -3.411874670e+01,  2.543383218e+00,
    	       -2.805120849e-05,  9.751641970e-09, -1.736855394e-12,
    	        1.246191931e-16,  2.171001786e+05,  4.063913515e+00
    }
  },
  { T_low  = 6000.0,
    T_high = 20000.0,
    coeffs = {  5.618135320e+05, -6.047058900e+03,  5.884541470e+00,
    	       -7.211894530e-04,  6.823484110e-08, -2.599878590e-12,
    	        3.633868358e-17,  2.581370458e+05, -2.280019759e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
C_plus.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.80124735e+00, B=0.17261643e+03, C=-0.69940019e+05, D=0.88364870e-01},
      {T_low=5000.0, T_high=15000.0, A=0.10344416e+01, B=0.31310924e+04, C=-0.45512020e+07, D=-0.23102402e+01},
      ref = 'from CEA2::trans.inp which cites Biolsi (1982)'
   }
}
C_plus.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.80224051e+00, B=0.17739617e+03, C=-0.72350849e+05, D=0.10329911e+01},
      {T_low=5000.0, T_high=15000.0, A=0.10355137e+01, B=0.31489830e+04, C=-0.45854028e+07, D=-0.13676372e+01},
      ref = 'from CEA2::trans.inp which cites Biolsi (1982)'
   }
}

-- Nonequilibrium data

C_plus.species_type = "monatomic"
C_plus.s_0 = {
   value = 12877.44,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
C_plus.h_f = {
   value = 150659589.69,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C_plus.I = {
   value = 195887126.33,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html'
}
C_plus.Z = {
   value = 1,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C_plus.electronic_levels = {
   -- n_levels = 42,
   n_levels = 3,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   comments = 'All the individual NIST levels expressed as multiplets',
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,       42.28,     6,    1,    1,    1,    1 },
   ilev_1   =  { 2,    43035.78,    12,   -1,    1,    2,    2 },
   ilev_2   =  { 2,    74931.11,    10,   -1,    2,    1,    2 },
   ilev_3   =  { 2,    96493.74,     2,   -1,    0,    1,    2 },
   ilev_4   =  { 2,   110651.76,     6,   -1,    1,    1,    2 },
   ilev_5   =  { 2,   116537.65,     2,    0,    0,    1,    2 },
   ilev_6   =  { 2,   131731.80,     6,    1,    1,    1,    1 },
   ilev_7   =  { 2,   142027.10,     4,   -1,    0,    2,    1 },
   ilev_8   =  { 2,   145550.13,    10,    2,    2,    1,    2 },
   ilev_9   =  { 2,   150463.62,    10,   -1,    2,    1,    1 },
   ilev_10  =  { 2,   157234.07,     2,    0,    0,    1,    2 },
   ilev_11  =  { 2,   162522.34,     6,    1,    1,    1,    1 },
   ilev_12  =  { 2,   167009.29,    12,    0,    1,    2,    1 },
   ilev_13  =  { 2,   168124.17,    10,    2,    2,    1,    2 },
   ilev_14  =  { 2,   168742.04,     6,   -1,    1,    1,    1 },
   ilev_15  =  { 2,   168978.34,    14,    3,    3,    1,    1 },
   ilev_16  =  { 2,   173347.84,     2,    0,    0,    1,    2 },
   ilev_17  =  { 2,   175292.30,     6,    1,    1,    1,    1 },
   ilev_18  =  { 2,   177787.22,     6,    0,    1,    1,    1 },
   ilev_19  =  { 2,   178495.47,    10,    2,    2,    1,    2 },
   ilev_20  =  { 2,   178955.94,    14,    3,    3,    1,    1 },
   ilev_21  =  { 2,   179073.05,    18,   -1,   -1,    1,    2 },
   ilev_22  =  { 2,   181264.24,     2,    0,    0,    1,    2 },
   ilev_23  =  { 2,   181741.65,    20,    1,    2,    2,    2 },
   ilev_24  =  { 2,   182036.89,     6,    1,    1,    1,    2 },
   ilev_25  =  { 2,   182993.52,     6,    1,    1,    1,    1 },
   ilev_26  =  { 2,   184075.00,    10,    2,    2,    1,    2 },
   ilev_27  =  { 2,   184376.06,    14,    3,    3,    1,    1 },
   ilev_28  =  { 2,   184449.27,    18,   -1,   -1,    1,    2 },
   ilev_29  =  { 2,   184466.50,    22,   -1,   -1,    1,    1 },
   ilev_30  =  { 2,   184690.98,     4,    1,    0,    2,    2 },
   ilev_31  =  { 2,   185732.93,     2,    0,    0,    1,    2 },
   ilev_32  =  { 2,   186452.13,    12,    1,    1,    2,    2 },
   ilev_33  =  { 2,   186746.17,     6,    1,    1,    1,    1 },
   ilev_34  =  { 2,   187353.00,    10,    2,    2,    1,    2 },
   ilev_35  =  { 2,   187641.60,    14,    3,    3,    1,    1 },
   ilev_36  =  { 2,   187691.40,    18,   -1,   -1,    1,    2 },
   ilev_37  =  { 2,   187701.00,    22,   -1,   -1,    1,    1 },
   ilev_38  =  { 2,   188601.54,    10,    1,    2,    1,    2 },
   ilev_39  =  { 2,   189794.20,    18,   -1,   -1,    1,    2 },
   ilev_40  =  { 2,   195786.71,    28,    2,    3,    2,    1 },
   ilev_41  =  { 2,   196572.80,    20,    2,    2,    2,    1 },
   -- ===========================================================
}
C_plus.spradian_electronic_levels = {
   n_levels = 12,
   ref = 'Spradian07::atom.dat',
   -- ===================================
   --  Level   n      E(cm-1)        g
   -- ===================================
   ilev_0  = { 2,        0.0,        2 },
   ilev_1  = { 2,       64.0,        4 },
   ilev_2  = { 2,    43033.0,       12 },
   ilev_3  = { 2,    74932.0,       10 },
   ilev_4  = { 2,    96494.0,        2 },
   ilev_5  = { 2,   110652.0,        6 },
   ilev_6  = { 3,   116538.0,        2 },
   ilev_7  = { 3,   131732.0,        6 },
   ilev_8  = { 2,   142024.0,        4 },
   ilev_9  = { 3,   145551.0,       10 },
   ilev_10 = { 2,   150460.0,       10 },
   ilev_11 = { 4,   157234.4,        2 },
   -- ===================================
}
