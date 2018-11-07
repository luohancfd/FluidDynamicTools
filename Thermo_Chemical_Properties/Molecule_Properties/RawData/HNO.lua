-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

HNO = {}
HNO.M = {
   value = 31.01404e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HNO.atomic_constituents = {H=1,N=1,O=1}
HNO.charge = 0
HNO.gamma = {
   value = 1.325,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HNO.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -6.854764860e+04,  9.551627200e+02, -6.000720210e-01, 
                7.995176750e-03, -6.547079160e-07, -3.670513400e-09,
                1.783392519e-12,  6.435351260e+03,  3.048166179e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -5.795614980e+06,  1.945457427e+04, -2.152568374e+01,
                1.797428992e-02, -4.976040670e-06,  6.397924170e-10,
               -3.142619368e-14, -1.104192372e+05,  1.818650338e+02
    }
  },
  ref="from CEA2::thermo.inp"
}
