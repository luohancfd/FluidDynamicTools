-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

HNO2 = {}
HNO2.M = {
   value = 47.01344e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HNO2.atomic_constituents = {H=1,N=1,O=2}
HNO2.charge = 0
HNO2.gamma = {
   value = 1.218,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HNO2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  8.591985060e+03,  1.203644046e+02,  9.412979120e-01,
                1.942891839e-02, -2.253174194e-05,  1.384587594e-08,
               -3.473550460e-12, -1.106337202e+04,  2.073967331e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  8.787904130e+05, -3.990455030e+03,  1.187349269e+01,
               -4.881900610e-04,  7.133636790e-08, -5.376303340e-12,
                1.581778986e-16,  1.246343241e+04, -4.608874688e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
