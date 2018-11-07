-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

HO2 = {}
HO2.M = {
   value = 33.00674e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HO2.atomic_constituents = {H=1,O=2}
HO2.charge = 0
HO2.gamma = {
   value = 1.312,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HO2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -7.598882540e+04,  1.329383918e+03, -4.677388240e+00,
                2.508308202e-02, -3.006551588e-05,  1.895600056e-08,
               -4.828567390e-12, -5.873350960e+03,  5.193602140e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -1.810669724e+06,  4.963192030e+03, -1.039498992e+00,
                4.560148530e-03, -1.061859447e-06,  1.144567878e-10,
               -4.763064160e-15, -3.200817190e+04,  4.066850920e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
