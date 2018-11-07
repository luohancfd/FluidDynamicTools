-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

HNO3 = {}
HNO3.M = {
   value = 63.01284e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HNO3.atomic_constituents = {H=1,N=1,O=3}
HNO3.charge = 0
HNO3.gamma = {
   value = 1.181,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HNO3.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  9.202869010e+03,  1.093774496e+02, -4.521042450e-01,
                2.984914503e-02, -3.190635500e-05,  1.720931528e-08, 
               -3.782649830e-12, -1.764048507e+04,  2.746644879e+01 
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = { -9.497809640e+04, -2.733024468e+03,  1.449426995e+01,
               -7.821868050e-04,  1.702693665e-07, -1.930543961e-11,
                8.870455120e-16, -4.882517780e+03, -5.928392985e+01 
    }
  },
  ref="from CEA2::thermo.inp"
}
