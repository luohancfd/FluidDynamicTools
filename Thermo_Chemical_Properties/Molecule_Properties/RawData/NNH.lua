-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

NNH = {}
NNH.M = {
   value = 29.021340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
NNH.atomic_constituents = {N=2,H=1}
NNH.charge = 0
NNH.gamma = {
   value = 1.3152e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
NNH.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.344692700e00, -4.849707200e-03, 2.005945900e-05, -2.172646400e-08, 7.946953900e-12, 2.879197300e04, 2.977941000e00, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.766754400e00, 2.891508200e-03, -1.041662000e-06, 1.684259400e-10, -1.009189600e-14, 2.865069700e04, 4.470506700e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
