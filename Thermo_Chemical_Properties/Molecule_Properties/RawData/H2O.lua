-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

H2O = {}
H2O.M = {
   value = 18.01528e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
H2O.atomic_constituents = {H=2,O=1}
H2O.charge = 0
H2O.gamma = {
   value = 1.329,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
H2O.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=373.2, T_high=1073.2, A=0.50019557e+00, B=-0.69712796e+03, C=0.88163892e+05, D=0.30836508e+01},
      {T_low=1073.2, T_high=5000.0, A=0.58988538e+00, B=-0.53769814e+03, C=0.54263513e+05, D=0.23386375e+01},
      {T_low=5000.0, T_high=15000.0, A=0.64330087e+00, B=-0.95668913e+02, C=-0.37742283e+06, D=0.18125190e+01},
      ref = 'from CEA2::trans.inp which cites Sengers and Watson (1986)'
   }
}
H2O.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=373.2, T_high=1073.2, A=0.10966389e+01, B=-0.55513429e+03, C=0.10623408e+06, D=-0.24664550e+00},
      {T_low=1073.2, T_high=5000.0, A=0.39367933e+00, B=-0.22524226e+04, C=0.61217458e+06, D=0.58011317e+01},
      {T_low=5000.0, T_high=15000.0, A=-0.41858737e+00, B=-0.14096649e+05, C=0.19179190e+08, D=0.14345613e+02},
      ref = 'from CEA2::trans.inp which cites Sengers and Watson (1986)'
   }
}
H2O.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -3.947960830e+04,  5.755731020e+02,  9.317826530e-01,
	        7.222712860e-03, -7.342557370e-06,  4.955043490e-09,
	       -1.336933246e-12, -3.303974310e+04,  1.724205775e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.034972096e+06, -2.412698562e+03,  4.646110780e+00,
                2.291998307e-03, -6.836830480e-07,  9.426468930e-11,
               -4.822380530e-15, -1.384286509e+04, -7.978148510e+00
    }
  },
  ref="from CEA2::thermo.inp"
}
H2O.T_c = {
   value = 647.14,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
H2O.p_c = {
   value = 220.64e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}

-- Nonequilibrium data

H2O.species_type = "nonlinear polar polyatomic"
H2O.s_0 = {
   value = 0.0,
   units = 'J/kg-K',
   description = 'Dummy standard state entropy at 1 bar',
   reference = 'NA'
}
H2O.h_f = {
   value = -13423382.82,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
H2O.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
H2O.I = {
   value = 0.0,
   units = 'J/kg',
   description = 'Dummy ground state ionization energy',
   reference = 'none'
}
H2O.electronic_levels = {
   n_levels = 1,
   ref = 'JANAF Tables 1985 (3rd edition)',
   -- ====================================================================================================================================
   --   n       Te         re       g       dzero      A0        B0          C0        sigma   sigma_rot   we[0]      we[1]      we[2]
   -- ====================================================================================================================================
   ilev_0  = {      0.00,  0.9587,  1,      0.0,       27.8847,  14.5118,    9.2806,   2,      0,          3651.1,    1594.7,    3755.9 }
   -- ====================================================================================================================================
}

-- Real gas data

H2O.reference_state = {
   description = 'Reference temperature, internal energy and entropy',
   reference = 'None',
   note = 'Ideal gas reference when U=0, S=0 for saturated liquid at T0',
   units = 'K, J/kg, J/kg.K',
   T0 = 273.16,
   u0 = 2.3741e6,
   s0 = 6.6945e3,
}
H2O.Cp0_coeffs = {
   description = 'Coefficients for Cp0 polynomial of form T^(-2) to T^4.',
   reference = 'Polt, A and G Maurer (1992). Fluid Phase Equilibria, 73: 27-38.',
   note = 'Polt published these coefficients in kJ/(kg.K)',
   units = 'J/kg.K', -- Resulting unit when using these coefficients in Cp0 equation.
   T_low  = 273.15, -- Range of validity.
   T_high = 923.91,
   G = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
        0.0,
        1.94077e+03,
       -9.67660e-01,
        3.08550e-03,
       -2.63140e-06,
        8.60950e-10}
}
H2O.Bender_EOS_coeffs = {
   description = 'Coefficients for Bender equation of state',
   reference = 'Polt, A and G Maurer (1992). Fluid Phase Equilibria, 73: 27-38.',
   units = 'Pa', -- Resulting unit when usng these coefficients in p-rho-T equation.
   A = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
       -1.3002084070e-07,   7.6533762730e-04,  -9.0745855500e-01,
       -2.9984514510e+02,  -1.6148196450e+04,   4.4890768230e-06,
       -6.3451282300e-03,   3.6353975160e+00,  -8.1460085550e-06,
        6.0971223530e-03,   1.0426877550e-05,  -1.0516986730e-02,
        2.7385365060e-03,  -2.2453231530e+03,   2.1342907950e+06,
       -3.1262434070e+08,   4.9365926510e+03,  -5.7492854100e+06,
        1.1440257530e+09,
        4.0e-6} -- Final element is "gamma", roughly 1/(rho_c^2) as used in the MBWR EOS.
}
