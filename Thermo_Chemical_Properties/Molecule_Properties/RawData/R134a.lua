-- Author: Peter Blyton
-- Date: 18-Apr-2012
-- also known as 1,1,1,2-Tetrafluoroethane CH2FCF3

R134a = {}
R134a.M = {
   value = 102.030e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'McLinden, MO et al (1989). Measurement and Formulation of the Thermodynamic Properties of Refrigerants 134a and 123',
}
R134a.atomic_constituents = {C=2,H=2,F=4}
R134a.charge = 0
R134a.T_c = {
   value = 374.179,
   units = 'K',
   description = 'critical temperature',
   reference = 'Huber, ML and MO McLinden (1992). Thermodynamic Properties of R134a',
}
R134a.p_c = {
   value = 4056.0e+3,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Huber, ML and MO McLinden (1992). Thermodynamic Properties of R134a',
}
R134a.rho_c = {
   value = 513.293, -- 5.0308 mol/L
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'Huber, ML and MO McLinden (1992). Thermodynamic Properties of R134a',
   note = 'To be used with the MBWR Coefficients of Huber and McLinden (1992)',
}
R134a.rho_c_disabled = {
   value = 515.3,
   units = 'kg/m^3',
   description = 'critical density',
   reference = 'DuPont (2004). Thermodynamic Properties of HFC-134a',
   note = 'To be used with the MBWR Coefficients of DuPont (2004)',
}
R134a.reference_state = {
   description = 'Reference temperature, internal energy and entropy',
   reference = 'None',
   note = 'IIR Ideal gas reference where H=200 kJ/kg, S=1 kJ/kg.K for saturated liquid at T0',
   units = 'K, J/kg, J/kg.K',
   T0 = 273.15,
   u0 = 3.83916e5,
   s0 = 1.9575e3,
}
R134a.Cp0_coeffs = {
   description = 'Coefficients for Cp0 polynomial.',
   reference = 'McLinden, MO et al (1989). Measurement and Formulation of the Thermodynamic Properties of Refrigerants 134a and 123',
   units = 'J/kg.K', -- Resulting unit when using these coefficients in Cp0 polynomial.
   T_low  = 169.85,
   T_high = 450.0,
   G = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
        0.0,
        1.90146e+02, -- McLinden published 1.94006e+01, 2.58531e-01, -1.29665e-04 in J/mol.K
        2.53387e+00,
       -1.27085e-03,
        0.0,
        0.0}
}
R134a.MBWR_EOS_coeffs = {
   description = 'Coefficients for MBWR equation of state',
   reference = 'Huber, ML and MO McLinden (1992). Thermodynamic Properties of R134a',
   units = 'bar', -- Resulting unit when using these coefficients in p-rho-T equation.
   b = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
        9.6520936222e-02,  -4.0182476889e+00,   3.9523953286e+01,
        1.3453286896e+03,  -1.3943974135e+06,  -3.0928135518e-03,
        2.9238151228e+00,  -1.6514661356e+03,   1.5070600312e+06,
        5.3497394831e-05,   5.4393331762e-01,  -2.1132604976e+02,
       -2.6819120385e-02,  -5.4106712595e-01,  -8.5173177940e+02,
        2.0518825365e-01,  -7.3305018809e-03,   3.8065596386e+00,
       -1.0583208759e-01,  -6.7924308442e+05,  -1.2699837860e+08,
       -4.2623443183e+04,   1.0197333823e+09,  -1.8669952678e+02,
       -9.3342632342e+04,  -5.7173520896e+00,  -1.7676273879e+05,
       -3.9728275231e-02,   1.4301684480e+01,   8.0308529426e-05,
       -1.7195907355e-01,   2.2623838566e+00}
}
R134a.MBWR_EOS_coeffs_disabled = {
   description = 'Coefficients for MBWR equation of state',
   reference = 'DuPont (2004). Thermodynamic Properties of HFC-134a',
   units = 'bar', -- Resulting unit when using these coefficients in p-rho-T equation.
   b = {0.0, -- First element zero to align array indexes with equation coefficient numbers.
       -6.5455235227e-02,   5.8893751817e+00,  -1.3761788409e+02,
        2.2693168845e+04,  -2.9262613296e+06,  -1.1923776190e-04,
       -2.7214194543e+00,   1.6295253680e+03,   7.2942203182e+05,
       -1.1724519115e-04,   8.6864510013e-01,  -3.0660168246e+02,
       -2.5664047742e-02,  -2.4381835971e+00,  -3.1603163961e+02,
        3.4321651521e-01,  -1.0154368796e-02,   1.1734233787e+00,
       -2.7301766113e-02,  -6.6338502898e+05,  -6.4754799101e+07,
       -3.7295219382e+04,   1.2614735899e+09,  -6.4742200070e+02,
        1.2362450399e+05,  -1.5699196293e+00,  -5.1848932204e+05,
       -8.1396321392e-02,   3.0325168842e+01,   1.3399042297e-04,
       -1.5856192849e-01,   9.0679583743e+00}
}
