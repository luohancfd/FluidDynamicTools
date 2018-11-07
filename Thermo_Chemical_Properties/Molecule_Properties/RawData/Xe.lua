-- Author: Daniel F. Potter
-- Date: 10-Jan-2013

Xe = {}
Xe.M = {
   value = 131.2930000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
Xe.atomic_constituents = {Xe=1}
Xe.charge = 0
Xe.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}

-- Nonequilibrium data

Xe.species_type = "monatomic"
Xe.s_0 = {
   value = 1292.38,
   units = 'J/kg-K',
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
Xe.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
Xe.I = {
   value = 8.9140494e6,
   units = 'J/kg',
   description = 'Ground state ionization energy',
   reference = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html'
}
Xe.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
Xe.electronic_levels = {
   n_levels = 13,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html, 2013-01-10 14:46:17.823710',
   -- This reduced level set was created by grouping individual levels with energies
   -- within 5% of each other
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity
   -- ===========================================================
   ilev_0 = {  5,        0.00,     1,   -1,    0,    0,    2 },
   ilev_1 = {  5,    67434.15,     8,    0,   -1,    1,    1 },
   ilev_2 = {  5,    79173.40,    48,   -1,   -1,   -1,   -1 },
   ilev_3 = {  5,    82935.23,    28,   -1,   -1,   -1,   -1 },
   ilev_4 = {  5,    90746.74,   224,   -1,   -1,   -1,   -1 },
   ilev_5 = {  5,    96203.77,  1407,   -1,   -1,   -1,   -1 },
   ilev_6 = {  5,   102109.35,    79,   -1,   -1,   -1,   -1 },
   ilev_7 = {  5,   106557.37,   156,   -1,   -1,   -1,   -1 },
   ilev_8 = {  5,   168985.00,     3,    1,   -1,   -1,    1 },
   ilev_9 = {  5,   186080.06,    54,   -1,   -1,   -1,   -1 },
   ilev_10 = {  5,   189432.50,    12,   -1,   -1,   -1,   -1 },
   ilev_11 = {  5,   216828.36,    33,   -1,   -1,   -1,   -1 },
   ilev_12 = {  5,   229255.31,    48,   -1,   -1,   -1,   -1 },
   -- ===========================================================
}
