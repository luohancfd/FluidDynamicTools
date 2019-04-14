#!/usr/bin/env python3
# encoding: utf-8

from ReadThermo import ReadThermoFile
from periodic import GetMass
from scipy import constants
import numpy as np

u0 = constants.u
Rgas = constants.R  #J/mol/K

THERMO,_ = ReadThermoFile('./database/therm.yaml')

def S2Rate(formula, stickCoeff=1.0, stoich=1, siteOccupy=1, siteDensity=1.0826715E-9):
    '''
    Calculate preexponential factor of surface reaction rate based on
    Eq4-8 in chemkin theory manual

    The default site density is for SiC epitaxal growth, check doi: 10.1149/1.2085688
    SiC : 6.52E14 molecules/cm^2 = 6.52 molecules/nm^2 doi: 10.1149/1.2085688

    Arguments:
        formula: the formula of specie attached to the surface
        stickCoeff: sticking coefficient
        stoich: the sum of all the stoichiometric coefficients of reactants that are surface species
        siteOccupy: production of [the number of sites that surface species in reactants occupies] rasied to [the correspoding stochimetric]
                    stiteOccupy can be:
                        - a single number, which is the value of these
                        - a list like the following: [ [1,1], [2,3]], which means one surface specie occupies 1 site has a reaction order 1,
                                                     the other one occupies 2 site and has a reaction order 3. YOU NEED TO MAKE SURE 1+3 = stoich
    Note:
        The reaction order of gas phase specie must be one


    Example:
    1.  IF the surf.inp is defined as the following
            SITE /SIC/  SDEN /1.0826715E-9/
                SIH2(S) SIH(S) CH(S)
                SI(S) C(S)
            END
        For the following reactions
            SiH + C(S) => C(b) + SiH(s),  sticking coefficient = 0.94
                S2Rate({"Si":1,"H":1}, stickCoeff = 0.94, siteDensity=1.0826715E-9)
            C2H5 + 2Si(S) => C(s) + CH(s) + 2H2 + 2Si(b), sticking coefficient = 1.0
                S2Rate({"C":2,"H":5}, stickCoeff = 1.0, stochi = 2,siteDensity=1.0826715E-9)

    '''
    if isinstance(formula, float) or isinstance(formula, int):
        mass = formula / 1000
    elif isinstance(formula, dict):
        mass = GetMass(formula) / 1000  # kg/mol
    else:
        mass = THERMO[formula]['m0']/1000

    if isinstance(siteOccupy, list):
        _stoich = sum([i[1] for i in siteOccupy])
        if _stoich != stoich:
            raise(ValueError('Inconsistent siteOccupy'))
        else:
            a = 1.0
            for [i,j] in siteOccupy:
                a *= i ** j
            siteOccupy = a

    rate = np.sqrt(Rgas/2/np.pi/mass) * 100 * siteOccupy / siteDensity ** stoich * stickCoeff  # cc/mol/s
    return rate

if __name__ == "__main__":
    sp = {'CH4': [5e-5, 1],
          'CH3': [0.01, 1],
          'CH2': [0.01, 1],
          'CH': [0.01, 1],
          'C2H6': [0.0016, 2],
          'C2H5': [0.03, 2],
          'C2H4': [0.0016, 2],
          'C2H3': [0.03, 2],
          'C2H2': [0.02, 2],
          'C2H': [0.03, 2],
          'SIH3': [0.5, 1],
          'SIH2': [0.7, 1],
          'SIH': [0.94, 1],
          'SI': [1.0, 1],
          'SI2C': [1.0, 1],
          'SI2': [1.0, 2],
          'SI3': [1.0, 3]}
    for key,val in sp.items():
        r = S2Rate(key, val[0], val[1])
        print('{:4s}   A = {:e}   log10(A) = {:2f}'.format(key,r,np.log10(r)))
