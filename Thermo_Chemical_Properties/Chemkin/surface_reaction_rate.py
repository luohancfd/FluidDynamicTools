#!/usr/bin/env python3
# encoding: utf-8

from ReadThermo import ReadThermoFile
from periodic import GetMass
from scipy import constants
import numpy as np

u0 = constants.u
Rgas = constants.R  #J/mol/K

THERMO,_ = ReadThermoFile('./database/therm.yaml')

def S2Rate(formula, stickCoeff=1.0, nsite=1, siteOccupy=1, siteDensity=1.0826715E-9):
    '''
    Calculate preexponential factor of surface reaction rate based on
    Eq4-8 in chemkin theory manual

    The default site density is for SiC epitaxal growth, check doi: 10.1149/1.2085688
    SiC : 6.52E14 molecules/cm^2 = 6.52 molecules/nm^2 doi: 10.1149/1.2085688



    Only support:
      A + nsite*S => ....
    where S is the site

    Arguments:
        formula: eg: {'H':1}
    '''
    if isinstance(formula, dict):
        mass = GetMass(formula) / 1000  #kg/mol
    else:
        mass = THERMO[formula]['m0']/1000

    rate = np.sqrt(Rgas/2/np.pi/mass) * 100 * (siteOccupy / siteDensity)**nsite * stickCoeff  # cc/mol/s
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
