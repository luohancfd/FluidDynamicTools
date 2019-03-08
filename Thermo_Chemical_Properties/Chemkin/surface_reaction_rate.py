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

    Only support:
      A + nsite*S => ....
    where S is the site

    Arguments:
        formula: eg: {'H':1}
    '''
    mass = GetMass(formula) / 1000  #kg/mol
    rate = np.sqrt(Rgas/2/np.pi/mass) * 100 * (siteOccupy / siteDensity)**nsite * stickCoeff  # cc/mol/s
    return rate

if __name__ == "__main__":
    print('{:e}'.format(S2Rate({'H':1})))
    print('{:e}'.format(S2Rate({'Si':1})))