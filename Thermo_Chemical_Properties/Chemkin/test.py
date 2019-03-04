#!/usr/bin/env python3
# encoding: utf-8
import numpy as np
from ReadThermo import ReadThermoFile
from ReadChem import EquilibriumConstant
from scipy.constants import N_A
import matplotlib.pyplot as plt

thermo,_ = ReadThermoFile('database\\therm.yaml')

T = np.arange(300, 2500, 100)
Kp, Kc = EquilibriumConstant('SIH3+SIH4=SI2H6+H', T, thermo)
fr = 8.67E-24 * T**3.0 * np.exp(-8799 / T)  # doi 10.1021/jp1082196  cm^3 /molecule /s
# reaction order don't change, no need to change unit
br = fr/Kc

br2 = 2.4e-10 * np.exp(-1250 / T)

plt.semilogy()


