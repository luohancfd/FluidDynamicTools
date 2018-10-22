#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:53:32 2018

@author: Han Luo
"""

## read capitelli data
#http://arc.aiaa.org/doi/10.2514/2.6517
import numpy as np
import re
import json

mass = {'N':14.007,'O':15.999,'H':1.008,'C':12.011,'Ar':39.948,
        'N2':14.007*2, 'O2':31.9988, 'NO':30.0061}

with open('DATA_Capitelli.txt', 'r') as f:
    lines = iter(f.readlines())
re1 = re.compile('([eA-Z]{1,2}\d?[+-]?)-([A-Z]{1,2}\d?[+-]?)\s*(.*)')

ref = '''Capitelli, M., C. Gorse, S. Longo, and D. Giordano. "Collision Integrals of High-Temperature Air Species." Journal of Thermophysics and Heat Transfer 14, no. 2 (April 2000): 259–268. doi:10.2514/2.6517.'''
doi = "10.2514/2.6517"


data = []
while True:
    try:
        line = next(lines)
        temp = re1.findall(line)
        sp1 = temp[0][0];
        sp2 = temp[0][1];
        # analysze mass
        sp1 = sp1.replace('+','_plus').replace('-','_minus')
        sp2 = sp2.replace('+','_plus').replace('-','_minus')
        omega = []
        omega.append([float(i) for i in temp[0][2].split()])
        for j in range(3):
            line = next(lines)
            omega.append([float(i) for i in line.split()]);
        data.append({'i':sp1,
                     'j':sp2,
                     'reference':ref,
                     'doi': doi,
                     'model': 'Capitelli et al. Heavy Particle',
                     'mathform': '(a1+a2*T^a3)/(a4+a5*T^a6)',
                     'parameters':{
                             'T_low':50.0,
                             'T_high':100000.0,
                             'Omega11': omega[0],
                             'Omega12': omega[1],
                             'Omega13': omega[2],
                             'Omega22': omega[3]}})
        print('find %s - %s\n'%(sp1,sp2));
    except StopIteration:
        break


# load change-neutral
with open('Heavy-charge.txt', 'r') as f:
    lines = iter(f.readlines())

while True:
    try:
        line = next(lines)
        temp = re1.findall(line)
        sp1 = temp[0][0];
        sp2 = temp[0][1];
        # analysze mass
        sp1 = 'e_minus'
        sp2 = sp2.replace('+','_plus').replace('-','_minus')
        omega = []
        omega.append([float(i) for i in temp[0][2].split()])
        for j in range(5):
            line = next(lines)
            omega.append([float(i) for i in line.split()]);
        data.append({'i':sp1,
                     'j':sp2,
                     'reference':ref,
                     'doi': doi,
                     'model': 'Capitelli et al. Heavy Particle',
                     'mathform': 'check paper',
                     'parameters':{
                             'T_low':50.0,
                             'T_high':100000.0,
                             'Omega11': omega[0],
                             'Omega12': omega[1],
                             'Omega13': omega[2],
                             'Omega14': omega[3],
                             'Omega15': omega[4],
                             'Omega22': omega[5]}})
        print('find %s - %s\n'%(sp1,sp2));
    except StopIteration:
        break

output = json.dumps(data, sort_keys=False, indent=2)
with open('Capitelli.json', 'w', encoding='utf8') as f:
    iline = 0
    w = iter(output.split('\n'))
    while True:
        try:
            line = next(w)
            if 'Omega' in line:
                f.write(line)
                while ']' not in line:
                    line = next(w).strip()
                    f.write(line)
                f.write('\n')
            else:
                f.write(line+'\n')
        except StopIteration:
            break