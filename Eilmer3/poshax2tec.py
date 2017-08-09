#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 19:40:14 2017

@author: Han Luo

Convert poshax3 result to tecplot format
"""
import sys
import numpy as np


uinf = float(sys.argv[1])
column = list()
data = list()
with open('output.data') as f:
    line = f.readline();
    while line:
        if 'Columns' in line:
            line = f.readline()
            while '#' in line:
                column.append(line.split(':')[1][:-1])
                line = f.readline()
            while line:
                data.append(line)
                line = f.readline()
        line = f.readline()

m = len(data)
time = [None,]*m
for i,v in enumerate(data):
    time[i] = float(v.split()[0])/uinf*1.e6
#n = data.shape[1]+1

with open('Result.dat','w') as f:
    f.write('TITLE = "Shock simulation u = %8.5f km/s"\n'%(uinf/1000,))
    f.write('VARIABLES = ')
    f.write('"time (<greek>m</greek>s)" ')
    for i in column:
        f.write('"'+i+'" ')
    f.write('\n')
    f.write('ZONE I=%d T="u=%8.5f km/s"\n'%(m,uinf/1000,))
    for i in range(m):
        f.write(str(time[i])+'   ')
        f.write(data[i])







