#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 16:32:07 2017

@author: luo160
"""
import os
rootfolder = next(os.walk('.'))[1]
T = [float(i.split('=')[1]) for i in rootfolder]
with open(os.path.split(os.getcwd())[1]+'_Rate.dat','w') as f:
    f.write('VARIABLES = "Vibrational Temperature (K)", "Actual Vibrational Temperature (K)",'+
            '"Rate (cm<sup>2</sup>mol<sup>-1</su>)", "Zv"\n')
    for i,tt in enumerate(T):
        Tdir = rootfolder[i]
        Tvdir = next(os.walk(Tdir))[1]
        Tv = [float(i.split('=')[1]) for i in Tvdir]
        # sort by temperature
        index = [b[0] for b in sorted(enumerate(Tv),key=lambda i:i[1])]
        Tv = [Tv[i] for i in index]
        Tvdir = [Tvdir[i] for i in index]
        ieq = Tv.index(tt)
        with open(os.path.join('.',Tdir,Tvdir[ieq],'DS1REAC.DAT'),'r') as f1:
            last_line = f1.readlines()[-1]
            eqrate = float(last_line.split()[-2])      

        f.write('ZONE I = %d, T = "%s"\n'%(len(Tvdir),Tdir))
        for j,ctvdir in enumerate(Tvdir):
            with open(os.path.join('.',Tdir,ctvdir,'DS1REAC.DAT'),'r') as f1:
                last_line = f1.readlines()[-1]
            rate = float(last_line.split()[-2])
            ATv = float(last_line.split()[3])
            f.write('%9.2f   %9.2f   %e   %e\n'%(Tv[j],ATv,rate,rate/eqrate))
        f.write('\n')
        