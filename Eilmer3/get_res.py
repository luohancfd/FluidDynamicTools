#!/usr/bin/env python
# created by: Han Luo
# usage: get_red.py A.log
# generate a pdf file A.log.pdf and A.log.DAT
import sys
FILE_IN = sys.argv[1]
FILE_OUT = FILE_IN+'.DAT'
cmd = '''awk 'BEGIN {print "step\\ttime\\tmass_res\\tenergy_res"}
     /mass global/ {printf $7"\\t"$9"\\t"$5"\\t"}
     /energy global/ {printf $5"\\n"}' '''+FILE_IN
import subprocess,shlex
cmd = shlex.split(cmd)
outp = open(FILE_OUT,'w')
p = subprocess.Popen(cmd,stdout=outp)
p.wait()
outp.flush()

import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt(FILE_OUT,skiprows=1)
x = data[:,1]*1.0e9
plt.plot(x,data[:,2],'k',label='mass_res')
plt.plot(x,data[:,3],'b',label='energy_res')
plt.ylabel('Residue')
plt.xlabel('Time [ns]')
if len(sys.argv) > 2:
   plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig(FILE_IN+'.pdf',format='pdf')
subprocess.call(['evince',FILE_IN+'.pdf'])
