# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 02:47:19 2018

@author: Han
"""

import json
import re

with open('gupta_us3d.db', 'r') as f:
    content = f.readlines()


reftxt = content[1848:-1]
### parse ref first
refiter = iter(reftxt)
ref = {}
refhead = re.compile('^\((\w+)\)\s*(.*)$')
while True:
    try:
        line = next(refiter).strip()
        w = refhead.findall(line)
        if w:
            refid = w[0][0]
            refnote = w[0][1]
            line = next(refiter).strip()
            while line:
                refnote += '\n'+line
                line = next(refiter).strip()
            while not line:
                line = next(refiter).strip()
                
                
            refpaper = line
            line = next(refiter).strip()
            while line:
                refpaper += ' '+line
                line = next(refiter).strip()
            while not line:
                line = next(refiter).strip()               
            
            refnote2 = ' '
            if not re.match('^---', line):
                refnote2 += line
                while line:
                    refnote2 += '\n'+line
                    line = next(refiter).strip()
                while not line:
                    line = next(refiter).strip()   
                    
            #print([refid, refpaper,refnote])
            
            ref[refid] = {'reference': refpaper, 'refnote':refnote}
            if refid == 'zz':
                break
           # ref[refid] = [refpaper,refnote]

    except StopIteration:
        break 
    


param = content[7:1783]
paramiter = iter(param)
data = []
paramhead = re.compile('^\s*(\w+(?:-|\+)?)-(\w+(?:-|\+)?)\s*(\d)\s*(\d)\s*(\d)\s*\((\w+)\)')
ipair = -1
while True:
    try:
        line = next(paramiter).strip()
        w = paramhead.findall(line)
        if line:
            ipair += 1
            sp1 = w[0][0]
            sp2 = w[0][1]
            data.append({'i':sp1,
                         'j':sp2,
                         'model': 'GuptaYos2',
                         })
            if w[0][5] in ref.keys():
                thisref = ref[w[0][5]]['reference']
            else:
                thisref = ''
            data[ipair]['reference'] = thisref
            
            if w[0][2] == '1':
                line = next(paramiter).strip()
                piomega11 = [float(i) for i in line.split()]
                data[ipair]['Pi_Omega_11'] = piomega11
            
            if w[0][3] == '1':
                line = next(paramiter).strip()
                piomega22 = [float(i) for i in line.split()]
                data[ipair]['Pi_Omega_22'] = piomega22
    
                
            if w[0][4] == '1':
                line = next(paramiter).strip()
                bij = [float(i) for i in line.split()]
                data[ipair]['Bij'] = bij
    except StopIteration:
        break 
    
    
output = json.dumps(data, sort_keys=False, indent=2)
with open('Gupta_us3d.json', 'w', encoding='utf8') as f:
    iline = 0
    w = iter(output.split('\n'))
    while True:
        try:
            line = next(w)
            if 'Omega' in line or '"D"' in line or '"Bij"' in line:
                f.write(line)
                while ']' not in line:
                    line = next(w).strip()
                    f.write(line)
                f.write('\n')
            else:
                f.write(line+'\n')
        except StopIteration:
            break    
        