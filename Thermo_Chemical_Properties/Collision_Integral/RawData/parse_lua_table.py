# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 01:20:47 2018

@author: Han
"""
from slpp import slpp as lua
import json
import re


file = 'Wright.txt'
content = []
with open(file, 'r') as f:
    for line in f:
        if re.match('}\s*CI',line):
            content.append('}\n')
            content.append('CI'+line.split('CI')[1])
        else:
            content.append(line)

sps_range = []
isp = -1
for i, line in enumerate(content):
    if 'CI' in line:
        isp += 1
        sps_range.append( [i, 0])
        if isp > 0:
            sps_range[isp-1][1] = i - 1
        content[i] = '{'
    elif 'i =' in line or 'j =' in line:
            w = line.split('=')
            content[i] = ''.join([w[0],' = '+w[1].replace('\'\'','"')])

sps_range[isp][1] = i
        
        
sps = []
for i,j in sps_range:
    w = ''.join(content[i:j+1])
    w = w.replace('\n','')
    w = lua.decode(w)
    if w['parameters']:
        if isinstance(w['parameters'], list):
            if len(w['parameters']) == 1:
                w['parameters'] = w['parameters'][0]
        sps.append(w)
    else:
        print('No param for %s - %s'%(w['i'],w['j']))


output = json.dumps(sps, sort_keys=False, indent=2)
with open(file.replace('txt','json'), 'w', encoding='utf8') as f:
    iline = 0
    w = iter(output.split('\n'))
    while True:
        try:
            line = next(w)
            if 'Omega' in line or '"D"' in line:
                f.write(line)
                while ']' not in line:
                    line = next(w).strip()
                    f.write(line)
                f.write('\n')
            else:
                f.write(line+'\n')
        except StopIteration:
            break 