# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 02:12:38 2018

@author: Han
"""

import glob
import json

sps = glob.glob('RawData/*.lua')
sps = [i for i in sps if 'lua_post.lua' not in i]
files = [sp.replace('RawData\\','') for sp in sps]
sps = [sp.replace('.lua', '') for sp in files]

with open('lua_post.lua', 'w') as f:
    f.write('require("json")\n')    
    f.write('sps = {}\n')
    
    for i, ifile in enumerate(files):
        f.write('dofile("%s")\n'%(ifile))
        f.write('sps.%s = %s\n'%(sps[i], sps[i]))
    
    f.write('file = io.open ("sps.json", "w")\n')
    f.write('io.output(file)\n')
    f.write('io.write(json.encode(sps))\n')
    f.write('io.close(file)\n')
 

# run the lua script lua_post.lua to get json file
with open('sps.json', 'r') as f:
    json_data = json.load(f)

output = json.dumps(json_data, sort_keys=False, indent=2)
with open('sps.json', 'w', encoding='utf8') as f:
    iline = 0
    w = iter(output.split('\n'))
    while True:
        try:
            line = next(w)
            if '"coeffs"' in line or '"ilev' in line:
                f.write(line)
                
                line = next(w)
                while ']' not in line:
                    if ',' in line:
                        f.write(line.strip()+' ')
                    else:
                        f.write(line.strip())
                    line = next(w)
                f.write(line.strip()+'\n')
            else:
                f.write(line+'\n')
        except StopIteration:
            break 
