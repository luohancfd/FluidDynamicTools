# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 18:35:18 2018

@author: Han
"""

import json
with open('thermo.json','r') as f:
    data = json.load(f)

output = json.dumps(data, sort_keys=False, indent=2)
with open('thermo.json', 'w', encoding='utf8') as f:
    iline = 0
    w = iter(output.split('\n'))
    while True:
        try:
            line = next(w)
            if '"Trange"' in line or '"Hdiff"' in line or '"coeff"' in line or '"coeffb"' in line:
                f.write(line)
                if '[' in line:
                    level = 1
                    while True:
                        line = next(w)
                        if '[' in line:
                            level += 1
                        if ']' in line:
                            level -= 1
                        if level > 0:
                            f.write(line.strip())
                        else:
                            f.write(line.strip()+'\n')
                            break       
                else:
                    f.write('\n')
            else:
                f.write(line+'\n')
        except StopIteration:
            break  