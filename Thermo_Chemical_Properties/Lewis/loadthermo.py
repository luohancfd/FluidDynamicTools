# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 18:35:18 2018

@author: Han
"""

import json

def PrettyPrint(file='thermo.json', print=False):
    with open(file,'r') as f:
        data = json.load(f)

    if print:
        output = json.dumps(data, sort_keys=False, indent=2)
        with open(file, 'w', encoding='utf8') as f:
            iline = 0
            w = iter(output.split('\n'))
            while True:
                try:
                    line = next(w)
                    if '"Trange"' in line or '"Hdiff"' in line or '"coeff"' in line:
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
                                    f.write(line.strip() + '\n')
                                    break
                        else:
                            f.write('\n')
                    elif '"formula"' in line:
                        f.write(line)
                        while True:
                            line = next(w)
                            if '}' in line:
                                f.write(line.strip() + '\n')
                                break
                            elif ',' in line:
                                f.write(line.strip() + ' ')
                            else:
                                f.write(line.strip())
                    else:
                        f.write(line+'\n')
                except StopIteration:
                    break
    return data