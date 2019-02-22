#!/usr/bin/env python3
# encoding: utf-8
# Read thermo.dat in Chemkin format
# Written by https://github.com/luohancfd
import warnings
import re
from periodic import GetMass, GetAtomicNumber
from collections import OrderedDict
import os

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

def ReadThermoEntry(lines, H0=None, ref=None, S0=None, CAS=None):
    """
    Read a thermodynamics `lines` for one species in a Chemkin file. Return
    a dict containing all the information
    """

    species = lines[0][0:16].strip()
    refcode = lines[0][18:24].strip()
    formula0 = {}
    for i in range(24,40,5):
        element,count = lines[0][i:i+2].strip(), lines[0][i+2:i+5].strip()
        if element:
            try:
                formula0[element]=int(count)
            except ValueError:
                # Chemkin allows float values for the number of atoms, so try this next.
                formula0[element]=int(float(count))
    formula = OrderedDict(sorted(list(formula0.items()), key=lambda x: GetAtomicNumber(x[0])))
    mass = GetMass(formula)
    phase = lines[0][44]
    if phase.upper() != 'G':
        warnings.warn('Specie %s is in %s phase' % (species, phase))
    # Extract the NASA polynomial coefficients
    # Remember that the high-T polynomial comes first!
    try:
        Tmin = float(lines[0][45:55].strip())
        Tmax = float(lines[0][55:65].strip())
        Tint = float(lines[0][65:75].strip())
        coeffHigh = [float(lines[1][i*15:(i+1)*15]) for i in range(5)]
        coeffHigh += [float(lines[2][i*15:(i+1)*15]) for i in range(2)]
        coeffLow = [float(lines[2][i*15:(i+1)*15]) for i in range(2,5)]
        coeffLow += [float(lines[3][i*15:(i+1)*15]) for i in range(4)]
    except (IndexError, ValueError):
        raise ChemkinError('Error while reading thermo entry for species {0}'.format(species))

    data = {'sp':   species,
            'nT':   2,
            'refcode':  refcode,
            'formula': formula,
            'phase':    phase,
            'm0':     mass,
            'Trange': [[Tmin, Tint], [Tint, Tmax]],
            'Tcommon': Tint,
            'coeff': [coeffLow, coeffHigh]
    }

    if H0:
        data['H0'] = H0
    if ref:
        data['ref'] = ref
    if CAS:
        data['CAS'] = CAS
    if S0:
        data['S0'] = S0

    return data

def ReadThermoBlock(content):
    """
    Read THERMO block of Chemkin

    content: a list of str, it should start with 'THERMO' and end with 'END'

    """
    c = content
    n = len(content)
    i = 0
    regCommentLine = re.compile('^!')
    regEndLine = re.compile('^END')
    # Loop to the start of 'THERMO' block
    while regCommentLine.match(c[i]):
        i += 1
    if not re.match('^THERMO', c[i]):
        raise ChemkinError('This is not THERMO Block')
    else:
        i += 1
        Trange = [float(i) for i in c[i].split()]
    print(Trange)

    # Loop all data block
    i += 1
    nsp = 0
    regRef = re.compile('Source:(.*)$')
    regH0 = re.compile('H0\(298K\)\s*=\s*([+-]?\d+\.?\d*)')
    regS0 = re.compile('S0\(298K\)\s*=\s*([+-]?\d+\.?\d*)')
    regCAS = re.compile('CAS="(.*)"')
    regC = [regRef, regH0, regS0, regCAS]

    data = []
    while not regEndLine.match(c[i]):
        commentBlock = []
        dataBlock = []
        # comment lines
        while regCommentLine.match(c[i]):
            commentBlock.append(c[i])
            i += 1
        # data lines
        while not regCommentLine.match(c[i]):
            dataBlock.append(c[i])
            i += 1
            if regEndLine.match(c[i]):
                break


        # extract information from comment lines
        auxopt = [None for i in regC]
        for j in commentBlock:
            for k, regk in enumerate(regC):
                if not auxopt[k]:
                    w = regk.findall(j)
                    if w:
                        w = w[0].strip()
                        try:
                            auxopt[k] = float(w)
                        except ValueError:
                            auxopt[k] = w

        # process dataBlock
        dat = ReadThermoEntry(dataBlock, ref=auxopt[0], H0=auxopt[1], S0=auxopt[2], CAS=auxopt[3])
        print('Specie: %s' % (dat['sp']))
        nsp += 1
        data.append(dat)

    return data


def ReadThermoData(filepath):
    '''
    Read yaml or json format thermo data

    YAML format (latest):
        (CH2O)3:
            formula: {C: 3, H: 6, O: 3}
        (CH3)2SICH2:
            formula: {C: 3, H: 8, SI: 1}

    Json format
        [
            {
                "sp": "(CH2O)3"
                "formula": {C: 3, H: 6, O: 3}
            },
            {
                "sp": "(CH3)2SICH2"
                "formula": {C: 3, H: 8, SI: 1}
            }
        ]
    '''
    fname, fext = os.path.splitext(filepath)
    with open(filepath, 'r', encoding='utf-8') as f:
        if 'yaml' in fext or 'yml' in fext:
            import oyaml as yaml
            thermo = yaml.load(f)
        elif 'json' in fext:
            import json
            thermo = json.load(f)

    if isinstance(thermo, list):
        # old format
        spDatabase = {i['sp']:i for i in thermo}
        for i in spDatabase.values():
            i.pop('sp', None)
    else:
        spDatabase = thermo

    return spDatabase

if __name__ == "__main__":
    with open('database\\therm.dat', 'r', encoding='utf-8') as f:
        content = f.readlines()
    data = ReadThermoBlock(content)

    iWriteJson = True
    if iWriteJson:
        import json
        with open('database\\therm.json', 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        import sys
        import os
        path2 = '\\'.join(os.getcwd().split('\\')[0:-1] + ['Lewis'])
        sys.path.append(path2)
        from loadthermo import PrettyPrint
        z = PrettyPrint(os.path.join(os.getcwd(), 'database\\therm.json'), True)

    iWriteYAML = True
    if iWriteYAML:
        import oyaml as yaml
        yamlData = {i['sp']: i for i in data}
        for i in yamlData.values():
            i.pop('sp')
            #i['formula'] = dict(i['formula'])
        with open('database\\therm.yaml', 'w', encoding='utf-8') as f:
            yaml.dump(yamlData, f, encoding='utf-8')