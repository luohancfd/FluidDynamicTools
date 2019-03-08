#!/usr/bin/env python3
# encoding: utf-8
# Read thermo.dat in Chemkin format
# Written by https://github.com/luohancfd
# pylint: disable
import logging
import re
import os
import argparse
import numpy as np
from collections import OrderedDict
from periodic import GetMass, GetAtomicNumber
from scipy import constants

logging.basicConfig(level = logging.ERROR,
                    format = '%(asctime)s  %(levelname)-10s %(processName)s  %(name)s %(message)s')

def Chemkin7Par(sp, Tin, thermo, unit='si'):
    '''
    Evaluate Cp, H, S, G for given temperature and given specie

    Arguments:
        sp: name of specie
        T:  temperature list
        thermo: thermo database
                       | Cp          |  S         |  H        |  G
        unit:   si/J:  | J/mol/K     |  J/mol/K   |  kJ/mol   |  kJ/mol
                cal:   | cal/mol/K   |  cal/mol/K |  kcal/mol |  kcal/mol

    '''
    if isinstance(Tin, list):
        T = np.array(Tin, dtype=np.float64)
    elif isinstance(Tin, float) or isinstance(Tin, int):
        T = np.array([Tin], dtype=np.float64)
    elif not isinstance(Tin, np.ndarray):
        raise ValueError('Wrong dtype for T')
    else:
        T = Tin

    if unit == 'si' or unit == 'cal':
        Rgas = constants.R   # J/mol/K
    elif unit == 'cal':
        Rgas = constants.R / constants.calorie
    else:
        raise ValueError('Unknown unit')

    if sp in thermo.keys():
        data = thermo[sp]
    elif sp.upper() in thermo.keys():
        data = thermo[sp.upper()]
    elif sp.capitalize() in thermo.keys():
        data = thermo[sp.capitalize()]
    else:
        raise ValueError('{} is not in thermo database'.format(sp))

    Trange = data['Trange']
    coeff = data['coeff']
    if len(Trange) > 1:
        Tbin = [i[0] for i in Trange]
        Tbin.append(Trange[-1][1])
    else:
        Tbin = Trange[0]

    binIndex = np.digitize(T, Tbin, right=True)
    tIndex = [np.where(binIndex == i)[0].tolist() for i in range(len(Tbin)+1)]

    Cp = np.zeros(T.shape, dtype=np.float64)
    S = np.zeros(T.shape, dtype=np.float64)
    G = np.zeros(T.shape, dtype=np.float64)
    H = np.zeros(T.shape, dtype=np.float64)

    for i, index in enumerate(tIndex):
        gT = T[index]
        if i == 0:
            a = coeff[0]
        elif i == len(tIndex)-1:
            a = coeff[-1]
        else:
            a = coeff[i-1]
        if len(gT) > 0:
            gTons = np.ones(gT.shape, np.float64)
            gTlog = np.log(gT)
            gTm1 = 1.0 / gT
            gT2 = gT * gT
            gT3 = gT2 * gT
            gT4 = gT3 * gT

            pT = np.vstack((gTons, gT, gT2, gT3, gT4))

            nCoeff = np.array(a, dtype=np.float64)

            gCp = Rgas * np.dot(nCoeff[:5], pT)
            gH  = Rgas * gT * (np.dot(nCoeff[:5]/np.arange(1,6), pT) + nCoeff[-2]*gTm1)

            pT = np.vstack((gTlog, gT, gT2, gT3, gT4))
            gS = Rgas * (np.dot(nCoeff[:5]/np.array([1,1,2,3,4]), pT) + nCoeff[-1])

            gG = gH - gT * gS

            gH /= 1000
            gG /= 1000

            Cp[index] = gCp
            S[index] = gS
            G[index] = gG
            H[index] = gH

    return Cp, H, G, S


def ReadThermoEntry(linesIn, H0=None, ref=None, S0=None, CAS=None, name=None):
    """
    Read a thermodynamics `lines` for one species in a Chemkin file. Return
    a dict containing all the information
    """
    lines = [i for i in linesIn if not re.match('^\s*!', i)]
    species = lines[0][0:16].strip()
    refcode = lines[0][18:24].strip()
    formula0 = {}
    for i in range(24, 40, 5):
        element,count = lines[0][i:i+2].strip(), lines[0][i+2:i+5].strip()
        if element:
            try:
                formula0[element] = int(count)
            except ValueError:
                # Chemkin allows float values for the number of atoms, so try this next.
                formula0[element] = int(float(count))
    formula = OrderedDict(sorted(list(formula0.items()), key=lambda x: GetAtomicNumber(x[0])))
    mass = GetMass(formula)
    try:
        phase = lines[0][44]
        if phase.upper() != 'G':
            logging.warning('Specie %s is in %s phase' % (species, phase))
    except Exception as err:
        logging.error(err)
        print(lines[0])
    # Extract the NASA polynomial coefficients
    # Remember that the high-T polynomial comes first!
    Tmin = float(lines[0][45:55].strip())
    Tmax = float(lines[0][55:65].strip())
    Tint = float(lines[0][65:75].strip())
    coeffHigh = [float(lines[1][i*15:(i+1)*15]) for i in range(5)]
    coeffHigh += [float(lines[2][i*15:(i+1)*15]) for i in range(2)]
    coeffLow = [float(lines[2][i*15:(i+1)*15]) for i in range(2,5)]
    coeffLow += [float(lines[3][i*15:(i+1)*15]) for i in range(4)]

    data = {
        'sp': species,
        'nT': 2,
        'refcode': refcode,
        'formula': formula,
        'phase': phase,
        'm0': mass,
        'Trange': [[Tmin, Tint], [Tint, Tmax]],
        'Tcommon': Tint,
        'coeff': [coeffLow, coeffHigh]
    }
    if name is not None:
        data['name'] = name
    if H0 is not None:
        data['H0'] = H0
    if S0 is not None:
        data['S0'] = S0
    if ref is not None:
        data['ref'] = ref
    if CAS is not None:
        data['CAS'] = CAS

    return data


def ReadThermoBlock(content):
    """
    Read THERMO block of Chemkin

    Arguments:
        content: a list of str representing the content of file. content[0] should
                 be "THERMO" or "THERMO ALL". content[-1] should be "END"

    Return:
        data: a dict representing the information like the following
                Trange: [300.000, 1000.000, 5000.000]  // This is the temperature range
                (CH2O)3:
                    nT: 2
                    refcode: '70590'
                    formula: {H: 6, C: 3, O: 3}
                    phase: G
                    m0: 90.078
        Trange: Temperature ranges for 2 sets of coefficients: lowest T, common T, and highest T
                This one only exists if
                - content is from a pure thermo file
                - content is from gas phase or surf inp and content[0] == 'THERMO ALL'
    """
    c = content
    i = 0
    regCommentLine = re.compile('^!')
    regEndLine = re.compile('^END')
    # Loop to the start of 'THERMO' block
    while regCommentLine.match(c[i]):
        i += 1
    if not re.match('^THERMO', c[i]):
        raise ValueError('This is not THERMO Block')
    else:
        i += 1
        while regCommentLine.match(c[i]):
            i += 1
        w = re.findall('^\s*(\d+(?:\.\d*)?|\.\d+)\s*(\d+(?:\.\d*)?|\.\d+)\s*(\d+(?:\.\d*)?|\.\d+)\s*$',c[i])
        if w:
            Trange = [float(i) for i in w[0]]
            i += 1
            while regCommentLine.match(c[i]):
                i += 1
            i -= 1
        else:
            Trange = [300.000, 1000.000, 5000.000]  # default value

    # Loop all data block
    nsp = 0
    regRef = re.compile('Source:(.*)$')
    regH0 = re.compile('H0\(298K\)\s*=\s*([+-]?\d+\.?\d*)')
    regS0 = re.compile('S0\(298K\)\s*=\s*([+-]?\d+\.?\d*)')
    regCAS = re.compile('CAS="(.*)"')
    regNAME = re.compile('COMMON_NAME="(.*)"')
    regC = [regRef, regH0, regS0, regCAS, regNAME]

    data = dict()
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
        dat = ReadThermoEntry(dataBlock,
                              ref=auxopt[0], H0=auxopt[1], S0=auxopt[2], CAS=auxopt[3], name=auxopt[4])
        # print('Specie: %s' % (dat['sp']))
        nsp += 1

        spName = dat.pop('sp', None)
        data[spName] = dat

    return data, Trange


def ReadThermoFile(filepath):
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

    if 'inp' in fext or 'dat' in fext:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.readlines()
        data, Trange = ReadThermoBlock(content)
    else:
        with open(filepath, 'r', encoding='utf-8') as f:
            if 'yaml' in fext or 'yml' in fext:
                import oyaml as yaml
                thermo = yaml.load(f)
            elif 'json' in fext:
                import json
                thermo = json.load(f)
            else:
                raise TypeError('Unknown file type: %s' % (filepath))
        if isinstance(thermo, list):
            # old format
            data = {i['sp']:i for i in thermo}
            for i in data.values():
                i.pop('sp', None)
        else:
            data = thermo
        Trange = data.pop('Trange', None)
    return data, Trange

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reader for Chemkin Thermo file")

    parser.add_argument('-i', '--input', default='database\\therm.dat', help="Path of thermo input file")
    parser.add_argument('-y', '--yaml', default='database\\therm.yaml', help="Path of yaml output file")
    parser.add_argument('-j', '--json', default='database\\therm.json', help="Path of json output file")

    ARGS = parser.parse_args()

    data, Trange = ReadThermoFile(ARGS.input)

    if ARGS.json:
        # in json, we use old format
        import json
        jsonData = [{'Trange': Trange}]
        for key,value in data.items():
            jsonData.append({**{"sp":key}, **value})
        with open(ARGS.json, 'w', encoding='utf-8') as f:
            json.dump(jsonData, f, indent=2)
        import sys
        path2 = '\\'.join(os.getcwd().split('\\')[0:-1] + ['Lewis'])
        sys.path.append(path2)
        from loadthermo import PrettyPrint
        z = PrettyPrint(os.path.join(os.getcwd(), 'database\\therm.json'), True)

    if ARGS.yaml:
        import oyaml as yaml
        yamlData = {**{'Trange': Trange}, **data}
        with open(ARGS.yaml, 'w', encoding='utf-8') as f:
            yaml.dump(yamlData, f, encoding='utf-8')