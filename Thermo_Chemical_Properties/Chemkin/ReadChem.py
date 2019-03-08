#!/usr/bin/env python3
# encoding: utf-8
import argparse
import re
import numpy as np
import periodic
import logging
from scipy import constants
from ReadThermo import ReadThermoFile, ReadThermoBlock, Chemkin7Par


def ReadGasInput(filename, yamlOutput=None, thermoDefaultFile=True):
    '''
    Read GAS PHASE input file, The file usually has ELEMENTS,
    SPECIES, THERMO and REACTIONS blocks.

    Arguments:
        filename:  path of the gas phase input file
        yamlOutput (default: None): path of yaml file to output
        thermoDefaultFile: path of additional pure therm.inp

    Return:
        block: a dict containing information
    '''
    if thermoDefaultFile:
        thermoDefault,_ = ReadThermoFile(thermoDefaultFile)
    else:
        thermoDefault = None

    lines = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            # strip comment
            line = re.sub('!.*', '', line)
            if len(line.strip()) > 0:
                lines.append(line)

    # divide content into blocks
    n = len(lines)
    i = 0
    blockHead = ['ELEMENTS', 'SPECIES', 'THERMO', 'REACTIONS','SITE', 'BULK']
    blockHeadReg = [re.compile('^%s' % (i)) for i in blockHead]
    blockEndReg = re.compile('^END')
    block = {i: [] for i in blockHead}
    while i < n:
        for name, reg in zip(blockHead, blockHeadReg):
            if reg.match(lines[i]):
                print('Read in block: %s' % (name))
                block[name].append(lines[i])
                i += 1
                while not blockEndReg.match(lines[i]):
                    block[name].append(lines[i])
                    i += 1
                block[name].append(lines[i])
                print('Finish block: %s' % (name))
                i += 1
                break

    # analyze THERMO BLOCK
    thermoRaw = block['THERMO']
    if thermoRaw:
        thermo,_ = ReadThermoBlock(thermoRaw)
        block['THERMO'] = {'Data': thermo, 'Raw': thermoRaw}
    else:
        thermo = None

    # analyze REACTION BLOCK
    reacRaw = block['REACTIONS']
    if reacRaw:
        if thermo:
            if thermoDefault:
                reac = ReadReactionBlock(reacRaw[:], thermo={**thermo, **thermoDefault})
            else:
                reac = ReadReactionBlock(reacRaw[:], thermo)
        else:
            if thermoDefault:
                reac = ReadReactionBlock(reacRaw[:], thermoDefault)
            else:
                raise FileNotFoundError('No thermodynamic data found')

        elem = reac.pop('elements', None)
        sps = reac.pop('species', None)
        unit = reac.pop('unit', None)
        block['REACTIONS'] = {'elements': elem, 'species': sps, 'unit': unit, 'Data': reac, 'Raw': reacRaw}

    if isinstance(yamlOutput, str):
        import oyaml as yaml
        print('Write yaml file: %s' % (yamlOutput))
        with open(yamlOutput, 'w', encoding='utf-8') as f:
            yaml.dump({'unit': unit}, stream=f, encoding='utf-8')
            yaml.dump({'elements': elem}, stream=f, encoding='utf-8')
            yaml.dump({'species': sps}, stream=f, encoding='utf-8')
            yaml.dump(reac, stream=f, encoding='utf-8')

    return block


def ReadReactionBlock(lines, thermo=None):
    '''
    Read reaction block of chemkin inp
    '''
    if not thermo:
        thermo, Trange = ReadThermoFile('database\\therm.yaml')

    eUnits = ['CAL/MOLE', 'KCAL/MOLE', 'JOULES/MOLE', 'KJOULES/MOLE', 'KELVINS', 'EVOLTS']
    from scipy import constants
    e2cal = [1.0, 1.0e3, 1.0/constants.calorie, 1.0e3/constants.calorie,
             constants.k/constants.calorie*constants.N_A, constants.eV/constants.calorie*constants.N_A]
    eConvertFac = [[1.0*i/j for j in e2cal] for i in e2cal]
    eConvert = lambda x,i,j='CAL/MOLE': x*eConvertFac[eUnits.index(i.upper())][eUnits.index(j.upper())]

    # handle unit
    eUnit = 'CAL/MOLE'
    molUnit = 'MOLES'
    l1 = lines[0].split()
    if len(l1) == 3:
        molUnit = l1[2]
    if len(l1) > 2:
        eUnit = l1[1]
    try:
        e2cgs = eConvert(1.0, eUnit)
        if molUnit == 'MOLECULES':
            e2cgs /= constants.N_A
        elif molUnit != 'MOLES':
            raise ValueError('Wrong molecule unit')
    except:
        raise ValueError('Wrong units: %s and %s' % (eUnit, molUnit))

    regNumber = '[-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?'
    reacReg = re.compile('([\w\+\-\*()]+\s*(?:=|=>|<=>)\s*[\w\+\-\*()]+)\s+(%s)\s+(%s)\s+(%s)' % (regNumber, regNumber, regNumber))

    reac = dict()
    elem = set()
    sp = set()
    reac['unit'] = {'eUnit': eUnit, 'molUnit': molUnit, 'e2cgs': e2cgs}
    reac['species'] = []
    reac['elements'] = []
    for line in lines[1:]:
        m1 = reacReg.findall(line)
        if m1:
            reacName = m1[0][0]
            reac[reacName] = {'info': ParseChemkinReaction(reacName, thermo),
                              'Ar': [float(i) for i in m1[0][1:]]}
            thisElem = reac[reacName]['info'].pop('elem', None)
            elem.update(set(thisElem))

            thisSp = reac[reacName]['info'].pop('sp', None)
            sp.update(set(thisSp))

            lastRK = m1[0][0]
        else:
            nSlash = line.count('/')
            if nSlash % 2 == 0:
                propStr = [i.strip() for i in line.split('/') if i.strip()]
                for propName, value in zip(propStr[0::2], propStr[1::2]):
                    if propName in ['TROE', 'SRI', 'HIGH', 'LOW']:
                        reac[lastRK][propName] = [float(i) for i in value.split()]
                    else:
                        if 'T' not in reac[lastRK].keys():
                            reac[lastRK]['T'] = {}
                        reac[lastRK]['T'][propName] = value
    reac['species'] = sorted(list(sp))
    reac['elements'] = sorted(list(elem), key=lambda x: periodic.GetAtomicNumber(x))
    return reac


def ParseChemkinReaction(r, thermo):
    '''
    Parse chemkin format reaction

    Arguments:
       r:  a string of reaction formula
       thermo:  a dict of thermo database (used to find elements of reactants and products)
                like the following one in yaml format:
                (CH2O)3:
                    formula: {C: 3, H: 6, O: 3}
                (CH3)2SICH2:
                    formula: {C: 3, H: 8, SI: 1}

    '''
    if not isinstance(thermo, dict):
        raise ValueError('Wrong type of thermo')
    try:
        spDatabase = {i:thermo[i]['formula'] for i in thermo.keys()}
    except:
        raise ValueError('Wrong type of thermo')
    # remove third-body symbol
    if '(+M)' in r:
        thirdBody = 'PRESSURE'
        r = r.replace('(+M)', '')
    elif re.findall('\+M(?!\w)', r):
        thirdBody = 'PART'
        r = re.sub('\+\s*M(?!\w)', '', r)
    else:
        thirdBody = 'NONE'

    # split reactants and products
    reversible = True
    reacStr, prodStr = r.split('=')
    if '<=>' in r:
        reacStr = reacStr[:-1].strip()
        prodStr = prodStr[1:].strip()
    elif '=>' in r:
        prodStr = prodStr[1:].strip()
        reversible = False

    spStrs = [reacStr, prodStr]
    spDicts = []

    # replace
    for spStr in spStrs:
        # substitute plus sign in cation
        spStr2 = re.sub('(\d?)\+\s*$','\\1_plus', spStr)
        spStr2 = re.sub('(\d?)\+(?=\s*\+)', '\\1_plus',spStr2)

        spsStr = [i.strip().replace('_plus','+') for i in spStr2.split('+')]
        spsCount = [1 for i in spsStr]
        spsName = spsStr[:]
        for i, j in enumerate(spsStr):
            k = re.findall('^(\d(?:\.\d+)?)\s*', j)
            if k:
                spsCount[i] = float(k[0])
                spsName[i] = re.sub('^(\d(?:\.\d+)?\s*)', '', j)

        spsFormula = [None] * len(spsName)
        element = dict()
        for i,(iname,icount) in enumerate(zip(spsName,spsCount)):
            if iname in spDatabase.keys():
                spsFormula[i] = spDatabase[iname]
                for j,k in spsFormula[i].items():
                    if j in element.keys():
                        element[j] += icount * k
                    else:
                        element[j] = icount*k
            else:
                logging.warning('Specie %s is not in thermo database, it might be an etchant' % (iname))

        spDicts.append({'name': spsName,
                        'count': spsCount,
                        'formula': spsFormula,
                        'element': element})

    reac = spDicts[0]
    prod = spDicts[1]

    for i,j in enumerate(reac['element'].keys()):
        if reac['element'][j] != prod['element'][j]:
            raise ValueError('Reaction %s is not balanced' % (r))

    elem = list(reac['element'].keys())
    sp = list(set(reac['name'] + prod['name']))

    reac.pop('formula', None)
    prod.pop('formula', None)
    reac.pop('element', None)
    prod.pop('element', None)

    return {'reac': reac, 'prod': prod, 'elem': elem, 'sp': sp,
            'reversible': reversible, 'thirdBody': thirdBody}


def EquilibriumConstant(reaction, Tin, thermo, thermEvaluator='Chemkin7Par'):
    '''
    Calculate Equilibrium constant at temperature T
    '''
    if isinstance(Tin, list):
        T2 = np.array(Tin, dtype=np.float64)
    elif isinstance(Tin, (float, int)):
        T2 = np.array([Tin], dtype=np.float64)
    elif not isinstance(Tin, np.ndarray):
        raise ValueError('Wrong dtype for T')
    else:
        T2 = Tin

    if thermEvaluator == 'Chemkin7Par':
        thermEval = Chemkin7Par

    Rgas = constants.R / 1000  # kJ/mol/K
    kb = constants.k

    Kp = np.zeros((len(T2),), dtype=np.float64)
    Kc = np.zeros((len(T2),), dtype=np.float64)

    parsedReac = ParseChemkinReaction(reaction, thermo)
    reactant = parsedReac['reac']
    product = parsedReac['prod']

    Gdiff = np.zeros(T2.shape, np.float64)   # kJ/mol
    eta = np.zeros(T2.shape, np.float64)

    for sp, isto in zip(reactant['name'], reactant['count']):
        _,_,G,_ = thermEval(sp, T2, thermo, 'si')
        Gdiff -= G
        eta -= isto

    for sp, isto in zip(product['name'], reactant['count']):
        _,_,G,_ = thermEval(sp, T2, thermo, 'si')
        Gdiff += G
        eta += isto

    Kp = np.exp(-Gdiff / Rgas / T2)  # for pressure in atm

    # Kp = Pi_i   [(n_i kb T)*(Pa to atm)]^sto_i
    Kc = Kp / (kb * T2 / constants.atm)**eta   # concentation in molecule/m^3

    return Kp, Kc

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reader for Chemkin chemical input file')

    parser.add_argument('-i', '--input', default='database\\chem.inp', help="Path of chem input")
    parser.add_argument('-y', '--yaml', default='database\\reaction.yaml', help='Path of yaml output')
    parser.add_argument('-t', '--therm', default='database\\therm.yaml', help='Path of therm input')

    ARGS = parser.parse_args()
    w = ReadGasInput(ARGS.input, yamlOutput=ARGS.yaml, thermoDefaultFile=ARGS.therm)
