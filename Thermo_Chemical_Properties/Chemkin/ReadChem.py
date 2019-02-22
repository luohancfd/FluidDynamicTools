#!/usr/bin/python3
from chempy.chemistry import Species, Reaction
from chempy.util import periodic, parsing
import itertools
import warnings
import re
from enum import Enum


# ========  Generate gas-phase reaction data ======================
# class ReactionRate(object):
#     '''
#     class definition of reaction rate
#     '''
#     availableKeys = ['Ar', 'HIGH', 'LOW', 'TROE', 'SRI', 'T']
#     def __init__(self, rate):
#         '''
#         rate: a dict define reaction rates, examples
#              Ar: [1.698E+16, 0.0, 84840.0]
#              HIGH: [1.698E+16, 0.0, 84840.0]
#              LOW: [1.698E+16, 0.0, 84840.0]
#              TROE: [1, 2, 3]
#              SRI: [1, 2, 3]
#              T:                  !Third body
#                 - H2: 2
#         '''
#         for i in availableKeys:
#             if i in rate.keys():
#                 setattr(self, i, rate[i])
#         if 'e_unit' in rate.keys():
#             self.e_unit = rate.e_unit
#         else:
#             self.e_unit = 'cal/mol'

#     def __repr__(self):
#         '''
#         print auxilary reaction rate
#         '''




class ChemkinReaction(object):
    '''
    Definition of reaction class for chemkin gas phase
    '''
    def __init__(self, string, rate):
        '''
        Initiation of the class

        string: chemical formula of the reaction
        rate: a dict containing different method of description
        '''
        TBODY_NONE=0
        TBODY_PART=1
        TBODY_PRESSURE=2
        if '<=>' in string or '=' in string:
            self.backward = True
        else:
            self.backward = False

        # clean string
        chempyString = string.replace('<=>', '->').replace('=>','->').replace('=','->')
        if '(+M)' in chempyString:
            self.thirdBody = TBODY_PRESSURE
            chempyString.replace('(+M)', '')
        elif re.findall('\+M(?!\w)', chempyString):
            self.thirdBody = TBODY_PART
            chempyString = re.sub('\+\s*M(?!\w)', '', chempyString)
        else:
            self.thirdBody = TBODY_NONE

        self.rk = Reaction.from_string(chempyString)
        self.sp = [i for i in list(self.rk.keys())]
        self.elem = list(set(itertools.chain.from_iterable([list(Species.from_formula(i).composition.keys()) for i in self.sp])))
        self.elem.sort()



def ReadGasInput(filename, yamlOutput=None):
    '''
    Read gasphase input file
    '''
    lines = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = re.sub('!.*', '', line).strip()
            if len(line) > 0:
                lines.append(line)

    # divide content into blocks
    n = len(lines)
    i = 0
    blockHead = ['ELEMENTS', 'SPECIES', 'THERMO', 'REACTIONS']
    blockHeadReg = [re.compile('^%s'%(i)) for i in blockHead]
    blockEnd = ['END' for i in blockHead]
    blockEndReg = re.compile('^END')
    block = {i: [] for i in blockHead}
    while i < n:
        for name, reg in zip(blockHead, blockHeadReg):
            if reg.match(lines[i]):
                print('Readin block: %s'%(name))
                i += 1
                while not blockEndReg.match(lines[i]):
                    block[name].append(lines[i])
                    i += 1
                print('Finish block: %s'%(name))
                i += 1
                break

    # process reactions block
    regNumber = '[-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?'
    reacReg = re.compile('([\w\+\-\*()]+\s*(?:=|=>|<=>)\s*[\w\+\-\*()]+)\s+(%s)\s+(%s)\s+(%s)' % (regNumber, regNumber, regNumber))

    reac = dict()
    for line in block['REACTIONS']:
        m1 = reacReg.findall(line)
        if m1:
            reac[m1[0][0]] = {'Ar': [float(i) for i in m1[0][1:]]}
            lastRK = m1[0][0]
        else:
            nSlash = line.count('/')
            if nSlash % 2 == 0:
                propStr = [i.strip() for i in line.split('/') if i.strip()]
                for propName, value in zip(propStr[0::2], propStr[1::2]):
                    if propName in ['TROE', 'SRI', 'HIGH', 'LOW']:
                        reac[lastRK][propName] = [float(i) for i in value.split()]
                    else:
                        if not 'T' in reac[lastRK].keys():
                            reac[lastRK]['T'] = {}
                        reac[lastRK]['T'][propName] = value

    if yamlOutput:
        import yaml
        with open(yamlOutput, 'w', encoding='utf-8') as f:
            yaml.dump(reac, stream=f, encoding='utf-8')

    block['REACTIONS'] = reac
    return block
