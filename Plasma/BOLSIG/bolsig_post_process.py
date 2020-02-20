#!/usr/bin/env python3

import re
import numpy as np
import argparse
import sys
import os

# Counter counts the number of occurrences of each item
from collections import Counter
from itertools import tee, count


def indexOfSubstr(theList, substring):
    for i, s in enumerate(theList):
        if substring in s:
            return i
    return -1


def uniquify(seq, suffs=count(1)):
    """Make all the items unique by adding a suffix (1, 2, etc).

    `seq` is mutable sequence of strings.
    `suffs` is an optional alternative suffix iterable.
    """
    not_unique = [k for k, v in Counter(
        seq).items() if v > 1]  # so we have: ['name', 'zip']
    # suffix generator dict - e.g., {'name': <my_gen>, 'zip': <my_gen>}
    suff_gens = dict(zip(not_unique, tee(suffs, len(not_unique))))
    for idx, s in enumerate(seq):
        try:
            suffix = str(next(suff_gens[s]))
        except KeyError:
            # s was unique
            continue
        else:
            seq[idx] += suffix


def Matrix2Tecplot(varnames, data, zonename, style=None):
    '''
    export a string to write in tecplot file
    varnames: a list of variable names
    data: a numpy array of data
    style: a list of output style, the size should be the same as the number of columns in data
    zonename: name of the zone
    with open(fname, 'w', encoding='utf-8') as f:
    '''
    l = 'VARIABLES=' + ','.join(["\"{:s}\"".format(i)
                                 for i in varnames]) + '\n'
    try:
        m, n = data.shape
    except:
        m = len(data)
        n = len(data[0])
    l += 'ZONE I={:d} T="{:s}"\n'.format(m, zonename)
    if not style:
        style = ['{:e}'] * n
    else:
        for i, j in enumerate(style):
            if not j:
                style[i] = '{:e}'
    d = '\n'.join(['  '.join([s.format(i) for s, i in zip(style, j)])
                   for j in data])
    return l + d


class BOLSIG:
    '''
    Class to store BOLSIG+ output
    '''

    def __init__(self, filename):
        '''
        Read BOLSIG combined output file
        '''
        self.Erd = []
        self.muN = []
        self.ne = []
        self.DN = []
        with open(filename, 'r') as f:
            content = f.readlines()
        i = 0
        while i < len(content):
            if 'Transport coefficients' in content[i]:
                # scan variable names
                i += 1
                transportVar = []
                m = re.findall('^\s*A\d+\s+(.*?)\s*$', content[i])
                while m:
                    transportVar.append(m[0])
                    i += 1
                    if i < len(content):
                        m = re.findall('^\s*A\d+\s+(.*?)\s*$', content[i])
                    else:
                        break

                # rename variables with the same name
                uniquify(transportVar)

                # skip header line of transport table
                transportVar = ['R#', 'E/N (Td)'] + transportVar
                ErdIndex = transportVar.index('E/N (Td)')
                muNIndex = transportVar.index('Mobility *N (1/m/V/s)')
                i += 1

                # read transport table
                transport = []
                m = content[i].split()
                while len(m) == len(transportVar):
                    try:
                        m = [float(j) for j in m]
                    except:
                        mm = m[:]
                        m = []
                        for j in mm:
                            w = re.findall('([+-]?\d+(?:\.\d*)?)E?([+-]?\d+)?', j)
                            if w[0][1]:
                                m.append(float(w[0][0])*10**float(w[0][1]))
                            else:
                                m.append(float(w[0][0]))
                    self.Erd.append(m[ErdIndex])
                    self.muN.append(m[muNIndex])
                    self.ne.append(m[transportVar.index('Mean energy (eV)')])
                    self.DN.append(m[transportVar.index('Diffusion coefficient *N (1/m/s)')])
                    transport.append(m)
                    i += 1
                    if i < len(content):
                        m = content[i].split()
                    else:
                        break

                # add Te to the table, Te = 2/3*average_energy
                transportVar += ['Te (eV)']
                self.Te = []
                j = indexOfSubstr(transportVar, 'Mean energy')
                for ii, k in enumerate(transport):
                    transport[ii].append(k[j]*2/3)
                    self.Te.append(k[j]*2/3)

                transport = np.array(transport)

                # clean variable name
                transportVar = [re.sub('\s+', ' ', j) for j in transportVar]
            elif 'Rate coefficients' in content[i]:
                # scan variable names
                i += 1
                rateVar = []
                m = re.findall('^\s*C\d+\s+(.*?)\s*$', content[i])
                while m:
                    rateVar.append(m[0])
                    i += 1
                    if i < len(content):
                        m = re.findall('^\s*C\d+\s+(.*?)\s*$', content[i])
                    else:
                        break
                uniquify(rateVar)

                # skip header line of rate table
                rateVar = ['R#', 'E/N (Td)', 'Energy (eV)'] + rateVar
                i += 1

                # read transport table
                rate = []
                m = content[i].split()
                while len(m) == len(rateVar):
                    m = [float(j) for j in m]
                    rate.append(m)
                    i += 1
                    if i < len(content):
                        m = content[i].split()
                    else:
                        break

                # add Te
                rateVar += ['Te (eV)']
                j = indexOfSubstr(rateVar, 'Energy')
                for ii, k in enumerate(rate):
                    rate[ii].append(k[j]*2/3)
                rate = np.array(rate)

                # clean variable name
                rateVar = [re.sub('\s+', ' ', i) for i in rateVar]
            else:
                i += 1

        self.transport = {'varnames': transportVar, 'data': transport}
        self.rate = {'varnames': rateVar, 'data': rate}

    def writeRateTecplot(self, fp, zonename):
        s = [None] * len(self.rate['varnames'])
        s[0] = '{:.0f}'
        s[1] = '{:f}'
        s[2] = '{:f}'
        fp.write(Matrix2Tecplot(self.rate['varnames'], self.rate['data'], zonename, s))

    def writeTransportTecplot(self, fp, zonename):
        s = [None] * len(self.transport['varnames'])
        s[0] = '{:.0f}'
        s[1] = '{:f}'
        fp.write(Matrix2Tecplot(self.transport['varnames'], self.transport['data'], zonename, s))

    def writeMunCOMSOL(self, fp):
        fp.write('#Te(eV) muN(1/m/V/s)\n')
        prev = -1
        for i, j in zip(self.Te, self.muN):
            if '{:.6e}'.format(prev) != '{:.6e}'.format(i):
                prev = i
                fp.write('{:.6e} {:.6e}\n'.format(i, j))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Post process BOLSIG+ output")
    parser.add_argument('INPUT', help="BOLSIG+ output")
    ARGS = parser.parse_args()

    BOLSIG_file = ARGS.INPUT
    basename = os.path.basename(BOLSIG_file).split('.')[0]
    dirname = os.path.dirname(BOLSIG_file)
    if len(sys.argv) >= 3:
        zonename = sys.argv[2]
    else:
        zonename = basename

    bol = BOLSIG(BOLSIG_file)

    with open(os.path.join(dirname, '{:s}_transport.tec'.format(basename)), 'w', encoding='utf8') as f:
        bol.writeTransportTecplot(f, zonename)

    with open(os.path.join(dirname, '{:s}_rate.tec'.format(basename)), 'w', encoding='utf8') as f:
        bol.writeRateTecplot(f, zonename)

    with open(os.path.join(dirname, '{:s}_muN.txt'.format(basename)), 'w', encoding='utf8') as f:
        bol.writeMunCOMSOL(f)
