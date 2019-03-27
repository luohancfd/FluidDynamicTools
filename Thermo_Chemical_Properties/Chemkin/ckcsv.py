#!/usr/bin/env python
# coding: utf-8

import os
import csv
import re
from glob import glob
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl

def ReadCKCSVBlock(rows, nx, type=1, ny=1):
    if type == 1:
        varName = rows[0].strip()
        unit = rows[1].strip('() ')
        data = np.array([float(i) for i in rows[2:nx+2]])
        return varName,{'unit':unit, 'data': data}
    if type == 2:
        varName = rows[0][0].strip()
        unit = rows[0][1].strip('() ')
        data = []
        for i in range(ny):
            data.append([float(i) for i in rows[1 + i][:nx]])
        data = np.array(data, np.float64)
        return varName,{'unit':unit, 'data': data}



def Read1DCKCSV(fileName, delimiter=','):
    '''
    Read CKCSV file starting with 'label,...'

    Arguments:
        fileName: name of the ckcsv file
        delimiter: delimiter in ckcsv file, usually comma symbol
    '''

    _,ext = os.path.splitext(fileName)
    if ext == '':
        fileName += '.ckcsv'

    rawData = []
    with open(fileName, newline='') as f:
        csvReader = csv.reader(f, delimiter=delimiter, dialect='excel')
        for row in csvReader:
            rawData.append(row)

    if 'label' not in rawData[0][0]:
        raise ValueError('File {:s} is not loaded with correct loader'.format(fileName))

    result = []
    nLabel = 0
    nRow = len(rawData)
    iRow = 0
    while iRow < nRow:
        if rawData[iRow][0] == 'label':
            # label row
            name = rawData[iRow][1].strip()
            nx = int(rawData[iRow][2].strip())
            nLabel += 1
            prop = {'name': name, 'nx': nx, 'ntype': '1D'}
            iResult = {}
            iRow += 1
            while 'label' not in rawData[iRow][0]:
                varName, d = ReadCKCSVBlock(rawData[iRow], nx, type=1)
                iResult[varName] = d
                iRow += 1
                if (iRow == nRow):
                    break

            result.append({'prop': prop, 'data':iResult})
        elif rawData[iRow][0] == 'label_2D':
            # label 2D row
            nRun = int(rawData[iRow][1].strip())
            nRunL = int(rawData[iRow][3].strip())

            nx = int(rawData[iRow][2].strip())
            nxL = int(rawData[iRow][4].strip())

            name = rawData[iRow][5].strip()
            prop = {'name':name, 'nx': nx, 'ntype': '2D', 'nrun':nRun}

            iResult = {}

            iRow += 1
            for i in range(nRunL):
                varName, d = ReadCKCSVBlock(rawData[iRow], nRun)
                iResult[varName] = d
                iRow += 1

            for i in range(nxL):
                varName, d = ReadCKCSVBlock(rawData[iRow], nx)
                iResult[varName] = d
                iRow += 1

            while 'label' not in rawData[iRow][0]:
                varName, d = ReadCKCSVBlock(rawData[iRow:iRow+nRun+1], nx, type=2, ny=nRun)
                iResult[varName] = d
                iRow += nRun+1
                if (iRow == nRow):
                    break

            result.append({'prop': prop, 'data':iResult})
        else:
            raise ValueError('Wrong')
    return result, rawData


def ReadPSR(projectName, postfix='vs_parameter', wkdir='./'):
    fileName = '_'.join(['CKSoln', projectName, postfix])
    filePath = os.path.join(wkdir, fileName+'.ckcsv')
    data, rawData = Read1DCKCSV(filePath)
    data = [i['data'] for i in data]
    for d in data:
        d = {re.sub('_Soln#\d+','', key): val for key,val in d.items() }


    if postfix == 'vs_parameter':
        ddata = data[0]
        nRun = int(ddata['Run_number']['data'][-1])
        ldata = [{key: val['data'][i] for key, val in ddata.items()} for i in range(nRun)]
        units = {key: val['unit'] for key,val in ddata.items()}
        return ldata, units, ddata


def ReadPSRCVD(wkdir, jobName, varExtract=None, groupvar=None):
    ldata,unit,ddata = ReadPSR(jobName, wkdir=wkdir)

    if groupvar is 'S2V':
        for i,idata in enumerate(ldata):
            idata['S2V'] = idata['Internal_Surface_Area_C1__PSR_PSR_(C1)']/idata['Volume']
        ddata['S2V'] = {'data': ddata['Internal_Surface_Area_C1__PSR_PSR_(C1)']['data'] / ddata['Volume']['data'],
                        'unit': ''}

    # Write something we need to a new file
    fieldWrite = []
    for key in ldata[0].keys():
        if any(x in key for x in ['Run_number','VolumetricFlow', 'Growth_rate', 'Site_fraction', 'Residence','S2V',
                                  'Inlet_flow_rate', 'Exit_mass_flow', 'Internal_Surface','Volume']):

            fieldWrite.append(key)
    with open(os.path.join(wkdir, '{:s}_run_result.csv'.format(jobName)), 'w', newline='') as f:
        csvWriter = csv.DictWriter(f, fieldWrite,extrasaction='ignore')
        csvWriter.writeheader()
        csvWriter.writerow(unit)
        for idata in ldata:
            csvWriter.writerow(idata)


    if varExtract is None:
        varExtract = list(ddata.keys())
        if groupvar is not None:
            varExtract.remove(groupvar)

    if groupvar is not None and isinstance(groupvar,str):
        groupSet = np.unique(np.sort(ddata[groupvar]['data']))
        result = [{groupvar:i} for i in groupSet]
        for j,i in enumerate(groupSet):
            index = np.where(ddata[groupvar]['data'] == i)
            for var in varExtract:
                result[j][var] = ddata[var]['data'][index]
        return result
    else:
        return {key: val['data'] for key,val in ddata.items()}
