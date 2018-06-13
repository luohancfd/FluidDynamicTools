# -*- coding: utf-8 -*-
"""
Created on Sat May  5 15:56:56 2018

@author: Han
"""

import numpy as np
from openpyxl import Workbook


def py2dat(tdata,fname):
    if not isinstance(tdata,dict):
        raise TypeError
    else:
        with open(fname,'w',encoding='utf-8') as f:
            # variables
            f.write('VARIABLES = ')
            for var in tdata['varnames'][:-1]:
                f.write('"%s", '%(var,))
            f.write('"%s"\n'%(tdata['varnames'][-1],))
            
            for iline,line in enumerate(tdata['lines']):
                x = np.asarray(line['x'])
                lenx = x.size
                data = x.reshape(1,lenx)
                
                y = np.asarray(line['y']).reshape(1,lenx)
                data = np.vstack((data,y))
                
                if 'z' in line.keys():
                    if line['z'].size > 0:
                        z = np.asarray(line['z']).reshape(1,lenx)
                        data = np.vstack((data,z))
                if 'v' in line.keys():
                    if line['v'].size > 0:
                        v = np.asarray(line['v'])
                        if v.ndim == 1:
                            v = v.reshape(1,lenx)
                        data = np.vstack((data,v))
                
                if 'zonename' in line.keys():
                    if len(line['zonename']) == 0:
                        zonename = 'ZONE %d'%(iline,)
                    else:
                        zonename = line['zonename']
                else:
                    zonename = 'ZONE %d'%(iline,)
                    
                ivarloc = 0        
                f.write('ZONE I = %d T="%s" DATAPACKING=POINT\n'%(lenx,zonename))
                f.write(' '+np.array2string(data.T,threshold=np.nan,max_line_width=np.inf).replace('[','').replace(']','') )
                f.write('\n\n')
            
def py2xls(tdata, fname, creator='Han Luo'):
    wb = Workbook()
    k =-1
    for line in tdata['lines']:
        k = k+1
        if 'zonename' in line.keys():
            if len(line['zonename']) == 0:
                zn = 'Sheet %d'%(k+1,)
            else:
                zn = line['zonename']
        else:
            zn = 'Sheet %d'%(k+1,)
            
        ws=wb.create_sheet(zn,k)
        for i,j in enumerate(tdata['varnames']):
            cellnum = chr(65+i)
            ws[cellnum+'1'] = j
        data = np.vstack((line['x'],line['y'],line['z'],line['v'])).T
        for i,row in enumerate(data):
            for j,ele in enumerate(row):
                cellnum = chr(65+j)
                cellid = cellnum+'%d'%(i+2,)
                ws[cellid] = ele
        
    wb.properties.creator = creator   
    wb.save(fname)                  
                            
                