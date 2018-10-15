# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:43:37 2018

@author: Han
"""
from shutil import copy2
import os
import re
from scipy.constants import N_A,pi,Boltzmann
import os,scipy.constants
import numpy as np
from tec2py import ReadTecplotDat

eV = scipy.constants.physical_constants['electron volt'][0]
eV2K = eV/Boltzmann

ArrheniusConstant = {'N2+O':[8.9340e-7*N_A,  -0.4807,   113950.0, 14, 16, 3352,'N2','O',113950],
                     'N2+N': [6.6831e-6*N_A,  -0.6996,   113200.0, 14, 14, 3352,'N2','N',113200],
                     'O2+O': [      2.5e18,   -0.565,   60491.819, 16, 16, 2273.54,'O2','O',5.21275*eV2K],
                     'N2+N2': [   4.5e-6*N_A,  -0.675,    117000.0, 14, 14, 3353,'N2','N2',117000.],
                     'O2+O2': [   2.0e+21   ,  -1.5  ,     59500.0, 16, 16, 2273.54,'O2','O2',59380+(5.21275-5.1153)*eV2K],
                     'O2+N2':[8.132E-10*N_A*10   ,  -0.131  ,     59380.0, 16, 14, 2273.54,'O2','N2',5.21275*eV2K],
                     'N2+O2':[1.9240e+17,       -0.5000,1.1317e+05  , 14, 16, 3353 , 'N2','O2',9.82163*eV2K],
                     'NO+N2':[5.0E15   ,  0 ,     75500, 16, 14, 2739.72,'NO','N2',6.55879*eV2K]}

def nparray2string(val):
    '''
    convert a numpy array to string
    '''
    return ' ' + np.array2string(val, threshold=np.nan,
                                 max_line_width=np.inf,
                                 formatter={'int': lambda x: '%d' % (x)}).replace('[', '').replace(']', '') + '\n'



def CombineVrate(sps,spfolder):
    #spfolder = 'N2O'
    #sps=['N2+O']
    ArRate = []
    fdia = []
    react_sp=[]
    for sp in sps:
        ArRate.append(ArrheniusConstant[sp])
        partner = sp.split('+')
        if partner[0][1] == '2':
            ps = [partner[0][0], partner[0][0]]
            partner[0] = partner[0].replace('2','<sub>2</sub>')
        else:
            ps = [partner[0][0], partner[0][1]]

        if len(partner[1]) == 2:
            fdia.append(True)
            if partner[1][1] == '2':
                partner[1] = partner[1].replace('2','<sub>2</sub>')
        else:
            fdia.append(False)
        react_sp.append([partner[0],partner[1],ps[0],ps[1],partner[1]])
    # symbols of reaction
    # react_sp2 = [[j.replace('<sub>2</sub>','2')   for j in i ]for i in react_sp]



    rootfolder = next(os.walk(spfolder))[1]
    rootfolder = [i for i in rootfolder if re.match('T=\d+$',i)]
    # list of temperatures
    T = [float(i.split('=')[1]) for i in rootfolder]

    index = sorted(list(range(len(T))), key=lambda i: T[i])
    T = [T[i] for i in index]
    rootfolder = [os.path.join(spfolder,rootfolder[i]) for i in index]
    eqrate = np.zeros((len(T),len(sps)),dtype=np.float64)
    qctrate = eqrate.copy()

    with open(os.path.join(spfolder,'DSMCVrate%s.dat'%(spfolder)),'w',encoding='utf8') as f:
        for ifd, ifolder in enumerate(rootfolder):
            # load vrate
            (data,VARlist) = ReadTecplotDat(os.path.join(ifolder,'IMF_Vrate.dat'))
            nreac = len(data)
            reac_str = ['' for i in range(len(data))]
            with open(os.path.join(ifolder,'IMF_Vrate.dat'),'r') as f1:
                ireac = 0
                for i,line in enumerate(f1):
                    if i == 0:
                        headerc = line
                    if re.match('^#.*->.*$',line):
                        reac_str[ireac] = line
                        ireac += 1
            Tv = float(re.findall('VT:\s*(\d+(?:\.\d+)?)',headerc)[0])

            # load eqrate
            w = np.loadtxt(os.path.join(ifolder,'DS1REAC.DAT'))
            CT =  w[-1][1]
            for i in range(nreac):
                eqrate[ifd,i] = w[-1][ 7 + nreac + i]

            for i in range(nreac):
                qctrate[ifd,i] = ArRate[i][0]*CT**ArRate[i][1]*np.exp(-ArRate[i][2]/CT)

            for i,idata in enumerate(data):
                Vrate = idata['data'][:,VARlist.index('Rate (cm3/mol/s)')]
                idata['data'] = np.vstack((idata['data'].T,Vrate/eqrate[ifd,i]*qctrate[ifd,i])).T

            VARlist.append('Rate-scaled (cm3/mol/s))')
#            for ireac,idata in enumerate(data):
#                Ev = idata['data'][:,VARlist.index('Evib(eV)')]
#                Vrate = idata['data'][:,VARlist.index('Rate (cm3/mol/s)')]
#
#                parfun = np.exp(-Ev*eV2K/Tv)
#                eq_rate = Vrate*parfun
#                eq_rate = np.sum(eq_rate)/np.sum(parfun)

            # write to file
            if (ifd == 0):
                f.write('VARIABLES = ')
                for i in VARlist[:-1]:
                    f.write('"%s",'%(i))
                f.write('"%s"\n'%(VARlist[-1]))
            f.write('\n'+headerc)
            for i in range(len(data)):
                f.write(reac_str[i])
                f.write('ZONE I = %d, '%(data[i]['data'].shape[0]))
                f.write('T = "%s T=%dK %s"'%(sps[i],T[ifd],data[i]['title']))
                if data[i]['passivevarlist']:
                    f.write(', PASSIVEVARLIST = %s\n'%(data[i]['passivevarlist'].__repr__()))
                else:
                    f.write('\n')
                for iline in data[i]['data']:
                    f.write('%12.6f'%(iline[0]))
                    f.write('  %3d'%(iline[1]))
                    for j in range(2,len(iline)):
                        f.write('  %14.6E'%(iline[j]))
                    f.write('\n')

#            with open(os.path.join(ifolder,'IMF_Vrate.dat'),'r') as f1:
#                for iline in f1:
#                    if re.match('^\s*ZONE.*T=',iline):
#                        f.write(re.sub(r'T\s*=\s*["''](.*)[''"]',r'T="%s T=%dK \1"'%(spfolder,T[ifd]),iline))
#                    elif re.match('^\s*VARIABLES',iline):
#                        if ifd == 0:
#                            f.write(iline)
#                    else:
#                        f.write(iline)


def CalculateEqFromVrate(spfolder,fname=None):
    if not fname:
        fname = 'DSMCVrate%s.dat'%(spfolder)
    (data,VARlist) = ReadTecplotDat(os.path.join(spfolder,fname))
    for izone,szone in enumerate(data):
        T = float(re.findall('T=(\d+(?:\.\d+)?)',szone['title'])[0])
        Ev = szone['data'][:,VARlist.index('Evib(eV)')]
        Vrate = szone['data'][:,VARlist.index('Rate (cm3/mol/s)')]
        Parfun = np.exp(-Ev*eV2K/T)
        Eqrate = np.sum(Vrate*Parfun)/np.sum(Parfun)




if __name__ == "__main__":
   #  os.chdir('E:\\Research\\NonequilibriumGas\\Macheret-Fridmann\\DSMC_MF\\AHO')
     CombineVrate(sps=['N2+O'],spfolder='N2O')
     CombineVrate(sps=['O2+O'],spfolder='O3')
     CombineVrate(sps=['N2+N'],spfolder='N3')

     CombineVrate(sps=['O2+N2','N2+O2'],spfolder='O2N2')
     CombineVrate(sps=['NO+N2'],spfolder='NON2')
     CombineVrate(sps=['O2+O2'],spfolder='O4')
#
#     os.chdir('E://Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/Eq/')
#     ExtractEqrate(ArrheniusConstant,sps=['O2+N2','N2+O2'],spfolder='O2N2Double')
#     ExtractEqrate(ArrheniusConstant,sps=['NO+N2'],spfolder='NON2')
#     ExtractEqrate(ArrheniusConstant,sps=['O2+O2'],spfolder='O4_Byron')
#     ExtractEqrate(ArrheniusConstant,sps=['N2+N2'],spfolder='N4_2')
