# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 22:30:23 2017

@author: Han
"""
from numpy import exp,sqrt
from scipy.constants import N_A,pi,Boltzmann
import os,scipy.constants
import numpy as np
from shutil import copy2
from distutils.dir_util import copy_tree,remove_tree


eV = scipy.constants.physical_constants['electron volt'][0]

ArrheniusConstant = {'N2+O':[8.9340e-7*N_A,  -0.4807,   113950.0, 14, 16, 3352,'N2','O'],
                     'N2+N': [6.6831e-6*N_A,  -0.6996,   113200.0, 14, 14, 3352,'N2','N'],
                     'O2+O': [      2.5e18,   -0.565,   60491.819, 16, 16, 2273.54,'O2','O'],
                     'N2+N2': [   4.5e-6*N_A,  -0.675,    117000.0, 14, 14, 3353,'N2','N2'],
                     'O2+O2': [   2.0e+21   ,  -1.5  ,     59500.0, 16, 16, 2273.54,'O2','O2'],
                     'O2+N2':[8.132E-10*N_A*10   ,  -0.131  ,     59380.0, 16, 14, 2273.54,'O2','N2'],
                     'N2+O2':[1.9240e+17,       -0.5000,1.1317e+05  , 14, 16, 3353 , 'N2','O2'],
                     'NO+N2':[5.0E15   ,  0 ,     75500, 16, 14, 2739.72,'NO','N2']}

def TecplotSub(file,datafile):
    with open(file,'r') as f:
        content = f.readlines()
    with open(file,'w') as f:
        for line in content:
            if '$!VarSet |LFDSFN1|' in line:
                f.write('$!VarSet |LFDSFN1| = \'%s\'\n'%(datafile,))
            else:
                f.write(line)
                
def MFrate(T,Tv,A,n,D,m,M,theta):
    b = 2
    alpha = (m/(M+m))**2;
    Ta = alpha*Tv+(1-alpha)*T
    fac = (1-exp(-theta/Tv))/(1-exp(-theta/T))
    deltaD = 3*b*alpha**2*D
    Dstar = D-deltaD
    # old
    # L2 = 9/64*sqrt(pi*(1-alpha))*(1+5*(1-alpha)*T/2/Dstar)*(T/D)**(1-n)*sqrt(D/Dstar)*sqrt(12*pi*b*alpha*(1-alpha)*D/T);
    L2 = sqrt((1-alpha))/pi**1.5*(1+5*(1-alpha)*T/2./Dstar)*(T/D)**(1-n)*sqrt(D/Dstar)*sqrt(12*pi*b*alpha*(1-alpha)*D/T);
    k21 = L2*A*T**n*exp(-D/Ta+deltaD*(1/Ta-1/T));
    k2h = (1-L2)*A*T**n*fac*exp(-D/Tv);#+deltaD*(1/Tv-1/T));
    k2 = k21+k2h;
    E2l = alpha*Dstar*(Tv/Ta)**2
    E2h = D
    E2 = (E2l*k21+E2h*k2h)/k2
    keq = A*T**n*exp(-D/T)
    return (keq,k2,E2,E2/D)

def MFrateDia(T,Tv,A,n,D,m,M,theta):
    b = 2
    omega = 0.5
    alpha = (m/(M+m))**2;
    Ta = alpha*Tv+(1-alpha)*T
    fac = (1-exp(-theta/Tv))/(1-exp(-theta/T))
    deltaD = 3*b*alpha**2*D
    Dstar = D-deltaD
    L2 = 2*(1-alpha)/(pi*pi*alpha**(3/4))*D/Dstar*(T/D)**(2-omega-n)*(1+7*(1-alpha)*(1+sqrt(alpha))*T/2/Dstar);
    k21 = L2*A*T**n*exp(-D/Ta+deltaD*(1/Ta-1/T));
    k2h = (1-L2)*A*T**n*fac*exp(-D/Tv);#+deltaD*(1/Tv-1/T));
    k2 = k21+k2h;
    E2l = alpha*Dstar*(Tv/Ta)**2
    E2h = D
    E2 = (E2l*k21+E2h*k2h)/k2
    keq = A*T**n*exp(-D/T)
    return (keq,k2,E2,E2/D)

def ExtractNeqrate(ArrheniusConstant,sps,spfolder,eqfolder=None):
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
    
    react_sp2 = [[j.replace('<sub>2</sub>','2')   for j in i ]for i in react_sp]


    if not eqfolder:
        eqfolder = os.path.join('..','Eq2',spfolder)

    # obtain subfolders
    rootfolder = next(os.walk(spfolder))[1]
    rootfolder = [i for i in rootfolder if "T=" in i]
    T = [float(i.split('=')[1]) for i in rootfolder]
    rootfolder = [os.path.join(spfolder,i) for i in rootfolder]
    
    # sort
    index = [b[0] for b in sorted(enumerate(T),key=lambda i:i[1])]
    T = [T[i] for i in index]
    rootfolder = [rootfolder[i] for i in index]
    
    # check number of reaction
    with open(os.path.join(eqfolder,'T=%d/DS1REAC.DAT'%(T[0],)),'r') as f1:
        last_line = f1.readlines()[-1]
    nreac = int(float(len(last_line.split())-7)/3.0)
    if nreac !=  len(sps):
        raise Exception("The file doesn't match input for number of reactions")                

    filename = sps[0].replace('+','')+'_NeqRate.dat'


    with open(os.path.join(spfolder,filename),'w') as f:
        # write header of data file
        for i,part in enumerate(react_sp2):
            f.write('# Reac: %d   %s+%s -> %s+%s+%s\n'%(i+1,part[0],part[1],part[2],part[3],part[4]))
                
        f.write('VARIABLES = "Vibrational Temperature (K)","Actual Vibrational Temperature (K)",'+
                '"Rate (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)","Zv",'+
                '"Ev<sub>rem</sub>","Ev<sub>rem</sub>/D" \n')

        for i,tt in enumerate(T):
            Tdir = rootfolder[i]
             # copy dir for equilibrium
            remove_tree(os.path.join(Tdir,'Tv=%d'%(tt,)))
            copy_tree(os.path.join(eqfolder,'T=%d'%(tt,)),os.path.join(Tdir,'Tv=%d'%(tt,)))
                      
            # get all Tv
            Tvdir = next(os.walk(Tdir))[1]
            Tvdir = [i for i in Tvdir if "Tv=" in i]
            Tv = [float(i.split('=')[1]) for i in Tvdir]
            Tvdir = [os.path.join(Tdir,i) for i in Tvdir]

            # sort by temperature
            index = [b[0] for b in sorted(enumerate(Tv),key=lambda i:i[1])]
            Tv = [Tv[i] for i in index]
            Tvdir = [Tvdir[i] for i in index]

            # copy dir for equilibrium
            ieq = Tv.index(tt)

            with open(os.path.join(Tvdir[ieq],'DS1REAC.DAT'),'r') as f1:
                last_line = f1.readlines()[-1]
            
            last_line = last_line.split()

            # get equilibrium rates
            eqrate = [float(last_line[7+nreac+i]) for i in range(nreac)]
            
            for ireac in range(nreac):
                # ====================================write DSMC data ============================
                f.write('ZONE I = %d, T = "DSMC: %.2f Reac: %d"\n'%(len(Tvdir),tt,ireac+1))
                data = np.zeros([len(Tv),6])
                for j,ctvdir in enumerate(Tvdir):
                    # get rates
                    with open(os.path.join(ctvdir,'DS1REAC.DAT'),'r') as f1:
                        last_line = f1.readlines()[-1]
                    last_line = last_line.split()
                    last_line = [float(i) for i in last_line[1:]]   
                    
                    rate = last_line[6+nreac+ireac]
                    ATv = last_line[2]   #actual vibrational temperature
                    ctv = Tv[j]          #simulated vT

                    # get vibrational energy removal            
                    if ireac == 0:
                        fname = 'REAC_COUNT.DAT'
                    else:
                        fname = 'REAC_COUNT_%03d.DAT'%(ireac+1,)
                    with open(os.path.join(ctvdir,fname),'r') as f1:
                        content = f1.readlines()[-1]
                    content = content.split()
                    evrem = float(content[3])  # in ev
                    evremfrac = float(content[2]) 
                    
                    if eqrate[ireac] == 0:
                        Zv = 0
                    else:
                        Zv = rate / eqrate[ireac]

                    data[j] = np.array([ctv, ATv, rate, Zv, evrem, evremfrac])

                for j in data:
                    f.write('{:>9.2f} {:>9.2f}'.format(j[0],j[1]))
                    for i in j[2:]:
                        f.write('{:>13.5G} '.format(i))
                    f.write('\n')
                f.write('\n')
                
            
                # ====================================write CFD-MF data ============================
    
                Tvrange = np.arange(1000.0,20500.0,200.0)
                f.write('ZONE I = %d, T = "CFD: %.2f Reac: %d"\n'%(len(Tvrange),tt,ireac+1))

             
                (A,eta,D,m,M,thetaV) = ArRate[ireac][:-2]

                for ctv in Tvrange:
                    if fdia[ireac]:
                        MF = MFrateDia(tt,ctv,A,eta,D,m,M,thetaV)
                    else:
                        MF = MFrate(tt,ctv,A,eta,D,m,M,thetaV)
                    
                    f.write('{:>9.2f} {:>9.2f}'.format(ctv,ctv))
                    f.write('{:>13.5G} '.format(MF[1]))
                    f.write('{:>13.5G} '.format(MF[1]/MF[0]))
                    f.write('{:>13.5G} '.format(MF[2]*Boltzmann/1.6021766208e-19))
                    f.write('{:>13.5G}\n'.format(MF[3]))
                f.write('\n')

   # copy2(os.path.join(spfolder,sp+'_NeqRate.dat'), os.path.join(sp+'_NeqRate.dat'))


if __name__ == '__main__':
    os.chdir('E://Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/NEq/')
    ExtractNeqrate(ArrheniusConstant,sps=['N2+N',],spfolder='N3',eqfolder='..\\Eq2\\N3')
    ExtractNeqrate(ArrheniusConstant,sps=['N2+O'],spfolder='N2O',eqfolder='..\\Eq2\\N2O')
    ExtractNeqrate(ArrheniusConstant,sps=['O2+O'],spfolder='O3_wysong',eqfolder='..\\Eq2\\O3_wysong')