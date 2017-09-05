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
from distutils.dir_util import copy_tree


eV = scipy.constants.physical_constants['electron volt'][0]

ArrheniusConstant = {'N2O':[8.9340e-7*N_A,  -0.4807,   113950, 14, 16, 3352],
                     'N3': [6.6831e-6*N_A,  -0.6996,   113200, 14, 14, 3352],
                     'O3': [      2.5e18,   -0.565,60491.819, 16, 16, 2273.54]}

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

def ExtractEqrate(ArrheniusConstant,sp,spfolder=None):
    A = ArrheniusConstant[sp][0]
    eta = ArrheniusConstant[sp][1]
    D = ArrheniusConstant[sp][2]
    m = ArrheniusConstant[sp][3]
    M = ArrheniusConstant[sp][4]
    thetaV = ArrheniusConstant[sp][5]
    if not spfolder:
        spfolder = sp

    rootfolder = next(os.walk(spfolder))[1]

    header = ["Tt", "Tr", "Tv", "Tave", "P", "Number of reactions"]

    datafile = os.path.join(spfolder,sp+'_EqRate.dat')
    with open(datafile,'w') as f:
        with open(os.path.join(spfolder,rootfolder[0],'DS1REAC.DAT'),'r') as f1:
            last_line = f1.readlines()[-1]
        nreac = int(float(len(last_line.split())-7)/2.0)
        ievrem = int(os.path.isfile(os.path.join(spfolder,rootfolder[0],'REAC_COUNT.DAT')))


        f.write('VARIABLES = ')
        for i in header:
            f.write('"%s",'%(i,))
        for i in range(nreac):
            f.write('"Frac%d", '%(i+1,))
        for i in range(nreac):
            f.write('"Rate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)",'%(i+1,))
        for i in range(nreac-1):
            f.write('"ColRate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)",'%(i+1,))
        f.write('"ColRate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)"'%(nreac))
        f.write(',"Ev<sub>rem</sub>", "Ev<sub>rem</sub>/D" \n')

        f.write('ZONE I = %d, T = "%s" \n'%(len(rootfolder),'DSMC'))

        data = np.zeros([len(rootfolder),3*nreac+6+2])

        for j,i in enumerate(rootfolder):
            with open(os.path.join(spfolder,i,'DS1REAC.DAT'),'r') as f1:
                content = f1.readlines()[-1]
            content = content.split()
            content = content[1:]
            data[j][:len(content)] = np.array([float(k) for k in content])

            if ievrem == 1:
                with open(os.path.join(spfolder,i,'REAC_COUNT.DAT'),'r') as f1:
                    content = f1.readlines()[-1]
                content = [float(i) for i in content.split()]
                data[j][-2:] = np.array([content[1], content[2]])

        data = data[data[:,0].argsort()]

        for j in data:
            f.write('{:>9.2f} '.format(j[0]))
            for i in j[1:]:
                f.write('{:>13.5G} '.format(i))
            f.write('\n')


        f.write("\n")
        Trange = np.arange(1000.0,20500.0,500.0)
        f.write('ZONE I = %d, T = "%s" \n'%(len(Trange),'QCT'))
        f.write('PASSIVEVARLIST=[2-%d,'%(len(header),))
        if nreac == 1:
            f.write('%d,'%(len(header)+1,))
        else:
            f.write('%d-%d,'%(len(header)+1,len(header)+nreac))
        f.write('%d]\n'%(len(header)+nreac*2+1,))

        for T in Trange:
            a = MFrate(T,T,A,eta,D,m,M,thetaV)
            f.write('{:>9.2f} {:>13.5G} {:>13.5G} {:>13.5G} \n'.format(T,a[0],a[2],a[3] ))


#    if os.path.isfile('rate.lay'):
#        copy2(os.path.join('.','rate.lay'),os.path.join('.',sp+'_Rate.lay'))
#        TecplotSub(os.path.join('.',sp+'_Rate.lay'),sp+'_EqRate.dat')
#    if os.path.isfile('Evrem.lay'):
#        copy2(os.path.join('.','Evrem.lay'),os.path.join('.',sp+'_Evrem.lay'))
#        TecplotSub(os.path.join('.',sp+'_Evrem.lay'),sp+'_EqRate.dat')
    copy2(datafile,os.path.join('.',sp+'_EqRate.dat'))

def ExtractNeqrate(ArrheniusConstant,sp,spfolder=None,eqfolder=None):
    A = ArrheniusConstant[sp][0]
    eta = ArrheniusConstant[sp][1]
    D = ArrheniusConstant[sp][2]
    m = ArrheniusConstant[sp][3]
    M = ArrheniusConstant[sp][4]
    thetaV = ArrheniusConstant[sp][5]
    if not spfolder:
        spfolder = sp

    if not eqfolder:
        eqfolder = os.path.join('..','Eq2',spfolder)

    rootfolder = next(os.walk(spfolder))[1]
    T = [float(i.split('=')[1]) for i in rootfolder]

    rootfolder = [os.path.join(spfolder,i) for i in rootfolder]

    with open(os.path.join(spfolder,sp+'_NeqRate.dat'),'w') as f:
        f.write('VARIABLES = "Vibrational Temperature (K)","Actual Vibrational Temperature (K)",'+
                '"Rate (cm<sup>2</sup>mol<sup>-1</su>)","Zv",'+
                '"Ev<sub>rem</sub>","Ev<sub>rem</sub>/D" \n')


        for i,tt in enumerate(T):
            Tdir = rootfolder[i]
            Tvdir = next(os.walk(Tdir))[1]
            Tv = [float(i.split('=')[1]) for i in Tvdir]
            Tvdir = [os.path.join(Tdir,i) for i in Tvdir]

            # sort by temperature
            index = [b[0] for b in sorted(enumerate(Tv),key=lambda i:i[1])]
            Tv = [Tv[i] for i in index]
            Tvdir = [Tvdir[i] for i in index]

            f.write('ZONE I = %d, T = "DSMC: %.2f"\n'%(len(Tvdir),tt))

            # copy dir for eq
            ieq = Tv.index(tt)
            copy_tree(os.path.join(eqfolder,'T=%d'%(tt,)),Tvdir[ieq])

            # get equilibrium rate
            with open(os.path.join(Tvdir[ieq],'DS1REAC.DAT'),'r') as f1:
                last_line = f1.readlines()[-1]
            eqrate = float(last_line.split()[8])

            data = np.zeros([len(Tvdir),6])
            for j,ctvdir in enumerate(Tvdir):
                with open(os.path.join(ctvdir,'DS1REAC.DAT'),'r') as f1:
                    last_line = f1.readlines()[-1]
                rate = float(last_line.split()[8])
                ATv = float(last_line.split()[3])
                ctv = Tv[j]

                with open(os.path.join(ctvdir,'REAC_COUNT.DAT'),'r') as f1:
                    last_line = f1.readlines()[-1]
                evrem = [float(ii) for ii in last_line.split()[1:3]]

                if eqrate == 0:
                    Zv = 0
                else:
                    Zv = rate / eqrate

                data[j] = np.array([ctv, ATv, rate, Zv, evrem[0], evrem[1]])

            for j in data:
                f.write('{:>9.2f} {:>9.2f}'.format(j[0],j[1]))
                for i in j[2:]:
                    f.write('{:>13.5G} '.format(i))
                f.write('\n')
            f.write('\n')

            Tvrange = np.arange(1000.0,20500.0,200.0)
            f.write('ZONE I = %d, T = "CFD: %.2f"\n'%(len(Tvrange),tt))
            for ctv in Tvrange:
                MF = MFrate(tt,ctv,A,eta,D,m,M,thetaV)
                f.write('{:>9.2f} {:>9.2f}'.format(ctv,ctv))
                f.write('{:>13.5G} '.format(MF[1]))
                f.write('{:>13.5G} '.format(MF[1]/MF[0]))
                f.write('{:>13.5G} '.format(MF[2]))
                f.write('{:>13.5G}\n'.format(MF[3]))
            f.write('\n')

    copy2(os.path.join(spfolder,sp+'_NeqRate.dat'), os.path.join(sp+'_NeqRate.dat'))


if __name__ == '__main__':
    ExtractNeqrate(ArrheniusConstant,sp='N2O',spfolder='N2O',eqfolder=None)
    ExtractNeqrate(ArrheniusConstant,sp='O3',spfolder='O3',eqfolder=None)
    ExtractNeqrate(ArrheniusConstant,sp='N3',spfolder='N3',eqfolder=None)
