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
import re
from scipy.interpolate import interp1d
#import matlab.engine

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
            #         'O4': [ 1.9337e22, -1.7334, 5.21275*eV2K, 16 , 16, 2273.54, 'O2', 'O2']}
science_not = re.compile('(.*)([+-]{1}\d+)')

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

#def EqRateFromMatlab(T,Tv,sp1,sp2,D,eng):
#    a = eng.MatlabMFMCrate(T,Tv,sp1,sp2,1e6,D,2,nargout=4)
#    return a[0:3]

def myfloat(a):
    try:
        b = float(a)
    except ValueError:
        w = science_not.findall(a)[0]
        b = float(w[0])*10**float(w[1])
    except:
        b = 0
        raise
    return b

def ReadTecplotDat(datfile):
    IVARlist = False
    data = []
    
    # read data once first
    with open(datfile,'r') as f1:
        flines = f1.readlines()
        

    
    for iline, line in enumerate(flines):
        l = line.strip()
        if not l:
            continue  #jump out if this line is empty
        if l[0] != '#':
            if 'VARIABLES' in l:
                VARlist = re.findall('[''"](.*?)[''"]',l)
                NVAR = len(VARlist)
                IVARlist = True
                break
            
    flines = iter(flines[iline+1:])
       
    if not IVARlist:
        raise Exception("No variables found")
    else:
        nzone = 0
        while True:
            try:
                l = next(flines).strip()  #f1.readline().strip()
                if not l: continue
                if l[0] == '#': continue
                if 'ZONE' in l:
                    nzone = nzone + 1
                    # ===========load zone information===================
                    inum = re.findall('I\s*=\s*(\d+)',l)
                    inum = int(inum[0])
                    jnum = re.findall('I\s*=\s*(\d+)',l)
                    if not jnum:
                        jnum = int(jnum[0])
                    else:
                        jnum = 1
                    ndata = inum*jnum
                    title = re.findall('T\s*=\s*["''](.*)[''"]',l)
                    if not title:
                        title = re.findall('T\s*=\s*(.*)',l)
                    if title:
                        title = title[0]
                    else:
                        title = 'data %d'%(nzone,)
                    
                    # ============== load passivelist ==================
                    mcolum = NVAR;   
                    while True:
                        l = next(flines).strip()   #f1.readline().strip()
                        if not l: continue
                        if l[0] == '#': continue
                        break
                    
                    passive_list = []         
                    i_passive = False
                    if "PASSIVEVARLIST" in l:
                        i_passive = True
                        passive_str = re.findall('\[(.*)\]',l)[0]
                        passive_str = passive_str.split(',')
                        for index,i in enumerate(passive_str):
                            if i.isdigit():
                                passive_list.append(int(i))
                            else:
                                j = [int(k)     for k in i.split('-')]
                                for k in range(j[0],j[1]+1):
                                    passive_list.append(k)
                        mcolum -= len(passive_list)
                    
                    npassive_list = [ i for i in range(NVAR) if i+1 not in passive_list]
                    
                    #  ============== load data ===========================
                    zone_data = []
                    nelement = 0
                    if not i_passive:
                        # we alreday read a line of data
                        l2append = [myfloat(i) for i in l.split()]
                        nelement += len(l2append)
                        zone_data += l2append

                        
                    while nelement < ndata*mcolum:
                        while True:
                            l = next(flines).strip()
                            if not l: continue                        
                            if l[0] == '#': continue
                            break
                        l2append = [myfloat(i) for i in l.split()]
                        nelement += len(l2append)
                        zone_data += l2append
                        
                    zone_data = np.array(zone_data)
                    
                    if jnum == 1:
                        #this is x-y line
                        zone_data = zone_data.reshape((ndata,mcolum))
                    else:
                        # this is 2D data
                        # following tecplot doc p134
                        zone_data = zone_data.reshape((mcolum,jnum,inum))
                        # first index loop NVAR, second jnum, third inum
                    
                    if i_passive:
                        # add zero if there is passive_list
                        if jnum == 1:
                            zone_data2 = np.zeros((ndata,NVAR))
                            for i,j in enumerate(npassive_list):
                                zone_data2[:,j] = zone_data[:,i]
                        else:
                            zone_data2 = np.zeros((mcolum,jnum,inum))
                            for i,j in enumerate(npassive_list):
                                zone_data2[j] = zone_data[i]
                        
                        zone_data = zone_data2
                    
                    
                    data.append({'title':title,'data':zone_data,'passivevarlist':passive_list})
            except StopIteration:
                break
    return (data,VARlist)
                
            
def CalulateEvRem(datfile,T):
    (data,VARlist) = ReadTecplotDat(datfile)
    nreac = len(data)
    Evrem = []
    for ireac in range(nreac):
        dat = data[ireac]['data']        
        evindex = [i    for i,j in enumerate(VARlist) if j == 'Ev (eV)' ][0]

        # filter useless rows
        istop = -1
        for i,j in enumerate(dat[:,evindex]):
            if j == 0:
                istop = i
                break
            elif j>=1e3:
                istop = i
                break
        if istop > 0:
            dat = dat[:istop]
            
        Ev = dat[:,evindex]

            
        c1 =  [i for i,j in enumerate(VARlist) if j == 'Ev_N'][0]
        c2 =  [i for i,j in enumerate(VARlist) if j == 'EvR_N'][0]
           
        (m,n) = dat.shape
        Vprob = np.zeros(m)
        for i,(a,b) in enumerate(zip(dat[:,c1],dat[:,c2])):
            if a<1:
                Vprob[i] = -np.nan
            else:
                Vprob[i] = b/a
        
        # interpolate probability
        index = [i for i,j in enumerate(Vprob) if np.isnan(j)]
        if index:
            index2 = [i for i in range(len(Vprob)) if i not in index]
            fintep = interp1d(Ev[index2],np.log10(Vprob[index2]),fill_value="extrapolate")
            Vprob[index] = 10** fintep(Ev[index])
        
        CEvrem = 0.
        ReacParFun = 0.
        for i in range(len(Ev)):
            CEvrem += Ev[i]*Vprob[i]*np.exp(-Ev[i]*eV/Boltzmann/T)
            ReacParFun += Vprob[i]*np.exp(-Ev[i]*eV/Boltzmann/T)
                
       
        if ReacParFun > 0 :
            CEvrem = CEvrem/ReacParFun
        else:
            CEvrem = 0.0
        
        if np.isnan(CEvrem):
            CEvrem = 0.0
            
        Evrem.append(CEvrem)
       
    return Evrem
            

def ExtractEqrate(ArrheniusConstant,sps,spfolder):
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
                
    header = ["Tt", "Tr", "Tv", "Tave", "P", "Number of reactions"]

   
    filename = sps[0].replace('+','')+'.dat'
    
    rootfolder = next(os.walk(spfolder))[1]
    rootfolder = [i for i in rootfolder if "T=" in i]
    T = [float(i.split('=')[1]) for i in rootfolder]

    header = ["Tt", "Tr", "Tv", "Tave", "P", "Number of reactions"]

    datafile = os.path.join(spfolder,filename)
    with open(datafile,'w') as f:
        with open(os.path.join(spfolder,rootfolder[0],'DS1REAC.DAT'),'r') as f1:
            last_line = f1.readlines()[-1]
        nreac = int(float(len(last_line.split())-7)/3.0)
        if nreac !=  len(sps):
            raise Exception("The file doesn't match input for number of reactions")
        ievrem = int(os.path.isfile(os.path.join(spfolder,rootfolder[0],'REAC_COUNT.DAT')))

        # write header of data file
        for i,part in enumerate(react_sp2):
            f.write('# Reac: %d   %s+%s -> %s+%s+%s\n'%(i+1,part[0],part[1],part[2],part[3],part[4]))
        

            
        f.write('VARIABLES = ')
        for i in header:
            f.write('"%s",'%(i,))
        for i in range(nreac):
            f.write('"Frac%d", '%(i+1,))
        for i in range(nreac):
            f.write('"Rate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)",'%(i+1,))
        for i in range(nreac-1):
            f.write('"ColRate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)",'%(i+1,))
        f.write('"ColRate%d (cm<sup>3</sup>mol<sup>-1</sup>s<sup>-1</sup>)"'%(nreac,))
        if ievrem == 1:
            f.write(',')
            for i in range(nreac):
                f.write('"Ev<sub>rem</sub>%d(eV)", "Ev<sub>rem</sub>%d/D", '%(i+1,i+1,))
            for i in range(nreac-1):
                f.write('"Ev2<sub>rem</sub>%d(eV)", "Ev2<sub>rem</sub>%d/D", '%(i+1,i+1,))
            f.write('"Ev2<sub>rem</sub>%d(eV)", "Ev2<sub>rem</sub>%d/D" \n'%(nreac,nreac,))
        else:
            f.write('\n')    
        

        #=========================write zone of DSMC data==================================
        f.write('\nZONE I = %d, T = "%s" \n'%(len(rootfolder),'DSMC'))
        data = np.zeros([len(rootfolder),len(header)+7*nreac])

        for j,i in enumerate(rootfolder):
            with open(os.path.join(spfolder,i,'DS1REAC.DAT'),'r') as f1:
                content = f1.readlines()[-1]
            content = content.split()
            content = content[1:]  #we don't need time
            data[j][:len(content)] = np.array([float(k) for k in content])
            
            IMF_Ev = os.path.join(spfolder,i,'IMF_EV.DAT')
            Evrem = CalulateEvRem(IMF_Ev,T[j])  #calculate Evrem based on probabilty

            if ievrem == 1:
                for ii in range(1,nreac+1):
                    if ii == 1:
                        fname = 'REAC_COUNT.DAT'
                    else:
                        fname = 'REAC_COUNT_%03d.DAT'%(ii,)
                    with open(os.path.join(spfolder,i,fname),'r') as f1:
                        content = f1.readlines()[-1]
                    content = [float(i) for i in content.split()]
                    
                    kk = len(header)+3*nreac+(ii-1)*2
                    data[j][kk:kk+2] = np.array([content[3], content[2]])
                    
                    kk = len(header)+5*nreac+(ii-1)*2
                    data[j][kk] = Evrem[ii-1]
                    data[j][kk+1] = Evrem[ii-1]*eV/ArRate[ii-1][-1]/Boltzmann

                
        data = data[data[:,0].argsort()]
        for j in data:
            f.write('{:>9.2f} '.format(j[0]))
            for i in j[1:]:
                f.write('{:>13.5G} '.format(i))
            f.write('\n')
        f.write("\n")

        
        
        #=========================write zone of QCT/EQ data==================================
        Trange = np.arange(1000.0,20500.0,500.0)
        
        for i in range(nreac):
            zonename = 'QCT: %s+%s'%(react_sp[i][0],react_sp[i][1])
            f.write('ZONE I = %d, T = "%s" \n'%(len(Trange),zonename))
            f.write('PASSIVEVARLIST=[2-%d,'%(len(header)+nreac,)) # ignore temperature abd fraction
            if nreac != 1:
                active = [i for i in range(len(header)+nreac+1,len(header) + nreac*7+1)]
                number2pop =[]
                number2pop.append(len(header)+nreac+i+1)  # reaction rates
                number2pop.append(len(header)+nreac*3+2*i+1) #evrem
                number2pop.append(len(header)+nreac*3+2*i+2) #evrem/D
               
                passive = [i  for i in active if i not in number2pop]
                
                    
                for j in passive[:-1]:
                    f.write('%d,'%(j,))


                f.write('%d]\n'%(passive[-1],))     
            else:
                f.write('%d,%d,%d]\n'%(len(header)+nreac*2+1,12,13))
                
            (A,eta,D,m,M,thetaV) = ArRate[i][:-3]
            for T in Trange:
                if fdia[i]:
                    a = MFrateDia(T,T,A,eta,D,m,M,thetaV)
                else:
                    a = MFrate(T,T,A,eta,D,m,M,thetaV)
                f.write('{:>9.2f} {:>13.5G} {:>13.5G} {:>13.5G} \n'.format(T,a[0],a[2]*Boltzmann/eV,a[3] ))
            f.write("\n")

        # write zone of Matlab
 #       eng = matlab.engine.start_matlab()
#        Trange = [1000.0,2000.0,3000.0,5000.0,6250.0,7500.0,8750.0,10000.0,11250.0,12500.0,13750.0,15000.0,16250.0,17500.0,18750.0,20000.0]
#        f.write('ZONE I = %d, T = "%s" \n'%(len(Trange),'Matlab'))
#        f.write('PASSIVEVARLIST=[2-%d,'%(len(header),))
#        if nreac == 1:
#            f.write('%d,'%(len(header)+1,))
#        else:
#            f.write('%d-%d,'%(len(header)+1,len(header)+nreac))
#        f.write('%d]\n'%(len(header)+nreac*2+1,))
#        for T in Trange:
#            a = EqRateFromMatlab(T,T,sp1,sp2,D,eng)
#            f.write('{:>9.2f} {:>13.5G} {:>13.5G} {:>13.5G} \n'.format(T,a[0],a[1]*eV2K,a[2] ))
#        eng.quit()
#    if os.path.isfile('rate.lay'):
#        copy2(os.path.join('.','rate.lay'),os.path.join('.',sp+'_Rate.lay'))
#        TecplotSub(os.path.join('.',sp+'_Rate.lay'),sp+'_EqRate.dat')
#    if os.path.isfile('Evrem.lay'):
#        copy2(os.path.join('.','Evrem.lay'),os.path.join('.',sp+'_Evrem.lay'))
#        TecplotSub(os.path.join('.',sp+'_Evrem.lay'),sp+'_EqRate.dat')
   # copy2(datafile,os.path.join('.',sp+'_EqRate.dat'))



if __name__ == "__main__":
     os.chdir('E://Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Eq2/')
     ExtractEqrate(ArrheniusConstant,sps=['N2+O'],spfolder='N2O')
     ExtractEqrate(ArrheniusConstant,sps=['O2+O'],spfolder='O3')
     ExtractEqrate(ArrheniusConstant,sps=['N2+N'],spfolder='N3')

     os.chdir('E://Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/Eq/')
     ExtractEqrate(ArrheniusConstant,sps=['O2+N2','N2+O2'],spfolder='O2N2Double')
     ExtractEqrate(ArrheniusConstant,sps=['NO+N2'],spfolder='NON2')
     ExtractEqrate(ArrheniusConstant,sps=['O2+O2'],spfolder='O4_Byron')
     ExtractEqrate(ArrheniusConstant,sps=['N2+N2'],spfolder='N4_2')
