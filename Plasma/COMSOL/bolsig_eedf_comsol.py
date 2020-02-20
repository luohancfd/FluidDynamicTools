#!/usr/bin/env python3
# %%
import re
import numpy as np


def line2data(line):
    lineSplit = line.split()
    dat = float(lineSplit[-1])

    line2 = ' '.join(lineSplit[:-1])
    if '(' in line2:
        varname, unit = line2.split('(')
        varname = varname.strip()
        unit = unit[:-1].strip()
    else:
        varname = line2
        unit = ''
    return {varname: [dat, unit]}


filename = './BOLSIG/Phelps_1Torr_EEDF.out'

content = []
with open(filename, 'r') as fid:
    for i in fid:
        l = i.strip()
        if l:
            content.append(l)

recomp = [re.compile(i) for i in ['^R(\d+)$', '^\-+$', '[+-]\d\.\d+E[+-]\d+']]

# %% Split file by run ID
cases = []
i = 0
while i < len(content):
    if recomp[0].match(content[i]):
        runID = int(recomp[0].findall(content[i])[0])
        parts = []
        i += 1
        while not recomp[0].match(content[i]):
            if recomp[1].match(content[i]):
                c = []
                i += 1
                while not recomp[1].match(content[i]) and not recomp[0].match(content[i]):
                    c.append(content[i])
                    i += 1
                    if i >= len(content):
                        break
                parts.append(c)
                if i >= len(content):
                    break
        cases.append({'input': parts[0],
                      'output': parts[1],
                      'rate': parts[2],
                      'EEDF': parts[3]})
    else:
        i += 1

# %% Convert to true dict
for i in range(len(cases)):
    for key in ['input', 'output']:
        cases[i][key] = [line2data(j) for j in cases[i][key]]
        cases[i][key] = {k: v for j in cases[i][key] for k, v in j.items()}


# %% Convert EEDF to numpy array
for i in range(len(cases)):
    try:
        EEDF = np.array([[float(kk) for kk in k.split()] for k in cases[i]['EEDF'][1:]])
        cases[i]['EEDF'] = EEDF.T
    except:
        print(i)

EEDF = [i['EEDF'] for i in cases]

# %% Export EEDF
import py2tec
tdata = {'varnames': ['Energy (eV)', 'f(E) (eV<sup>-3/2</sup>)', 'E/e'],
         'lines': [{'data': np.vstack((c['EEDF'][0:2],c['EEDF'][0]/c['output']['Mean energy'][0])), 'zonename': 'Te={:f}'.format(c['output']['Mean energy'][0]/3*2)}   for c in cases]
         }
py2tec.py2tec(tdata, 'BOLSIG/EEDF.tec')

#%%
# with open('E_N_range.txt','w') as fid:
#     for i in range(len(cases)):
#         fid.write('{:f} [Td] '.format(cases[i]['input']['Electric field / N'][0]))

#%% plot E/V vs mean energy
import matplotlib.pyplot as plt
Erd = [i['input']['Electric field / N'][0] for i in cases]
meanEnergy = [i['output']['Mean energy'][0] for i in cases]

fig, ax = plt.subplots()
ax.plot(Erd,meanEnergy)

#%% Check how good is the nomarlization
EEDFint = np.zeros(len(cases))
for i,c in enumerate(cases):
    EEDFint[i] = np.trapz(c['EEDF'][1]*np.sqrt(c['EEDF'][0]), c['EEDF'][0])
plt.plot(EEDFint)

#%%
from numpy.polynomial import polynomial as poly
f_min = 1e-10
for i, c in enumerate(cases):

    # ----- Method 1 --------------
    # Extend EEDF to a point with 0 to strictly satisfy normalization
    # The method fails as new point is ~3E6 eV
    # if EEDFint[i] < 1.0:
    #     x_new = (1.0 - EEDFint[i]) * 2 / (y_end*np.sqrt(x_end)) + x_end
    #
    # ------ Method 2 --------------
    # Extend EEDF to f(e) = 1e-30
    cases[i]['EEDF_ext'] = c['EEDF'][0:2]

    nsample = 5
    coeff = poly.polyfit(c['EEDF'][0][-nsample:], np.log(c['EEDF'][1][-nsample:]),1)

    if c['EEDF'][1][-1] > f_min:
        x_end = (np.log(f_min) - np.log(c['EEDF'][1][-1])) / coeff[1] + c['EEDF'][0][-1]
        y_end = np.exp(np.log(c['EEDF'][1][-1]) + coeff[1] * (x_end - c['EEDF'][0][-1]))
        print(y_end)
        cases[i]['EEDF_ext'] = np.append(cases[i]['EEDF_ext'], np.array([[x_end], [y_end]]), axis=1)

EEDFint_ext = np.zeros(len(cases))
for i,c in enumerate(cases):
    EEDFint_ext[i] = np.trapz(c['EEDF_ext'][1]*np.sqrt(c['EEDF_ext'][0]), c['EEDF_ext'][0])

fig, ax = plt.subplots()
ax.plot(EEDFint_ext,'b')
#ax.plot(EEDFint, 'k')



#%% Export the EEDF to comsol format
with open('Comsol_Input/comsol_eedf.txt', 'w', encoding='utf-8') as fid:
    fid.write('% {:s}  {:s} {:s}\n'.format('Energy (eV)', 'Mean Energy (eV)', 'EEDF (eV^(-3/2))'))
    for c in cases:
        ebar = c['output']['Mean energy'][0]
        for x,y in zip(c['EEDF_ext'][0], c['EEDF_ext'][1]):
            fid.write('{:>12.4e} {:>12.4e} {:>12.4e}\n'.format(x,ebar,y))

#%% What hell is wrong with EEDF for low energy
index = [i for i, x in enumerate(meanEnergy) if x/3*2 < 0.66]
fig, ax = plt.subplots()
j = 150
ax.semilogy(cases[j]['EEDF'][0], cases[j]['EEDF'][1])


#%% understand how bolsig integrates the cross sections
from bolos import parser
from scipy import constants as cs
fname = 'Phelps_1Torr.txt'
with open(fname, 'r') as fp:
    lxcat = parser.parse(fp)
from bolos.process import Process
lxcat_processed = [Process(**i) for i in lxcat]
# N2 ionization
xsec = lxcat_processed[24].data.T
rate = np.zeros(len(cases))
for i, c in enumerate(cases):
    x = c['EEDF'][0]
    dist = c['EEDF'][1]
    x_xsec = np.interp(x, xsec[0],xsec[1])
    y = np.sqrt(cs.electron_volt*2/cs.electron_mass)*x*x_xsec*dist
    rate[i] = np.trapz(y,x)

fig, ax = plt.subplots()
ax.loglog([i/3*2 for i in meanEnergy], rate)


