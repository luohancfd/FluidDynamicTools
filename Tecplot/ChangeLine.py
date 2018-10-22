# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 23:20:52 2017

@author: Han
"""

#import logging
#logging.basicConfig(level=logging.DEBUG)
import tecplot as tp
import glob,re,os
from tecplot.constant import PlotType, Color, LinePattern

def save_eps(filename):
    command = r'''#!MC 1410
               $!PRINTSETUP PALETTE = COLOR
               $!EXPORTSETUP EXPORTFNAME = "%s"
               $!EXPORTSETUP EXPORTFORMAT = EPS
               $!EXPORT  EXPORTREGION = CURRENTFRAME
              '''%(filename,)
   # print(command)
    tp.macro.execute_command(command)



cwd =  os.getcwd().split('/')[-1]
print('CurrentFolder : %s'%(cwd,))
infiles = glob.glob('*/*T*.lay',recursive=True)
infiles = [i for i in infiles if 'VibDistT=10000.lay' not in i]
infiles = [os.path.join(cwd,i) for i in infiles]
Temperature = [re.search('T[=_]?(\d+)',i).group(1)  for i in infiles]
outfiles = [os.path.join(os.path.dirname(i),'T='+j+'.eps') for i,j in zip(infiles,Temperature)]

for infile,outfile in zip(infiles,outfiles):
    tp.load_layout(infile)
    frame = tp.active_frame()
    plt = frame.plot()
    for lmapid in plt.active_linemap_indices:
        lmap = plt.linemap(lmapid)
        zonename = lmap.zone.name
        line = lmap.line
        if 'ROT' in zonename or 'rot' in zonename:
            print(zonename)
            line.color = Color.Green
            line.line_pattern = LinePattern.LongDash
            line.pattern_length = 0.8
        elif 'MF' in zonename:
            print(zonename+'  2')
            line.color = Color.Blue
            line.line_pattern = LinePattern.Dashed
            line.pattern_length = 2.0
    save_eps(outfile)

    epsdir = os.path.split(os.path.dirname(outfile))[-1]
    fname = os.path.basename(outfile)
    filename2 = re.sub('[\+=]','',epsdir+fname)
    filename2 = os.path.join(cwd,'Figs',filename2)
    save_eps(filename2)

    print('%s is exported'%(outfile,))
    tp.save_layout(infile)

