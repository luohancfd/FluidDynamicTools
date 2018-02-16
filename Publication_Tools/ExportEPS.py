#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG)
import tecplot
import glob,re,os,shutil

def save_eps(filename):
    command = r'''#!MC 1410
               $!PRINTSETUP PALETTE = COLOR
               $!EXPORTSETUP EXPORTFNAME = "%s"
               $!EXPORTSETUP EXPORTFORMAT = EPS
               $!EXPORT  EXPORTREGION = CURRENTFRAME
              '''%(filename,)
   # print(command)
    tecplot.macro.execute_command(command)

## Export equilibrium reaction rates

infiles = ["N2N2.lay",]

name = [os.path.basename(i) for i in infiles]
outfiles = [os.path.join("/home/Figure",i.replace(".lay","_EqRate.eps")) 
            for i in name]
outfiles2 = [os.path.join("/home/Draft/Figs",i.replace(".lay","_EqRate.eps")) 
            for i in name]

for i,j in enumerate(infiles):
    tecplot.load_layout(j)
    save_eps(outfiles[i])
    shutil.copyfile(outfiles[i],outfiles2[i])
 