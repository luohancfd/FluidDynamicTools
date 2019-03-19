#!/usr/bin/env python
import os
import shutil
from glob import glob
import argparse


parser = argparse.ArgumentParser(description='Reader for DS1VD output')
parser.add_argument('folder', nargs='+', help='folder of DS1FP files')
ARGS=parser.parse_args()

files = []
for f in ARGS.folder:
    files += glob('{:s}\\DS1FP00_*.DAT'.format(f))
files = ['"{:s}"'.format(i) for i in files]

tecFile = 'ShockMovie.lay'

with open(tecFile, 'r') as f:
    data = f.readlines()

with open(tecFile, 'w', encoding='utf-8') as f:
    for i in data:
        if '$!VarSet'  in i and 'LFDSFN1' in i:
            f.write('$!VarSet |LFDSFN1| = \'{:s}\'\n'.format(' '.join(files)))
        else:
            f.write('{:s}'.format(i))