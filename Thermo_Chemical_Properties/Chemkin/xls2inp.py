#!/usr/bin/python3
# -*- coding: utf-8 -*-
import csv
import argparse
import platform
from openpyxl import load_workbook

if platform.system() == 'Windows':
    CHARSET = 'utf-8-sig'
else:
    CHARSET = 'utf-8'

AUXOPTION = ["CHEB", "EXCI", "FIT1", "FORD", "HIGH",
             "JAN", "LOW", "LT", "MOME", "PCHEB", "PLOG",
             "REV", "RLT", "RORD", "SRI", "TCHEB", "TDEP",
             "TROE", "UNITS", "USRPROG", "XSMI"]
numForm = '{:<30s} {:>10.3E} {:>8.3f} {:>10.1f}\n'
strForm = '{:<30s} {:>10s} {:>8s} {:>10s}\n'

def Xls2Inp(ARGS, csvParams=None):
    with open(ARGS.OUTPUT, 'w', encoding='utf-8') as f:
        f.write('REACTIONS\n')
        if '.xls' in ARGS.INPUT or '.xlsx' in ARGS.INPUT:
            fr = load_workbook(filename=ARGS.INPUT, read_only=True)
            rowIter = fr.active.rows
            for row in rowIter:
                if row[0].value in AUXOPTION:
                    f.write('   {:s} / '.format(row[0].value))
                    for i in row[1:-1]:
                        f.write('{:>10.3E} '.format(i.value))
                    f.write('{:>10.3E} /\n'.format(row[-1].value))
                elif row[0].value == 'COMMENT':
                    for i in row[1:]:
                        if i.value:
                            f.write('! {}\n'.format(i.value))
                else:
                    data = [i.value for i in row]
                    f.write(numForm.format(*data))
        else:
            fr = open(ARGS.INPUT, 'r', encoding=CHARSET)
            if csvParams:
                rowIter = csv.reader(fr, **csvParams)
            else:
                rowIter = csv.reader(fr)
            for row in rowIter:
                if row[0] in AUXOPTION:
                    f.write('   {:s} / '.format(row[0]))
                    for i in row[1:-1]:
                        f.write('{:>10.3E} '.format(float(i)))
                    f.write('{:>10.3E} /\n'.format(float(row[-1])))
                elif row[0] == 'COMMENT':
                    for i in row[1:]:
                        if i:
                            f.write('! {}\n'.format(i))
                else:
                    data = [i for i in row]
                    f.write(strForm.format(*data))
        fr.close()
        f.write('END')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert csv to chem.inp reaction bloack')

    parser.add_argument('INPUT', help='input file in csv format')
    parser.add_argument('OUTPUT', help='output file')

    ARGS = parser.parse_args()
    Xls2Inp(ARGS)
