#!/usr/bin/env python3
import sys,re,os
f90file = sys.argv[1]

with open(f90file, 'r') as f:
    c = f.readlines()

fname = os.path.basename(f90file)
fname,fext = fname.split('.')

if not os.path.isdir(fname):
    os.mkdir(fname)

patternGroup = [['program', 'end program'],
                ['module', 'end module'],
                ['subroutine', 'end subroutine'],
                ['function', 'end function'],
                ['real function', 'end function'],
                ['double function', 'end function'],
                ['integer function', 'end function']]

pCs = [ dict() for i in patternGroup]
pRes = [[re.compile('^\s*{:}\s+(\w+)\s*(?:\s!.*)?'.format(i),re.IGNORECASE)  for i in j] for j in patternGroup]

inModule = False
i = 0

buffer = []
while i < len(c):
    for j, pRe in enumerate(pRes):
        pMstart = pRe[0].findall(c[i])
        if pMstart:
            pName = pMstart[0]
            print('Line {:d}  {:s} {:s}'.format(i+1, patternGroup[j][0].upper(), pName))
            buffer.append(c[i])
            i += 1
            while i < len(c):
                buffer.append(c[i])
                pMend = pRe[1].findall(c[i])
                if pMend:
                    print('Line {:d}  END {:s} {:s}'.format(i, patternGroup[j][0].upper(), pName))
                    i += 1
                    break

                pMend = re.findall('^\s*end\s*(?:\s!.*)?\n$', c[i], re.IGNORECASE)
                if pMend:
                    print('Line {:d}  END {:s}'.format(i, pName))
                    i += 1
                    break

                i += 1

            if pName in pCs[j]:
                raise ValueError('Duplicated {:s} {:s}'.format(patternGroup[j][0].upper(), pName))
            pCs[j][pName] = buffer
            buffer = []
            break
    if not pMstart:
        buffer.append(c[i])
        i += 1

for i in pCs:
    for j in i:
        with open('{:s}/{:s}.{:s}'.format(fname,j,fext), 'w', encoding='utf8') as f:
            for k in i[j]:
                f.write(k)
