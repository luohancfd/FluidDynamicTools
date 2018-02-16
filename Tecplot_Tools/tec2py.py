#!/usr/bin/env python

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