# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 11:24:27 2018

@author: Han
"""
from types import MethodType
import struct

def writeb (self,a,atype='i'):
    self.write(struct.pack(atype,a))

def writefloat32(self,a,flen=64):
    self.write(struct.pack('f',a))
    
def writestring (self,a):
    for i in a:
        self.write(struct.pack('I',ord(i)))
    self.write(struct.pack('I',0))

def writechar (self,a):
    for i in a:
        self.write(struct.pack('B',ord(i)))
    
if ( __name__ == '__main__' ):
    with open('test.plt','wb') as fh:
        # add method
        fh.writeb = MethodType(writeb, fh)
        fh.writestring = MethodType(writestring, fh)
        fh.writechar = MethodType(writechar, fh)

        # magic header
        fh.writechar('#!TDV112')
        # write one
        fh.writeb(1) 
        # full 
        fh.writeb(0)
        # title
        fh.writestring('abc')
        # number of variables
        var=['x','y','v']
        x= [1,2,3,1,2,3,4,4]
        y= [0,0,0,1,1,1,0,1]
        v = [1,3,2]
        conn=[[1,2,5,4],[2,3,6,5],[3,7,8,6]]
        NumPts = len(x)
        NumElements = len(v)
        Nvar = len(var)
        fh.writeb(Nvar) 
        for i in var:
            fh.writestring(i)
            
        #============= ZONE
        fh.writeb(299.0,'f')
        # zone name
        fh.writestring('ZONE1')
        # no parent zone
        fh.writeb(-1)
        # no strandID
        fh.writeb(-2)
        # solution time
        fh.writeb(0.0,'d')
        # not used
        fh.writeb(-1)
        # zone type
        fh.writeb(3)
        # specify var location
        fh.writeb(1) 
        # var location
        fh.writeb(0)
        fh.writeb(0)
        fh.writeb(1)
        #raw local 1-to-1 face neighbors not supplied
        fh.writeb(0)
        #No. of miscellaneous user defined face neighbor connections
        NoOfUserDefinedNeighbourConn = 0; # 0 for no, >0 for yes
        fh.writeb(NoOfUserDefinedNeighbourConn)
        # NumPts
        fh.writeb(NumPts)
        # NumElements
        fh.writeb(NumElements)
        #Reserved by tecplot for future use; set to zero here
        ICellDim=0
        JCellDim=0
        KCellDim=0
        fh.writeb(ICellDim)
        fh.writeb(JCellDim)
        fh.writeb(KCellDim)
        # no aux        
        aux_data=0
        fh.writeb(aux_data)
        
        #This EOH_MARKER separates HEADER SECTION and DATA SECTION in tecplot
        #binarly data format
        #here EOH means '\eOH' beginning-of-line bindkey
        #here EOF means '\eOF' end-of-line bindkey 
        EOH_MARKER = 357.0
        fh.writeb(EOH_MARKER,'f')
        # FiniteElement zone
        fh.writeb(299.0,'f')
        #var_prec  is calculated above  using tdata.vformat
        # i.e. its value should be 1 for float
        #                          2 for double
        #                          3 for longInt
        #                          4 for shortInt
        #                          5 for Byte
        #                          6 for Bit
        for i in range(Nvar):
            fh.writeb(1)
        # No passive variable
        has_passive = 0
        fh.writeb(has_passive)
        has_var_share = 0
        fh.writeb(has_var_share)
        zone_number_to_share_connectivity=-1
        fh.writeb(zone_number_to_share_connectivity)
        
        # find min and max pair of each variable and output them
        fh.writeb(min(x),'d')
        fh.writeb(max(x),'d')
        fh.writeb(min(y),'d')
        fh.writeb(max(y),'d')            
        fh.writeb(min(v),'d')
        fh.writeb(max(v),'d')
        
        # write zone
        for i in x+y+v:
            fh.writeb(i,'f')
        # Zone Connectivity  minus 1 to get to 0-based 
        for i in conn:
            for j in i:
                fh.writeb(j-1)
        