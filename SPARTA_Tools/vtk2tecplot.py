import vtk,struct
from types import MethodType
# Changelog:
# - March 2018: Modified by luohancfd to enable export binary tecplot file

# define a contex manager
class BinaryFile():
    def __init__(self,filename,mode):
        self.filename = filename
        self.mode = mode
    def __enter__(self): #this is given to with
        self.fid = open(self.filename,self.mode)
        return self
    def __exit__(self,*args):
        self.fid.close()
    def writeb (self,a,atype='i'):
        self.write(struct.pack(atype,a))
    def writestring (self,a):
        for i in a:
            self.fid.write(struct.pack('I',ord(i)))
        self.fid.write(struct.pack('I',0))
    def writechar (self,a):
        for i in a:
            self.fid.write(struct.pack('B',ord(i)))
    def write(self,*args):
        self.fid.write(*args)
        
def write_points(points,
                 fh,
                 index):
    for idx, i in enumerate(range(0, points.GetNumberOfPoints())):
        x = points.GetPoint(i)
        fh.write(str(x[index]))
        if idx % 5 == 0 and idx != 0:
            fh.write("\n")
        else:
            fh.write(" ")
    fh.write("\n")

def write_array(array,
                fh):
    for idx, i in enumerate(range(0, array.GetNumberOfTuples())):
        fh.write(str(array.GetTuple1(i)))
        if idx % 5 == 0 and idx != 0:
            fh.write("\n")
        else:
            fh.write(" ")
    fh.write("\n")

def write_tecplot_ascii(ug,
                        fh,
                        title,
                        time,
                        chunk_id,
                        var_names):
                         
    if not ug.IsHomogeneous():
        print ("WARNING: Unable to convert non-homogenous (cells) VTK unstructured grid")
        return
 
    if not ug.GetNumberOfCells():
        return

    ct = ug.GetCell(0).GetCellType()
    is3d = True
    if ct == vtk.VTK_LINE: 
        tecplot_type = "FELINESEG"
        is3d = False
    elif ct == vtk.VTK_TRIANGLE:
        tecplot_type = "FETRIANGLE"
    elif ct == vtk.VTK_QUAD:
        tecplot_type = "FEQUADRILATERAL"
        is3d = False
    elif ct == vtk.VTK_HEXAHEDRON:
        tecplot_type = "FEBRICK"
    else:
        print ("WARNING: Unknown cell type in VTK unstructured grid, cannot convert")

    array_names = []
    for i in range(0,ug.GetCellData().GetNumberOfArrays()):
      array = ug.GetCellData().GetArray(i)
      array_names.append(array.GetName())

    fh.write("TITLE = " + '"' + title + '"' + "\n")
    if is3d:
        fh.write('VARIABLES = "x", "y", "z"')
    else:
        fh.write('VARIABLES = "x", "y"')

    if var_names:
        if len(array_names) == len(var_names):
            for name in var_names:
                fh.write(', "' + name + '"')
            fh.write("\n")
        else:
            print("WARNING: number of variables in rname file doesn't match datafile") 
            for name in array_names:
                fh.write(', "' + name + '"')
            fh.write("\n")
    else:
        for name in array_names:
            fh.write(', "' + name + '"')
        fh.write("\n")
            
    

    points = ug.GetPoints()
    fh.write("ZONE NODES = " + str(points.GetNumberOfPoints()) + 
             ", ELEMENTS = " + str(ug.GetNumberOfCells()) + 
             ", DATAPACKING = BLOCK" + ", TITLE=\"ZONE %d\""%(chunk_id,))
    
    if array_names:
       fh.write(", VARLOCATION = ([")
       if is3d:
           fh.write(str(4) + "-" + str(3+len(array_names)))
       else:
           fh.write(str(3) + "-" + str(2+len(array_names)))
       fh.write("] = CELLCENTERED)")

    fh.write(", ZONETYPE = " + tecplot_type + ", SOLUTIONTIME = " + str(time) + "\n")

    write_points(points, fh, 0)
    write_points(points, fh, 1)
    if is3d:
        write_points(points, fh, 2)

    for i in range(0,ug.GetCellData().GetNumberOfArrays()):
        array = ug.GetCellData().GetArray(i)
        write_array(array, fh)

    for i in range(0,ug.GetNumberOfCells()):
        pids = ug.GetCell(i).GetPointIds()        
        for j in range(0,pids.GetNumberOfIds()):
            fh.write(str(pids.GetId(j)+1) + " ")
        fh.write("\n")
    
def write_tecplot_bin(ug,
                        filepath,
                        title,
                        time,
                        chunk_id,
                        var_names):
                         
    if not ug.IsHomogeneous():
        print ("WARNING: Unable to convert non-homogenous (cells) VTK unstructured grid")
        return
 
    if not ug.GetNumberOfCells():
        return
     
    
    ct = ug.GetCell(0).GetCellType()
    is3d = True
    if ct == vtk.VTK_LINE: 
        tecplot_type = 1   # FELINESEG
        is3d = False
    elif ct == vtk.VTK_TRIANGLE:
        tecplot_type = 2   # FETRIANGLE
    elif ct == vtk.VTK_QUAD:
        tecplot_type = 3   # FEQUADRILATERAL
        is3d = False
    elif ct == vtk.VTK_HEXAHEDRON:
        tecplot_type = 5   # FEBRICK
    else:
        print ("WARNING: Unknown cell type in VTK unstructured grid, cannot convert")
    # extract data and varnames
    array_names = []
    data = [] # data is stored here
    if ug.GetCellData().GetNumberOfArrays() >0 :
        idata = True
    else:
        idata = False
    if idata:
        for i in range(0,ug.GetCellData().GetNumberOfArrays()):
          array = ug.GetCellData().GetArray(i)
          data0 =  [array.GetTuple1(j)  for j in range(array.GetNumberOfTuples())]
          data.append(data0)
          array_names.append(array.GetName())
    if is3d:
        varnames = ["X", "Y", "Z"]
        varloc = [0]*3
        x = []; y= []; z=[];
    else:
        varnames = ["X", "Y"]  
        varloc = [0]*2
        x = []; y= []
        
    if idata:
        varloc   = varloc+[1]*len(array_names)
        if var_names:
            if len(array_names) == len(var_names):
                varnames = varnames + var_names[:]
            else:
                print("WARNING: number of variables in rname file doesn't match datafile") 
                varnames = varnames + array_names[:]
        else:
            varnames = varnames + array_names[:]
            
    Nvar = len(varnames)
    
    # extract grid information
    points = ug.GetPoints()
    NumPts = points.GetNumberOfPoints()

    if is3d:
        for i in range(points.GetNumberOfPoints()):
            x0 = points.GetPoint(i)
            x.append(x0[0])
            y.append(x0[1])
            z.append(x0[2])
    else:
        for i in range(points.GetNumberOfPoints()):
            x0 = points.GetPoint(i)
            x.append(x0[0])
            y.append(x0[1])
    NumElements = ug.GetNumberOfCells()
        
    
    # other settings
    # DataPacking = 0  # 0 for block, 1 for point, FEM only use block
    SolutionTime = time
    
    # ======================== HEADER SECTION  ====================
    with BinaryFile(filepath,'wb') as fh:
        fh.writechar('#!TDV112')
        fh.writeb(1) 
        fh.writeb(0)
        fh.writestring(title)
        fh.writeb(Nvar) 
        for name in varnames:
            fh.writestring(name)
        
        fh.writeb(299.0,'f')
        # zone name
        fh.writestring('ZONE %d'%(chunk_id,))
        # no parent zone
        fh.writeb(-1)
        # no strandID
        fh.writeb(-2)
        # solution time
        fh.writeb(SolutionTime,'d')
        fh.writeb(-1)
        fh.writeb(tecplot_type)
        # specify var location
        fh.writeb(1) 
        for i in varloc:
            fh.writeb(i)
        #raw local 1-to-1 face neighbors not supplied
        fh.writeb(0)
        NoOfUserDefinedNeighbourConn = 0 
        fh.writeb(NoOfUserDefinedNeighbourConn)
        fh.writeb(NumPts)
        fh.writeb(NumElements)
        #print('Number of Elements',NumElements)
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
        EOH_MARKER = 357.0
        fh.writeb(EOH_MARKER,'f')
        
        # ======================== DATA SECTION  ====================   
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
        if is3d:
            fh.writeb(min(z),'d')
            fh.writeb(max(z),'d')        
        for data0 in data:  
            fh.writeb(min(data0),'d')
            fh.writeb(max(data0),'d')
        
        # write zone
        for i in x:
            fh.writeb(i,'f')
        for i in y:
            fh.writeb(i,'f')
        if is3d:
            for i in z:
                fh.writeb(i,'f')
        if idata:    
            for data0 in data:
                for dat in data0:
                    fh.writeb(dat,'f')
         
        # Zone Connectivity  0-based 
        for i in range(ug.GetNumberOfCells()):
            pids = ug.GetCell(i).GetPointIds()        
            for j in range(pids.GetNumberOfIds()):
                fh.writeb(pids.GetId(j))
