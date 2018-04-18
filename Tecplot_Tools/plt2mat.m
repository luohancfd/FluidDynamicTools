%% Author: Han Luo
%% Read binary format tecplot dat
fname = 'bc.plt';
clc
fclose(fid)
fid = fopen(fname,'rb');
tdata = [];
magic_number = plt_read_string(fid,8,'int8');
int_type = fread(fid,1,'int32');
if int_type ~= 1
    error('Int type is not correct');
end
filetype = fread(fid,1,'int32');
tdata.title = plt_read_string(fid);
tdata.Nvar =  fread(fid,1,'int32');
tdata.varnames = cell(1,tdata.Nvar);
fprintf('=============================\nVariable names:\n=============================\n');
for i  = 1:tdata.Nvar
    tdata.varnames{i} = plt_read_string(fid);
end
for i = 1:tdata.Nvar
    fprintf('   %s\n',tdata.varnames{i});
end
%% zones section
[icheck, isec] = plt_check_section(fid,[299,399,499,599,357.0]);
% dummy, do nothing
nzone = 0;
nlines = 0; nsurfs = 0; ncubes = 0;
nFElines = 0; nFEsurfaces = 0; nFEvolumes=0;
while true
    fseek(fid,-4,'cof'); % jump back to check header of zone
    
    zonecheck = fread(fid,1,'float32');
    if zonecheck ~= 299
        error('Something wrong with starter of zone section');
    end
    
    nzone = nzone +1;
    izone = nzone;
    fprintf('=============================\nReading info of zone %d...\n=============================\n',nzone);
    zonename = plt_read_string(fid);
    fprintf('    Zonename: %s\n',zonename);
    parent_zone = fread(fid,1,'int32');
    strandID = fread(fid,1,'int32');
    fprintf('    StrandID: %d\n',strandID);
    solutiontime = fread(fid,1,'float64');
    fprintf('    Solutiontime: %f\n',solutiontime);
    dummy = fread(fid,1,'int32');
    [zonetype,izonetype] = plt_read_zonetype(fid);
    fprintf('    Zonetype:  %s\n',zonetype);
    % +-----------+
    % | INT32     |       Specify Var Location.  0 = Don't specify, all data is
    % +-----------+       located at the nodes.  1 = Specify
    %
    %   0 ----  Not specifying
    %   1 ----  Specify
    have_varloc = logical(fread(fid,1,'int32'));
    varloc = zeros(1,length(tdata.Nvar));
    if have_varloc
        varloc = fread(fid,tdata.Nvar,'int32');
        raw_1_to_1 = fread(fid,1,'int32');
    end
    fprintf('    Varloc: ')
    for i = 1:length(varloc)
        fprintf('%d:%d ',i,varloc(i));
    end
    fprintf('\n')
    % if “number of miscellaneous user-defined
    NoOfUserDefinedNeighbourConn = fread(fid,1,'int32');
    if (NoOfUserDefinedNeighbourConn ~=0)
        user_defined_face_neighbor_mode = fread(fid,1,'int32');
        if itype >0
            finite_element_face_neighbors_mis = fread(fid,1,'int32');
        end
    end
    %if Ordered Zone:
    Imax = 0; Jmax = 0; Kmax = 0;
    if izonetype == 0
        Imax = fread(fid,1,'int32');
        Jmax = fread(fid,1,'int32');
        Kmax = fread(fid,1,'int32');
        fprintf('    Imax: %d, Jmax: %d, Kmax: %d\n',Imax,Jmax,Kmax);
        if Jmax == 1 && Kmax == 1
            nlines = nlines +1; nID = nlines; 
            tdata.lines(nlines).x = zeros(1,Imax);
            ztname = 'lines';
        elseif Kmax == 1
            nsurfs = nsurfs +1; nID = nsurfs;
            tdata.surfaces(nsurfs).x = zeros(1,Imax);  
            tdata.surfaces(nsurfs).y = zeros(1,Jmax);  
            ztname = 'surfaces';
        else
            ncubes = ncubes +1; nID = ncubes;
            tdata.cubes(ncubes).x = zeros(1,Imax);  
            tdata.cubes(ncubes).y = zeros(1,Jmax);             
            tdata.cubes(ncubes).z = zeros(1,Kmax); 
            ztname = 'cubes';
        end
        tdata.(ztname)(nID).zonename = zonename;
        tdata.(ztname)(nID).varloc = varloc;
        tdata.(ztname)(nID).solutiontime = solutiontime;
        tdata.(ztname)(nID).strandID = strandID;
    end
    %if FE Zone:
    if izonetype > 0
        NumPts = fread(fid,1,'int32');
        if strcmp(zonetype,'FEPOLYGON') || strcmp(zonetype,'FEPOLYHEDRON')
            NumFaces = fread(fid,1,'int32');
            NumNodes = fread(fid,1,'int32');
            NumBoundSurf = fread(fid,1,'int32');
            NumBoundConn = fread(fid,1,'int32');
            error('Not support FEPOLYGON and FEPOLYHEDRON yet');
        end        
        NumElements = fread(fid,1,'int32');
        fprintf('    NumPts: %d, NumElements: %d\n',NumPts,NumElements);
        ICellDim = fread(fid,1,'int32');
        JCellDim = fread(fid,1,'int32');
        KCellDim = fread(fid,1,'int32');
        
        if strcmp(zonetype,'FEQUADRILATERAL')
            ztname = 'FEsurfaces'; conn_n = 4;
            nFElines = nFElines+1; nID = nFElines;
        elseif strcmp(zonetype,'FELINESEG')
            ztname = 'FElines'; conn_n = 2;
            nFEsurfaces = nFEsurfaces+1; nID = nFEsurfaces;
        elseif strcmp(zonetype,'FEBRICK')
            ztname = 'FEvolumes'; conn_n = 8;
            nFEvolumes = nFEvolumes+1; nID = nFEvolumes;
        else
            error('Not support some FEM yet');
        end
        tdata.(ztname)(nID).e2n = zeros(NumElements,conn_n);
        tdata.(ztname)(nID).zonename = zonename;
        tdata.(ztname)(nID).varloc = varloc;
        tdata.(ztname)(nID).solutiontime = solutiontime;
        tdata.(ztname)(nID).strandID = strandID;
    end
    %For all zone types (repeat for each Auxiliary data name/value pair):
    have_auxdata = logical(fread(fid,1,'int32'));
    aux_names = {}; aux_values = {};
    if have_auxdata
        N_aux = 0;
        while true
            [icheck, isec] = plt_check_section(fid,[299,399,499,599,357.0]);
            if sum(icheck(2:4))>0
                error('Not support geom/text/custom section');
            elseif icheck(1) || icheck(5)
                break
            else
                N_aux = N_aux + 1;
                aux_name   = plt_read_string(fid);
                aux_format = fread(fid,1,'int32'); % should be 0
                aux_value  = plt_read_string(fid);
                aux_names{N_aux} = aux_name;
                aux_values{N_aux} = aux_value;
            end
        end
        tdata.(ztname)(nID).auxname = aux_names;
        tdata.(ztname)(nID).aux_value = aux_values;
    else
        [icheck, isec] = plt_check_section(fid,[299,399,499,599,357.0]);
        if sum(icheck(2:4))>0
            error('Not support geom/text/custom section');
        end
    end
    if icheck(5)
        break  %now go to data section
    else
        continue
    end
end 
    
fprintf('=============================\nFinish loading zone info, %d zones\n=============================\n',nzone);    
%% data section


