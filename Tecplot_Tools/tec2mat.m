function [zone,VARlist,auxdata] = tec2mat(fname,varargin)
%% load data file in Tecplot ASCII format
%% Author: Mach6
% Github repository: https://github.com/luohancfd/FluidDynamicTools/tree/master/Tecplot_Tools
% Method of calling:
%   [zone,VARlist] = tec2mat(fname,debug);
%   fname(in): name of data file
%   optional flags:
%          "debug" print zone information out (default: false)
%          "safe" check the existence of scientific notation like 100-100 (default: false)
%        For IJK ordered date only:
%          "nopassive" include passive variables (0) in the data field (default: true)
%          "passive"   not include passive variables in the data field
% -------------------------------------------------------------------------
% ChangeLog:
% v1: support IJK data without mixing data format (node with cell-centered)
%     support Finite Element grid
% v2: add a "safe" optional flag to handle number like 100-100 (100E-100)
%--------------------------------------------------------------------------
debug = false;
nopassive = true;
safe = false;
if nargin > 1
    for i = 1:length(varargin)
        if strcmp(varargin{i},'debug')
            debug = true;
        elseif strcmp(varargin{i},'nopassive')
            nopassive = true;
        elseif strcmp(varargin{i},'passive')
            nopassive = false;
        elseif strcmp(varargin{i},'safe')
            safe = true;
        else
            error('Unknown option: %s',varargin{i})
        end
    end
end
%% First test if the file has "strange" scientific notations like  2.3-100 (2.3E-100)
if ~safe
    IMY_FSCANF = false;
else
    filetext = fileread(fname);
    strange = regexp(filetext,'(?<!(?:(,|\[)))-?(\d+(\.\d+)?|\.\d+)([-+]\d+)','match','once');
    if isempty(strange)
        IMY_FSCANF = false;
    else
        IMY_FSCANF = true;
    end
end
%% load tecplot ASCII file and write result to zone and VARList
fid = fopen(fname,'r');
NVAR = 0;
%% read file untill we get a line of variables name
IVAR = false;
varline = '';
while ~feof(fid)
    pos = ftell(fid);  %line above ZONE
    line = fgetl(fid);
    line = strtrim(line);
    if ~isempty(line)
        if (~strcmp(line(1),'#'))
            if contains(line,'VARIABLES')
                varnames = regexp(line,'VARIABLES\s*=\s*(.*)','tokens');
                varline = [varline, varnames{1}{1}];
                IVAR = true;
            elseif (contains(line,'ZONE') || contains(line,'DATASETAUXDATA'))
                fseek(fid,pos,-1); %jump out at position above zone
                break
            elseif (IVAR)
                varline = [varline, line];
            end
        end
    end
end

if debug
    fprintf('===========================================\n');
    fprintf('Read file: %s',fname);
end
if IVAR
    varline0 = varline;
    % replace all delimiter by a special symbol
    varline = strrep(varline,'''','$#');
    varline = strrep(varline,'"','$#');
    varline = strrep(varline,',','$#');
    VARlist0 =strsplit(varline,'$#');
    VARlist = {};
    for i = 1:length(VARlist0)
        if ~isempty(strtrim(VARlist0{i}))
            VARlist = [VARlist, VARlist0{i}];
        end
    end
    
    NVAR = length(VARlist);
    if debug
        fprintf('The line defining variables:\n%s\n',varline0);
        fprintf('The variables we get:\n');
        fprintf('%s, ',VARlist{:});
        fprintf('\n');
    end
end

nzone = 0;
zone = [];
auxdata = {};
if ~IVAR
    error("there is no variable list")
else
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~isempty(line)
            if strcmp(line(1),'#')
                continue
            end
        end
        if contains(upper(line),'DATASETAUXDATA')
            if strcmpi(line(1:14),'DATASETAUXDATA')
                line=strtrim(line(15:end));
                temp = strsplit(line,'=');
                temp(2)=strtrim(temp(2));
                temp(2)=strrep(temp(2),'"','');
                temp(2)=strrep(temp(2),'''','');
                auxdata = [auxdata;{temp(1),temp(2)}];
            end
        end
            
        if contains(line,'ZONE')
            nzone = nzone + 1;
                           
            zoneline = regexp(line,'ZONE(.*)','tokens');
            zoneline = strtrim(zoneline{1}{1});
            
            while ~feof(fid)
                pos = ftell(fid);  %position above data
                line = strtrim(fgetl(fid));
                if ~isempty(line)
                    if strcmp(line(1),'#')
                        continue
                    end
                end
                if contains(line,'=')
                    zoneline = [zoneline,' ',line];
                else
                    % arrive at data
                    fseek(fid,pos,-1); %jump out at position above zone
                    break
                end
            end
            
            zoneline0 = zoneline;
            % preprocess zone property, replace equal sign in quotes with a
            % special symbol
            quoted =  regexp(zoneline,'"(.+?)"','tokenExtents');
            quoteds = regexp(zoneline,'''(.+?)''','tokenExtents');
            quoted = [quoted quoteds];
            
            zoneline = '';
            iendp = 1;
            for i = 1:length(quoted)
                istart = quoted{i}(1);
                iend = quoted{i}(2);
                if istart-1 >= iendp
                    zoneline = [zoneline,zoneline0(iendp:istart-1)];
                end
                zoneline = [zoneline,strrep(zoneline0(istart:iend),'=','$#')];
                if iend == length(zoneline0)
                    break
                end
                iendp = iend+1;
            end
            if iend ~= length(zoneline0)
                zoneline = [zoneline, zoneline0(iendp:end)];
            end
                
            % get zone properties
            [ZONEprop.name,ZONEprop.start,ZONEprop.end] = regexp(zoneline,'\w+(?=\s*=)\s*=','match','start','end');
            Nzoneprop = length(ZONEprop.name);
            ZONEprop.prop = cell(1,Nzoneprop);
            for i = 1:Nzoneprop-1
                ZONEprop.name{i} = strtrim(ZONEprop.name{i}(1:end-1)); %remove "=" sign
                ZONEprop.name{i} = upper(ZONEprop.name{i});  %upper case
                if strcmpi(ZONEprop.name{i},'Nodes')
                    ZONEprop.name{i} = 'N';
                end
                if strcmpi(ZONEprop.name{i},'Elements')
                    ZONEprop.name{i} = 'E';
                end
                
                ZONEprop.prop{i} = strtrim(zoneline(ZONEprop.end(i)+1:ZONEprop.start(i+1)-1));
                if strcmp(ZONEprop.prop{i}(end),',')
                    ZONEprop.prop{i} = ZONEprop.prop{i}(1:end-1);
                end
            end
            ZONEprop.name{Nzoneprop} = strtrim(ZONEprop.name{Nzoneprop}(1:end-1)); 
            ZONEprop.prop{Nzoneprop} = zoneline(ZONEprop.end(Nzoneprop)+1:end);
            
            
            % reformat properties
            Prop = [];
            for i = 1:Nzoneprop
                if isstrprop(ZONEprop.prop{i},'digit')
                    Prop.(ZONEprop.name{i}) = str2num(ZONEprop.prop{i});
                else
                    Prop.(ZONEprop.name{i}) = strtrim(strrep(ZONEprop.prop{i},'$#','='));
                end
            end
            clear ZONEprop;
            
            if ~isfield(Prop,'T')
                Prop.T = sprintf('ZONE %d',nzone);
            else
                if strcmp(Prop.T(1),'"') && strcmp(Prop.T(end),'"')
                    Prop.T = Prop.T(2:end-1);
                elseif strcmp(Prop.T(1),'''') && strcmp(Prop.T(end),'''')
                    Prop.T = Prop.T(2:end-1);
                end
            end
            Prop.title = Prop.T;
                
            localVARindex = 1:NVAR;
            if isfield(Prop,'PASSIVEVARLIST')
                list = regexp(Prop.PASSIVEVARLIST,'[(.*)\]','tokens');
                if ~isempty(list)
                    PASSIVEVARLIST = strl2list(list{1}{1});
                    Prop.PASSIVEVARLIST = PASSIVEVARLIST;
                    localVARindex = setdiff(localVARindex,Prop.PASSIVEVARLIST);
                    ZONEVARlist = VARlist(localVARindex);
                else
                    ZONEVARlist = VARlist;
                end
            else
                Prop.PASSIVEVARLIST =[];
                ZONEVARlist = VARlist;
            end
            localNVAR = length(ZONEVARlist);
            
            VARLOC = zeros(1,NVAR);  %0 for node, 1 for cell
            Prop.IVARLOC = false;
            if isfield(Prop,'VARLOCATION')
                varloc = regexp(Prop.VARLOCATION,'\((.*)\)','tokens');
                if ~isempty(varloc)
                    Prop.IVARLOC = true;
                    varloc = strsplit(varloc{1}{1},{','});
                    nodalvar = []; cellvar = [];
                    for i = 1:length(varloc)
                        if contains(varloc{i},'NODAL')
                            nodalvar =  regexp(varloc{i},'\[(.*)\]','tokens');
                            nodalvar = strl2list(nodalvar{1}{1});
                        elseif contains(varloc{i},'CELLCENTERED')
                            cellvar =  regexp(varloc{i},'\[(.*)\]','tokens');
                            cellvar = strl2list(cellvar{1}{1});
                        end
                    end
                    
                    for i = 1:NVAR
                        if ismember(i,cellvar)
                            VARLOC(i) = 1;
                        end
                    end
                end
            end
            Prop.VARLOC =VARLOC; 

            if isfield(Prop,'I')
                % IJK ordered ZONE
                if isfield(Prop,'ZONETYPE')
                    if ~strcmpi(Prop.ZONETYPE,'ORDERED')
                        error(['ZONE ',num2str(nzone),' type is not ordered']);
                    end
                else
                    Prop.ZONETYPE = 'ORDERED';
                end
                idata = Prop.I;
                jdata = 1; kdata = 1;
                if isfield(Prop,'J')
                    jdata = Prop.J;
                end                
                if isfield(Prop,'K')
                    kdata = Prop.K;
                end
                if debug
                    fprintf('ZONE %d is a IJK ordered zone\n',nzone);
                    fprintf('ZONE %d Size: %d x %d x %d NVar: %d\n',nzone,idata,jdata,kdata,localNVAR);
                end
                ndata = idata*jdata*kdata;
                if ~isfield(Prop,'DATAPACKING')
                    Prop.DATAPACKING = 'POINT';
                end
                % load zone              
                if Prop.IVARLOC 
                    error('Don''t support IJK with varloc now');
                end
                    
                % let's load the data   
                data = myfscanf(fid,'%f',ndata*localNVAR,IMY_FSCANF);
                if strcmpi(Prop.DATAPACKING,'POINT')
                    data = reshape(data,localNVAR,idata,jdata,kdata);
                    data2 = zeros(NVAR,idata,jdata,kdata);
                    if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                        for i = 1:length(localVARindex)
                            data2(localVARindex(i),:,:,:) = data(i,:,:,:);
                        end
                    end
                elseif strcmpi(Prop.DATAPACKING,'BLOCK')
                    data = reshape(data,idata,jdata,kdata,localNVAR);
                    data2 = zeros(idata,jdata,kdata,NVAR);
                    if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                        for i = 1:length(localVARindex)
                            data2(:,:,:,localVARindex(i)) = data(:,:,:,i);
                        end
                    end
                else
                    error('DATAPACKING not supported');
                end
                
                if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                    data = squeeze(data2);
                else
                    data = squeeze(data);
                end
                    
                if (jdata == 1 && kdata == 1)
                    data = data'; %if x-y line data, transpose the data for easier reading
                end

                Propfields = fieldnames(Prop);
                for i =1:length(Propfields)
                    zone(nzone).(Propfields{i}) = Prop.(Propfields{i});
                end
                zone(nzone).data = data;
                if nopassive
                    zone(nzone).ZONEVARlist = ZONEVARlist;
                end
                
                clear Prop;

            elseif ((isfield(Prop,'N') && isfield(Prop,'E')) || ( isfield(Prop,'NODES') && isfield(Prop,'ELEMENTS')))
                % load finite element zone
                if (isfield(Prop,'N') && isfield(Prop,'E') )
                    % load zone information
                    Nnode = Prop.N;
                    Ecell = Prop.E;
                elseif ( isfield(Prop,'NODES') && isfield(Prop,'ELEMENTS'))
                    Nnode = Prop.NODES;
                    Ecell = Prop.ELEMENTS;
                else
                    error("No node information found");
                end
                
                if ~isfield(Prop,'DATAPACKING')
                    Prop.DATAPACKING = 'BLOCK';
                end
                if debug
                    fprintf('ZONE %d is a finite element zone\n',nzone);
                    fprintf('ZONE %d size: Node: %d Element: %d\n',nzone,Nnode,Ecell);
                end

                % start loading data
                data = repmat(struct('name',[],'data',[]),localNVAR,1);
                if strcmpi(Prop.DATAPACKING,'BLOCK')
                    for i = 1:localNVAR
                        data(i).name = ZONEVARlist{i};
                        ivarloc = Prop.VARLOC(localVARindex(i));
                        if ivarloc == 0
                            ndata = Nnode;
                        else
                            ndata = Ecell;
                        end
                        data(i).data = myfscanf(fid,'%e',ndata,IMY_FSCANF);
                    end
                elseif (strcmpi(Prop.DATAPACKING,'POINT') && ~Prop.IVARLOC)
                    %ensure not VARLOCATION and in point mode
                    data2 = myfscanf(fid,'%e',Nnode*localNVAR,IMY_FSCANF);
                    data2 = reshape(data2,localNVAR,Nnode);
                    
                    for i = 1:localNVAR
                        data(i).name = ZONEVARlist{i};
                        data(i).data = data2(i,:);
                    end
                end
                
                % loading connect information
                if strcmpi(Prop.ZONETYPE,'FEQUADRILATERAL')
                    conn = fscanf(fid,'%d',Ecell*4);
                    conn = reshape(conn,4,[]);
                elseif strcmpi(Prop.ZONETYPE,'FETRIANGLE')
                    conn = fscanf(fid,'%d',Ecell*3);
                    conn = reshape(conn,3,[]);
                elseif strcmpi(Prop.ZONETYPE,'FEBRICK')
                    conn = fscanf(fid,'%d',Ecell*8);
                    conn = reshape(conn,8,[]);
                elseif strcmpi(Prop.ZONETYPE,'FETETRAHEDRON')
                    conn = fscanf(fid,'%d',Ecell*4);
                    conn = reshape(conn,4,[]);   
                else
                    warning('FEPOLYGON and FEPOLYHEDRAL are not verfied')
                    conn = cell(1,Ecell);
                    for i = 1:Ecell
                        line = fgetl(fid);
                        conn{i} = sscanf(line,'%d');
                        conn{i} = conn{i}';
                    end
                end
                conn = conn';
                
                Propfields = fieldnames(Prop);
                for i =1:length(Propfields)
                    zone(nzone).(Propfields{i}) = Prop.(Propfields{i});
                end
                zone(nzone).data = data;
                zone(nzone).conn = conn;      
                zone(nzone).ZONEVARlist = ZONEVARlist;
                clear Prop;
            else
                error('format not supported yet')
            end
        end
    end
end            

zone = zone(1:nzone);
                                
fclose(fid);
if debug
    fprintf('===========================================\n');
end
end

function [list ] = strl2list(temp)
temp = strsplit(temp,',');
list = [];
for i = 1:length(temp)
    if ~contains(temp{i},'-')
        list = [list, str2num(temp{i})];
    else
        c = strrep(temp{i},'-',' ');
        c = str2num(c);
        list = [list, c(1):c(2)];
    end
end
end      

function [temp2] = myfscanf(fid,format,n,IMY_FSCANF)
%% A modified fscanf inorder to read fortran format scientific notation number
% like 1.234-100 = 1.234E-100
if IMY_FSCANF
    temp = textscan(fid,'%s',n,'CommentStyle','#');
    temp = temp{1};
    temp2 = regexprep(temp,'(-?)(\d+(.\d+)?|.\d+)([-+]\d+)','$1$2E$3');
    temp2 = cellfun(@str2double,temp2,'UniformOutput', false);
    temp2 = cell2mat(temp2);
else
    temp2 = fscanf(fid,format,n);
end
end
