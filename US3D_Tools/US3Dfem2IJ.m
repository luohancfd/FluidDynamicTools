function [tdata,surfaces,VARlist] = US3Dfem2IJ (varargin)
%% Author: Han Luo
% conver US3D finite element format surf solution to IJ ordered data
cellfind = @(x,y) find(ismember(x,y));
%% calculate size of domain
if mod(length(varargin),2) ~= 0
    error('Wrong number of arguments')
end
h5grid = 'grid.h5';
FEMdat = 'bc.dat';
MultiZone = 0;
if length(varargin)>=2
    for i = 1:2:length(varargin)-1
        switch varargin{i}
            case 'h5grid'
                h5grid = varargin{i+1};
            case 'FEMdat'
                FEMdat = varargin{i+1};
            case 'BCtype'
                switch varargin{i+1}
                    case 'wall'
                        bctype = 3;
                    case 'symm'
                        bctype = 7;
                    otherwise
                        error('Wrong BCtype')
                end
            case 'MultiZone'
                MultiZone = varargin{i+1};
                % ----Option 1
                % this should be a matrix telling how the zones are
                % combined into the FEMdat
                % first row  -- imax
                % second row -- jmax
                % third row  -- direction
                % fourth row  -- give the order in x direction
                % for example the following one is combined by two zones
                % the left part is zone 1 (2*4), the right one is zone 2(3*4)
                % |---|---|---|---|---|
                % | 8 | 4 |10 |11 |12 |
                % |---|---|---|---|---|
                % | 7 | 3 | 7 | 8 | 9 |
                % |---|---|---|---|---|
                % | 6 | 2 | 4 | 5 | 6 |
                % |---|---|---|---|---|
                % | 5   1 | 1 | 2 | 3 |
                % |---|---|---|---|---|
                % the MultiZone should be 
                % [3 2; 4 4; 1 2; 2 1];
                % note: in the matrix, from left to right, it should follow
                % the cell orders
                % ----- Option 2
                % give the thickness of wall in z direction and guess Imax
                % from that
            otherwise
                error(['Wrong argument ',varargin{i}])
        end
    end
end
% load FEM data
[FEMzones,VARlist] = tec2mat(FEMdat);
zonename = strrep(FEMzones.T,'"','');
zonename = strrep(zonename,'''','');

% load H5Grid data
zoneinfos = h5read(h5grid,'/zones/zdefs');
zoneinfos = zoneinfos';
zonenames = h5read(h5grid,'/zones/znames');
zonenames = cellfun(@strtrim,zonenames,'UniformOutput',0);
izone = cellfind(zonenames,zonename);
izone_size = zoneinfos(izone,5)-zoneinfos(izone,4)+1;

% calculate Imax, Jmax
if size(MultiZone,1) == 1
    jj = MultiZone;
    % find Imax from wall
    iwall = cellfind(zonenames,'wall');
    MultiZone = zeros(3,1);
    MultiZone(1) = (zoneinfos(iwall,5)-zoneinfos(iwall,4)+1)/jj;
    MultiZone(2) = izone_size/MultiZone(1);
    if mod(izone_size,MultiZone(1)) ~= 0 
        error('You should use MultiZone, the program couldn''t guess');
    end
    MultiZone(3) = 1;
    MultiZone(4) = 1;
end
[~,Nmulti] = size(MultiZone);

% set i range for each zone
MultiZoneX = zeros(2,Nmulti);
ix = 1;
for i = 1:Nmulti
    ii = find(MultiZone(4,:) == i);
    MultiZoneX(1,ii) = ix;
    MultiZoneX(2,ii) = ix+MultiZone(1,ii)-1;
    ix = MultiZoneX(2,ii)+1;
end

MSUM = sum(MultiZone(1,:));
% Imax, Jmax is the number of cells in each direction
Imax = MSUM;
Jmax = MultiZone(2,1);

% Now they are number of points
Imax = Imax + 1;
Jmax = Jmax + 1;
surfaces.x = zeros(Jmax,Imax);
surfaces.y = zeros(Jmax,Imax);
surfaces.z = zeros(Jmax,Imax);
surfaces.v = zeros(Jmax-1,Imax-1,length(VARlist)-3);
%% convert to IJ order
tdata = [];
tdata.varnames = VARlist;
zfileds = {FEMzones.data.name};
ix = cellfind(zfileds,'x');
iy = cellfind(zfileds,'y');
iz = cellfind(zfileds,'z');

datav = zeros(length(VARlist)-3,(Imax-1)*(Jmax-1));
for i = 4:length(VARlist)
    datav(i-3,:) = FEMzones.data(i).data';
end
Nvar = length(VARlist)-3;

icell0 = 1;
for i = 1:Nmulti
    [ x,y,v,z] = fem2ij(FEMzones.conn,datav,Nvar,icell0,...
        FEMzones.data(ix).data,FEMzones.data(iy).data,FEMzones.data(iz).data,...
        MultiZone(1,i)+1,MultiZone(2,i)+1,MultiZone(3,i));
    ibegin = MultiZoneX(1,i); iend = MultiZoneX(2,i);
    surfaces.x(:,ibegin:iend+1) = x;
    surfaces.y(:,ibegin:iend+1) = y;
    surfaces.z(:,ibegin:iend+1) = z;
    surfaces.v(:,ibegin:iend,:) = v;
    icell0 = icell0 + MultiZone(1,i)*MultiZone(2,i);
end

%contourf(surfaces.x,surfaces.y,surfaces.v(:,:,1))

%% up to now, the direction of matrix is
%    ----->i
%    |
%    j
% now, we transpose it
tdata = [];
tdata.varnames = VARlist;
tdata.surfaces.x = surfaces.x';
tdata.surfaces.y = surfaces.y';
tdata.surfaces.z = surfaces.z';
tdata.surfaces.v = zeros(Nvar,Imax-1,Jmax-1);
for i = 1:length(VARlist)-3
    tdata.surfaces.v(i,:,:) = surfaces.v(:,:,i)';
end
tdata.surfaces.order = 3;
tdata.vformat = ones(1,length(VARlist)+3)*2;
tdata.surfaces.varloc = 1;
tdata.surfaces.zonename = FEMzones.T;
%mat2tecplot(tdata,'bc_ij.plt')
end




