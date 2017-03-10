function [ NewGrid,NewFlow ] = ConcateBlocks( map,Grid,varargin )
%% Concate multiblock together, assume the whole space is uniformely divided
% Default concate blocks with one overlap (skip=0)
% First: find the correct direction of i j and k, input blockid in the
% following way
% num_k = 1;
% num_j = 4;
% num_i = 5;
% map = zeros(num_i,num_j,num_k);
% for k = 1:num_k   
%     for j = 1:num_j
%        for i = 1:num_i 
%             map(i,j,k) = blockid
%             w = w+1;
% end; end; end; 
% ex: map = reshape(1:20,4,[])'
% if "Flow" is given, also concate flow, flow should have same size as grid
% if "skip" is given, skip = 1 for cell centered 0 for grid centered

flag = [0 0]; %default value: no flow in, vertex-centered grid
if length(varargin) >1
    for i = 1:2:length(varargin)
        if strcmp(varargin{i},'Flow')
            flag(1) = 1;
            Flow = varargin{i+1};
        elseif strcmp(varargin{i},'Skip')
            flag(2) = varargin{i+1};
        end
    end
end

skip = flag(2);
NewGrid = struct('Block',1,'Size',[],'X',[],'Y',[],'Z',[]);
NewFlow = struct('Block',1,'Size',[],'Dat',[]);
Nblock = length(Grid);
Corner = zeros(Nblock,3);   %min y, min x, max z for each block
for i = 1:Nblock
    Corner(i,:) = [ Grid(i).X(1,1) ,Grid(i).Y(1,1),Grid(i).Z(1,1)];
end
Range = zeros(Nblock,6);  %[imin, imax, jmin, jmax, kmin kmax] of each block


[isize,jsize,ksize] = size(map);
k0 = 1;
for kloop = 1:ksize  % sweep over blocks with same z
    kblock = map(:,:,kloop);  %index of blocks on same layer
    kblock = reshape(kblock,[],1);   % change to 1D index list
    s = Grid(kblock(1)).Size;  %assume same Size(3)= kmax-kmin+1 for all kblock
    Range(kblock,5) = k0;
    Range(kblock,6) = k0+s(3)-1;
    k0 = k0 + s(3) -1 + skip;
    
    j0 = 1;
   for jloop = 1:jsize     % sweep over blocks with same y
        jblock = map(:,jloop,kloop);
        s = Grid(jblock(1)).Size;  %assume same Size(1)=imax -imin+1 for all jblock
        Range(jblock,3) = j0;
        Range(jblock,4) = j0+s(2)-1;
        j0 = j0 + s(2) -1 + skip;
        
        i0 = 1;
        for iloop = 1:isize
            iblock = map(iloop,jloop,kloop);
            s = Grid(iblock).Size;
            Range(iblock,1)= i0;
            Range(iblock,2)= i0 + s(1)-1;
            i0 = i0 + s(1) -1 +skip;
        end
    end
end
Imax = max(Range(:,4));
Jmax = max(Range(:,2));
Kmax = max(Range(:,6));
NewGrid.Size = [Imax,Jmax,Kmax];
[NewGrid.X, NewGrid.Y,NewGrid.Z ] = deal(zeros(NewGrid.Size));

for i = 1:Nblock
    Irange = Range(i,1) : Range(i,2);
    Jrange = Range(i,3) : Range(i,4);
    Krange = Range(i,5) : Range(i,6);
    NewGrid.X(Irange,Jrange,Krange) = Grid(i).X;
    NewGrid.Y(Irange,Jrange,Krange) = Grid(i).Y;
    NewGrid.Z(Irange,Jrange,Krange) = Grid(i).Z;
end
    

% concate flow
if flag(1)
    s = size(Flow(1).Dat);
    Nvar = s(end);
    NewFlow.Size = [Imax,Jmax,Kmax,Nvar];
    NewFlow.Dat = zeros(Imax,Jmax,Kmax,Nvar);
    Nblock = length(Flow);
    for i = 1:Nblock
        Irange = Range(i,1) : Range(i,2);
        Jrange = Range(i,3) : Range(i,4);
        Krange = Range(i,5) : Range(i,6);
        NewFlow.Dat(Irange,Jrange,Krange,:) = Flow(i).Dat;
    end
else
    NewFlow = false;
end
% mesh(NewGrid.X,NewGrid.Y,NewGrid.Z)
% hold on
% view([0 90])
% axis equal
end