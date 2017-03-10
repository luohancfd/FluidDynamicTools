function [Grid,nblock, dimen] = ReadGrid(filename)
%% Read Plot3D Grid with node centered
% Created by Han Luo to Read plot3d format grid exported by Eilmer3 
% This file is used to read multiple grid xyz_file
% I have only tested for 2d grid, should also work for 3D grid, use it at
% your own responsbility
% filename = '0/plot/sphere.t0216.grd'
Grid = struct('Block',[],'Size',[],'X',[],'Y',[],'Z',[]);
fid=fopen(filename,'r');
cline = fgetl(fid);
nblock=sscanf(cline,'%d',1);
Grid = repmat(Grid,[nblock,1]);
for i = 1:nblock
    Grid(i).Block = i;
    cline = fgetl(fid);
    cline = sscanf(cline,'%d');
    Imax = cline(1);
    Jmax = cline(2);
    if length(cline) == 3
        Kmax = cline(3);
        dimen = 3;
    else
        Kmax = 1 ;
        dimen = 2;
    end
    Grid(i).Size = [Imax,Jmax,Kmax];
end
for i = 1:nblock
    Kmax = Grid(i).Size(3); Jmax = Grid(i).Size(2); Imax = Grid(i).Size(1);
    g = zeros(Imax,Jmax,Kmax,3);
    for loop = 1:dimen
        for kloop = 1:Kmax
            for jloop = 1:Jmax
                for iloop = 1:Imax
                    w = fgetl(fid);
                    g(iloop,jloop,kloop,loop) = sscanf(w,'%e',1);
                end
            end
        end
    end
    Grid(i).X = g(:,:,:,1);
    Grid(i).Y = g(:,:,:,2);
    Grid(i).Z = g(:,:,:,3);
end
fclose(fid);
% for i=1:nblock
%     mesh(Grid(i).X,Grid(i).Y,Grid(i).Z)
%     hold on
% end
% view([0 90])
% axis equal
end