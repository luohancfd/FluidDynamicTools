function [Grid,nblock] = ReadCellGrid( filename,export_flag)
%filename = '0/plot/sphere.t0216.g';
Grid = struct('Block',[],'Size',[],'X',[],'Y',[],'Z',[]);
fid=fopen(filename,'r');
nblock = fscanf(fid,'%d\n',1);
Grid = repmat(Grid,[nblock,1]);
for i = 1:nblock
    Grid(i).Block = i;
    cline = fgetl(fid);
    cline = strsplit(cline);
    Imax = str2double(cline{1});
    Jmax = str2double(cline{2});
    if length(cline) == 3
        Kmax = str2double(cline{3});
    else
        Kmax =1 ;
    end
    Grid(i).Size = [Imax,Jmax,Kmax];
    g = zeros(Imax,Jmax,Kmax,3);
    for dimen = 1:3
        for kloop = 1:Kmax
            for jloop = 1:Jmax
                for iloop = 1:Imax
                    g(iloop,jloop,kloop,dimen) = fscanf(fid,'%e\n',1);
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
%% export correct plot3d format
if nargin == 1
    export_flag = 0;
end

if export_flag
    f_file = [];
    w = strsplit(filename,'.');
    for i = 1:length(w)-1
        f_file = [f_file , w{i},'.'];
    end
    f_file = [f_file,'xyz'];
    fid = fopen(f_file,'w');
    fprintf(fid, '%d\n',nblock);
    for i=1:nblock
        fprintf(fid, '%d ',Grid(i).Size);
        fprintf(fid,'\n');
    end
    for Nloop = 1:nblock
        fprintf(fid,'%e\n',Grid(Nloop).X);
        fprintf(fid,'%e\n',Grid(Nloop).Y);
        fprintf(fid,'%e\n',Grid(Nloop).Z);
    end
    fclose(fid);
end
end

