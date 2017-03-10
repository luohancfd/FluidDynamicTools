function [ Flow,nblock,Nvar ] = ReadFlowEilmer3( filename,export_flag )
%% Read Plot3D Grid function file
% Read plot3d function file exported by Eilmer 3
% Note the one by Eilmer 3 is not same as the official plot3d file
% set export_flag to 1 if you want to get correct plot3d file
Flow = struct('Block',[],'Size',[],'Dat',[]);
%filename = '0/plot/sphere.t0216.f';
fid=fopen(filename,'r');
nblock = fscanf(fid,'%d\n',1);
Flow = repmat(Flow,[nblock,1]);
dimen = 3;
for i = 1:nblock
    Flow(i).Block = i;
    cline = fgetl(fid);  %get number of grid and number of field
    cline = strsplit(cline);
    Imax = str2double(cline{1});  %Imax
    Jmax = str2double(cline{2});   %Jmax
    Nvar = str2double(cline{end}); %Nvar ,number of variable, same as nam file
    if length(cline) == 4
        Kmax = str2double(cline{3});  %if it's 3D file
    else
        Kmax = 1 ;
        dimen = 2;
    end
    Flow(i).Size = [Imax,Jmax,Kmax,Nvar];
    Flow(i).Dat = zeros([Imax,Jmax,Kmax,Nvar]);
    for Nloop = 1:Nvar
        for kloop = 1:Kmax
            for jloop = 1:Jmax
                for iloop = 1:Imax
                    Flow(i).Dat(iloop,jloop,kloop,Nloop) = fscanf(fid,'%e\n',1);
                end
            end
        end
    end
end
fclose(fid);
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
    f_file = [f_file,'fun'];
    fid = fopen(f_file,'w');
    fprintf(fid, '%d\n',nblock);
    for i=1:nblock
        fprintf(fid, '%d %d %d %d\n',[Flow(i).Size]);
    end
    for i = 1:nblock
        fprintf(fid,'%e\n',Flow(i).Dat);
    end
    fclose(fid);
end
end

