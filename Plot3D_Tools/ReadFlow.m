function [ Flow,nblock,Nvar ] = ReadFlow( filename )
%% Read Plot3D Grid function file
% Read ordinary plot3d function file 
% Note the one by Eilmer 3 is not same as the official plot3d file
% set export_flag to 1 if you want to get correct plot3d file
Flow = struct('Block',[],'Size',[],'Dat',[]);
%filename = '0/plot/sphere.t0216.f';
fid=fopen(filename,'r');
nblock = fscanf(fid,'%d\n',1);
Flow = repmat(Flow,[nblock,1]);
for i=1:nblock
     Flow(i).Block = i;
     cline = fgetl(fid);  %get number of grid and number of field
     cline = strsplit(cline);
     Imax = str2double(cline{1});  %Imax
     Jmax = str2double(cline{2});   %Jmax
     Kmax = str2double(cline{3});
     Nvar = str2double(cline{end}); %Nvar ,number of variable, same as nam file
     Flow(i).Size = [Imax, Jmax, Kmax, Nvar];
     Flow(i).Dat = zeros(Imax,Jmax,Kmax,Nvar);
end
for i = 1:nblock
    Flow(i).Dat= reshape(fscanf(fid,'%e\n',Imax*Jmax*Kmax*Nvar),size(Flow(i).Dat));
end
fclose(fid);
end

