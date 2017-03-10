function [nblock,Nvar] = ExportFUNCTION( Flow,filename )
%% export correct plot3d multiblock format,
% Flow must be 3D  i.e. length(size) = 3
% doesn't necessay to be multiblock grid
% keep in mind, X, Y , Z rows is for i, colum is for j, page is for k
nblock = length(Flow);
Nvar = Flow(1).Size(end);
fid = fopen(filename,'w');
fprintf(fid, '%d\n',nblock);
for i=1:nblock
    fprintf(fid, '%d %d %d %d\n',[Flow(i).Size]);
end
for i = 1:nblock
    fprintf(fid,'%e\n',Flow(i).Dat);
end
fclose(fid);
end


