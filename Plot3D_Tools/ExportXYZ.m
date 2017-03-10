function [nblock] = ExportXYZ( Grid,filename )
%% export correct plot3d multiblock format,
% Grid must be 3D grid i.e. length(size) = 3
% doesn't necessay to be multiblock grid
% keep in mind, X, Y , Z rows is for i, colum is for j, page is for k
nblock = length(Grid);
fid = fopen(filename,'w');
fprintf(fid, '%d\n',nblock);

for i=1:nblock
    fprintf(fid, '%d ',Grid(i).Size);
    fprintf(fid,'\n');
end

for Nloop = 1:nblock
    fprintf(fid,'%e\n',Grid(Nloop).X);
    fprintf(fid,'%e\n',Grid(Nloop).Y);
    fprintf(fid,'%e\n',Grid(Nloop).Z);
    %  The above have the same effect as the following one.
    %    [imax, jmax, kmax ] = size(Grid(Nloop).X);
    %     for k=1:kmax
    %         for j = 1:jmax
    %             for i = 1:imax
    %                     fprintf(fid,'%e\n',Grid(Nloop).X(i,j,k));
    %             end
    %         end
    %     end
end
fclose(fid);
end


