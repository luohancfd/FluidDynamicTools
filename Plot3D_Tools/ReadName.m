function [Name]  = ReadName(filename)
%% Read 
Name = cell(100,1);
% filename = '0/plot/sphere.t0216.nam';
fid = fopen(filename,'r');
nfield = 0;
cline = fgetl(fid);
while ischar(cline)
    nfield = nfield + 1;
    Name{nfield} = strtrim(cline);
    cline = fgetl(fid);
end
fclose(fid);
Name = Name(1:nfield,:);
end