function [ output_args ] = mat2dat(tdata,file)
%% export ascii format tecplot file
fid = fopen(file,'w');
fprintf(fid,'VARIABLES = ');
for i = 1:length(tdata.varnames)-1
    fprintf(fid,'"%s",',tdata.varnames{i});
end
i = i+1;
fprintf(fid,'"%s"\n',tdata.varnames{i});

lines = tdata.lines;
iz = isfield(lines, 'z');
iv = isfield(lines, 'v');
idimen = 2;
if iz
    idimen = idimen+1;
end

if iv
    maxv = 0;
    for z = 1:length(lines)
        if ~isempty(lines(z).v)
            [l,~]=size(lines(z).v);
            if l > maxv
                maxv = l;
            end
        end
    end
    idimen = idimen + maxv;
end
    
for z = 1:length(lines)
    len = length(lines(z).x);
    fprintf(fid,'ZONE I = %d T="%s"\n',len,lines(z).zonename);
    data = zeros(len,idimen);
    data(:,1) = reshape(lines(z).x,len,1);
    data(:,2) = reshape(lines(z).y,len,1);
    if iz
        if ~isempty(lines(z).y)
            data(:,3) = reshape(lines(z).z,len,1);
        end
    end
    if iv
        if ~isempty(lines(z).v)
            [l,~] = size(lines(z).v);
            data(:,4:4+l-1) = lines(z).v';
        end
    end
    for i = 1:len
        for j = 1:idimen - 1
            fprintf(fid,'%.6G  ',data(i,j));
        end
        j = j+1;
        fprintf(fid,'%.6G  \n',data(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

