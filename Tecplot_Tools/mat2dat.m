function [ output_args ] = mat2dat(tdata,file)
%% export ascii format tecplot file
fid = fopen(file,'w');
if isfield(tdata,'auxdata')
    if isstruct(tdata.auxdata)
        for i = 1:length(tdata.auxdata)
            if isnumeric(tdata.auxdata(i).value)
                tdata.auxdata(i).value = sprintf('%.6G',tdata.auxdata(i).value);
            end
            fprintf(fid,'DATASETAUXDATA %s = "%s"\n',tdata.auxdata(i).name,tdata.auxdata(i).value);
        end
    end
end

fprintf(fid,'VARIABLES = ');
for i = 1:length(tdata.varnames)-1
    fprintf(fid,'"%s",',tdata.varnames{i});
end
i = i+1;
fprintf(fid,'"%s"\n',tdata.varnames{i});

lines = tdata.lines;
iz = isfield(lines, 'z');
iv = isfield(lines, 'v');
idimen = 1;
if iz
    idimen = idimen+1;
end

if iv
    maxv = 0;
    for ii = 1:length(lines)
        if ~isempty(lines(ii).v)
            [l,~]=size(lines(ii).v);
            if l > maxv
                maxv = l;
            end
        end
    end
    idimen = idimen + maxv;
end

for ii = 1:length(lines)
    lenx = length(lines(ii).x);
    if isfield(lines(ii),'zonename')
        if ~isempty(lines(ii).zonename)
            fprintf(fid,'ZONE I = %d T="%s" ',lenx,lines(ii).zonename);
        else
            fprintf(fid,'ZONE I = %d T="ZONE %d" ',lenx,ii);
        end
    else
        fprintf(fid,'ZONE I = %d T="ZONE %d" ',lenx,ii);
    end

    if isfield(lines(ii),'varloc')
        if lines(ii).varloc == 0
            ivarloc = 0;
        else
            ivarloc = 1;
        end
    else
        ivarloc = 0;
    end

    datax = lines(ii).x;
    if ivarloc == 0
        fprintf(fid,'DATAPACKING=POINT\n');
        lendata = lenx;
    else
        fprintf(fid,'DATAPACKING=BLOCK, VARLOCATION=([2-%d]=CELLCENTERED)\n',idimen);
        lendata = lenx-1;
    end

    % zone aux data
    if isfield(lines(ii),'auxname')
        if ~isempty(lines(ii).auxname)
            aux_values = lines(ii).auxval;
            if iscell(aux_values)
                for iiii = 1:length(aux_values)
                    if isnumeric(aux_values{iiii})
                        aux_values{iiii} = sprintf('%.6G',aux_values{iiii});
                    end
                end
            else
                if isnumeric(aux_values)
                    aux_values = {sprintf('%.6G',aux_values)};
                end
            end
            for jj = 1:length(lines(ii).auxname)
                fprintf(fid,'AUXDATA %s = "%s"\n',lines(ii).auxname{jj},aux_values{jj});
            end
        end
    end
    data = zeros(lendata,idimen-1);
    data(:,1) = reshape(lines(ii).y,lendata,1);
    if iz
        if ~isempty(lines(ii).z)
            data(:,2) = reshape(lines(ii).z,lendata,1);
        end
    end
    if iv
        if ~isempty(lines(ii).v)
            [l,~] = size(lines(ii).v);
            data(:,3:3+l-1) = lines(ii).v';
        end
    end

    if ivarloc == 0
        for i = 1:lenx
            fprintf(fid,'%.6G  ',datax(i));
            for j = 1:idimen
                fprintf(fid,'%.6G  ',data(i,j));
            end
            fprintf(fid,'\n');
        end
    else
        for i = 1:lenx-1
            fprintf(fid,'%.6G  ',datax(i));
        end
        fprintf(fid,'%.6G\n',datax(lenx));
        for j = 1:idimen -1
            for i = 1:lendata-1
                fprintf(fid,'%.6G  ',data(i,j));
            end
            fprintf(fid,'%.6G\n',data(lendata,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

