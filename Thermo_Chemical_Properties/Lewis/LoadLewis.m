% script to load lewis style data
% thermo.inp can be obtained from NASA CEA
% or RPA(http://www.propulsion-analysis.com/downloads/thermo/thermo.inp)
% or anyone of the following one
% NASA Technical Memorandum 4513
% NASA Technical Memorandum 4647
% NASA/TP-2001-210959/REV1
% NASA/TP-2002-211556
clear
clc
fclose all
fid = fopen('thermo.inp','r');

%% skip all line with ! at the beginning
while ~feof(fid)
    line = fgetl(fid);
    if ~strcmp(line(1),'!')
        break
    end
end
%% skip another two lines of comment
line = fgetl(fid);
i = 0;
data = repmat(struct('sp',[],'ref',[],'nT',0,'refcode',[],'formula',[],'phase',1,...
    'm0',0,'H0',0,'Trange',[],'Hdiff',[],'coeff',[]),100,1);
%% species
while ~feof(fid)
    % record 1
    line = fgetl(fid);
    if contains(line,'END REACTANTS') || contains(line,'END PRODUCTS')
        break
    end
    i = i+1;
    data(i).sp = strtrim(line(1:16));
    data(i).ref = strtrim(line(19:end));
    % record 2
    line = fgetl(fid);
    data(i).nT = sscanf(line(1:2),'%d');
    data(i).refcode = strtrim(line(4:9));
    chemformStr = strtrim(line(11:50));
    chemform = [];
    for j = 1:5
        elem = strtrim(chemformStr(1+(j-1)*8: 2+(j-1)*8));
        if numel(elem) > 0
            n = sscanf(chemformStr(3+(j-1)*8: 8+(j-1)*8), '%f');
             chemform.(elem) = n;
        end
    end
    data(i).formula = chemform;
    
    data(i).phase = sscanf(line(51:52),'%d');
    data(i).m0 = sscanf(line(53:65),'%f');
    data(i).H0 = sscanf(line(66:80),'%f');
    % record 3
    data(i).Trange = zeros(data(i).nT,2);
    data(i).Hdiff = zeros(data(i).nT,1);
    data(i).coeff = zeros(data(i).nT,9);
    coeffb = zeros(data(i).nT,2);
    
    for j = 1:data(i).nT
        line = fgetl(fid);
        Trange = sscanf(line(1:22),'%f');
        data(i).Trange(j,:) = reshape(Trange,1,[]);
        ncoeff = sscanf(line(23),'%d');
        try
            data(i).Hdiff(j) = sscanf(line(66:80),'%f');
        catch
            Hdiff = 0;
        end
        line = fgetl(fid);
        data(i).coeff(j,1:5) = sscanf(strrep(line,'D','e'),'%e');
        line = fgetl(fid);
        data(i).coeff(j,6:7) = sscanf(strrep(line(1:32),'D','e'),'%e');
        coeffb(j,:) = sscanf(strrep(line(49:80),'D','e'),'%e');
    end
    data(i).coeff(:,8:9) = coeffb;
end
nspecie = i;
fclose(fid);

unit = struct('mass','u','H0','J/mol', 'Cp', '1/R*J/mol');

fid=fopen('thermo.json','w');
fprintf(fid,jsonencode(data));
fclose(fid);

fid=fopen('unit.json','w');
fprintf(fid,jsonencode(unit));
fclose(fid);
