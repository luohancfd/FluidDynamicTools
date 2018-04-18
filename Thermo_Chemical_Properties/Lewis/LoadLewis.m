% script to load lewis style data 
% thermo.inp can be obtained from NASA CEA
% or RPA(http://www.propulsion-analysis.com/downloads/thermo/thermo.inp)
% or anyone of the following one
% NASA Technical Memorandum 4513
% NASA Technical Memorandum 4647
% NASA/TP—2001-210959/REV1
% NASA/TP—2002-211556
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
data = repmat(struct('sp',[],'ref',[],'nT',0,'refcode',[],'chemform',[],'type',1,...
    'mass',0,'heat',0,'Trange',[],'Hdiff',[],'coeff',[],'coeffb',[]),100,1);
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
    data(i).chemform = strtrim(line(11:50));
    data(i).type = sscanf(line(51:52),'%d');
    data(i).mass = sscanf(line(53:65),'%f');
    data(i).heat = sscanf(line(66:80),'%f');
    % record 3
    data(i).Trange = zeros(data(i).nT,2);
    data(i).Hdiff = zeros(data(i).nT,1);
    data(i).coeff = zeros(data(i).nT,7);
    data(i).coeffb = zeros(data(i).nT,2);

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
        data(i).coeffb(j,:) = sscanf(strrep(line(49:80),'D','e'),'%e');
    end
end
nspecie = i;
fclose(fid);
