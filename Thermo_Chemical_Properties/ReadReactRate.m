function [sp,Arrates,i_reac,i_source,i_sink,mass,MassDB] = ReadReactRate(fname)
%% Read DPLR or US3D style reaction rates
% fname = 'air-11sp-park90-correct2.dat';
fid = fopen(fname,'r');
Arrates = repmat(struct('type',[],'r1','','r2','','r3','','p1','','p2','','p3','',...
    'Cfm',[],'eta',[],'Ea',[],'irxon',[],...
    'txf',[],'txb',[],'name',[]),20,1);

igeneral = zeros(20,11);
ngeneral = 0;
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line)
        continue
    end
    if strcmp(line(1),'!')
        continue
    end
    if strcmp(line,'[GAS_SPECIES]')
        sp = strtrim(fgetl(fid));
        %      sp = strrep(sp,'+','p');
        sp = strsplit(sp,',');
        nsp = length(sp);
        igeneral = zeros(20,nsp);
    end
    if strcmp(line,'[GAS_REACTIONS]')
        while true
            line = strtrim(fgetl(fid));
            if isempty(line)
                continue
            end
            if ~strcmp(line(1),'!')
                break
            end
        end
        nreac = 0;
        while ~contains(line,'[/GAS_REACTIONS]') && ~strcmp(line(1),'!')
            reac_line = strsplit(line,'|');
            reac_line = cellfun(@strtrim,reac_line,'UniformOutput',false);
            % reaction type
            switch reac_line{1}
                case 'diss'
                    type = 10;
                case 'cs_diss'
                    type = 11;
                case 'recomb'
                    type = 20;
                case 'cs_recomb'
                    type = 21;
                case 'exch'
                    type = 3;
                case 'ion'
                    type = 4;
                case '23'
                    type = 5;
                case '32'
                    type = 6;
                otherwise
                    type = 7;
            end
            
            % species
            reac_sp = cellfun(@strtrim,strsplit(reac_line{2},'<=>'),'UniformOutput',false);
            sp0 = cellfun(@strtrim,strsplit(reac_sp{1},'+'),'UniformOutput',false);
            sp1 = cellfun(@strtrim,strsplit(reac_sp{2},'+'),'UniformOutput',false);
            sp0 = cellfun(@(x) strrep(x,'p','+'),sp0,'UniformOutput',false);
            sp1 = cellfun(@(x) strrep(x,'p','+'),sp1,'UniformOutput',false);
            
            % Arrhenius parameter
            Arpar = cellfun(@strtrim,strsplit(reac_line{3},','),'UniformOutput',false);
            Arpar = cellfun(@(x)strrep(x,'d','e'),Arpar,'UniformOutput',false);
            Arpar = cell2mat(cellfun(@str2double,Arpar,'UniformOutput',false));
            if length(Arpar) == 5
                Arpar = [0 Arpar];
            end
            % Put back to struct
            if (strcmp(sp0{end},'M') && mod(type,10)==0 ) %general reaction
                ngeneral = ngeneral+1;
                if type/10 ==1 %general dissociation
                    for i = 1:nsp
                        nreac = nreac+1;
                        igeneral(ngeneral,i) = nreac;
                        Arrates(nreac).type = type;
                        Arrates(nreac).r1 = sp0{1};
                        Arrates(nreac).r2 = sp{i};
                        Arrates(nreac).p1 = sp1{1};
                        Arrates(nreac).p2 = sp1{2};
                        Arrates(nreac).p3 = sp{i};
                        Arrates(nreac).Cfm = Arpar(1); Arrates(nreac).eta = Arpar(2);
                        Arrates(nreac).Ea  = Arpar(3); Arrates(nreac).irxon = Arpar(4);
                        Arrates(nreac).txf = Arpar(5); Arrates(nreac).txb = Arpar(6);
                        Arrates(nreac).name = sprintf('%3s+%3s->%3s+%3s+%3s',sp0{1},sp{i},sp1{1},sp1{2},sp{i});
                    end
                elseif type/10 == 2 %general recombination
                    for i = 1:nsp
                        nreac = nreac+1;
                        igeneral(ngeneral,i) = nreac;
                        Arrates(nreac).type = type;
                        Arrates(nreac).r1 = sp0{1};
                        Arrates(nreac).r2 = sp0{2};
                        Arrates(nreac).r3 = sp{i};
                        Arrates(nreac).p1 = sp1{1};
                        Arrates(nreac).p2 = sp{i};
                        Arrates(nreac).Cfm = Arpar(1); Arrates(nreac).eta = Arpar(2);
                        Arrates(nreac).Ea  = Arpar(3); Arrates(nreac).irxon = Arpar(4);
                        Arrates(nreac).txf = Arpar(5); Arrates(nreac).txb = Arpar(6);
                        Arrates(nreac).name = sprintf('%3s+%3s+%3s->%3s+%3s',sp0{1},sp0{2},sp{i},sp1{1},sp{i});
                        
                    end
                end
            else
                if length(sp0) == 2
                    sp0{3} = '';
                end
                if length(sp1) == 2
                    sp1{3} = '';
                end
                if type == 11
                    for i = 1:nreac
                        if strcmp(sp0{1},Arrates(i).r1) && strcmp(sp0{2},Arrates(i).r2) && Arrates(i).type == 10
                            Arrates(i).type = type;
                            Arrates(i).r1 = sp0{1};
                            Arrates(i).r2 = sp0{2};
                            Arrates(i).p1 = sp1{1};
                            Arrates(i).p2 = sp1{2};
                            Arrates(i).p3 = sp1{3};
                            
                            Arrates(i).Cfm = Arpar(1); Arrates(i).eta = Arpar(2);
                            Arrates(i).Ea  = Arpar(3); Arrates(i).irxon = Arpar(4);
                            Arrates(i).txf = Arpar(5); Arrates(i).txb = Arpar(6);
                            Arrates(nreac).name = sprintf('%3s+%3s->%3s+%3s+%3s',sp0{1},sp0{2},sp1{1},sp1{2},sp1{3});
                            
                            break
                        end
                    end
                elseif type == 21
                    for i = 1:nreac
                        if strcmp(sp0{1},Arrates(i).r1) && strcmp(sp0{2},Arrates(i).r2)&& Arrates(i).type == 20
                            Arrates(i).type = type;
                            Arrates(i).r1 = sp0{1};
                            Arrates(i).r2 = sp0{2};
                            Arrates(i).r3 = sp0{3};
                            Arrates(i).p1 = sp1{1};
                            Arrates(i).p2 = sp1{2};
                            
                            Arrates(i).Cfm = Arpar(1); Arrates(i).eta = Arpar(2);
                            Arrates(i).Ea  = Arpar(3); Arrates(i).irxon = Arpar(4);
                            Arrates(i).txf = Arpar(5); Arrates(i).txb = Arpar(6);
                            Arrates(nreac).name = sprintf('%3s+%3s+%3s->%3s+%3s',sp0{1},sp0{2},sp0{3},sp1{1},sp1{2});
                            break
                        end
                    end
                else
                    nreac = nreac+1;
                    Arrates(nreac).type = type;
                    Arrates(nreac).r1 = sp0{1};
                    Arrates(nreac).r2 = sp0{2};
                    Arrates(nreac).r3 = sp0{3};
                    Arrates(nreac).p1 = sp1{1};
                    Arrates(nreac).p2 = sp1{2};
                    Arrates(nreac).p3 = sp1{3};
                    
                    Arrates(nreac).Cfm = Arpar(1); Arrates(nreac).eta = Arpar(2);
                    Arrates(nreac).Ea  = Arpar(3); Arrates(nreac).irxon = Arpar(4);
                    Arrates(nreac).txf = Arpar(5); Arrates(nreac).txb = Arpar(6);
                    
                    if isempty(sp0{3}) && isempty(sp1{3})
                        Arrates(nreac).name = sprintf('%3s+%3s->%3s+%3s',sp0{1},sp0{2},sp1{1},sp1{2});
                    elseif isempty(sp0{3})
                        Arrates(nreac).name = sprintf('%3s+%3s->%3s+%3s+%3s',sp0{1},sp0{2},sp1{1},sp1{2},sp1{3});
                    elseif isempty(sp1{3})
                        Arrates(nreac).name = sprintf('%3s+%3s+%3s->%3s+%3s',sp0{1},sp0{2},sp0{3},sp1{1},sp1{2});
                    else
                        Arrates(nreac).name = sprintf('%3s+%3s+%3s->%3s+%3s+%3s',sp0{1},sp0{2},sp0{3},sp1{1},sp1{2},sp1{3});
                    end                    
                end
                
            end
            line = strtrim(fgetl(fid));
        end
        Arrates = Arrates(1:nreac);
        
        % clean the table
        
    end
    
    if strcmp(line,'[GAS_REACTION_PARAMETERS]')
        while true
            line = strtrim(fgetl(fid));
            if isempty(line)
                continue
            end
            if ~strcmp(line(1),'!')
                break
            end
        end
        kcount = 0;
        while ~contains(line,'[/GAS_REACTION_PARAMETERS]') && ~strcmp(line(1),'!')
            w= regexp(line,'(\d+)(?:\~(\d+))?\:(\d\.\d+[deE][+-]\d+)','match');
            if ~isempty(w)
                kcount = kcount+1;
                for i = 1:length(w)
                    Apar = strsplit(w{i},':');
                    Cfm = str2double(strrep(Apar{2},'d','e'));
                    sp0 = cellfun(@str2double,strsplit(Apar{1},'~'),'UniformOutput',false);
                    sp0 = cell2mat(sp0);
                    if length(sp0) == 2
                        for j = sp0(1):sp0(2)
                            index =igeneral(kcount,j);
                            Arrates(index).Cfm = Cfm;
                        end
                    else
                        index = igeneral(kcount,sp0);
                        Arrates(index).Cfm = Cfm;
                    end
                end
            end
            line = strtrim(fgetl(fid));
        end
    end
end
fclose(fid);

%% assign numerical id
cellfind=@(x,y) find(cellfun(@(z)strcmp(z,y),x));
pair_reac = cell(nsp,nsp);
i_reac = cell(nsp,1);
i_source = cell(nsp,1);
i_sink = cell(nsp,1);
for i = 1:nreac
    r3n = zeros(1,3);
    for j = 1:3
        w = cellfind(sp,Arrates(i).(['r',num2str(j)]));
        if ~isempty(w)
            r3n(j) = w;
        end
    end
    Arrates(i).reac = r3n;
    r3n = sort(r3n(1:2));  %!ASSUME ONLY BINARY COLLISION, NO RECOMBINATION
    
    % count pair
    pair_reac{r3n(1),r3n(2)} = [pair_reac{r3n(1),r3n(2)} i];
    i_reac{r3n(1)} = [i_reac{r3n(1)}, i];
    i_reac{r3n(2)} = [i_reac{r3n(2)}, i];
    
    p3n = zeros(1,3);
    for j = 1:3
        w = cellfind(sp,Arrates(i).(['p',num2str(j)]));
        if ~isempty(w)
            p3n(j) = w;
        end
    end
    Arrates(i).prod = p3n;
    
    % calculate the net result of this reaction
    npre = zeros(1,nsp); npost = npre;
    for ii = 1:3
        if ii <= length(r3n)
            if r3n(ii) ~= 0
                npre(r3n(ii)) = npre(r3n(ii))+1;
            end
        end
        if ii<= length(p3n)
            if p3n(ii) ~= 0
                npost(p3n(ii)) = npost(p3n(ii))+1;
            end
        end
    end
    net_change = npost - npre;
    
    for ii = 1:length(net_change)
        if net_change(ii)>0
            i_source{ii} = [i_source{ii}, [i;net_change(ii)]];
        elseif net_change(ii) < 0
            i_sink{ii} = [  i_sink{ii}, [i;-net_change(ii)]];
        end
    end
end
%% clean name
for i = 1:nreac
    Arrates(i).name = regexprep(Arrates(i).name,'\+(?=[-\+])','<sup>+</sup>');
    Arrates(i).name = regexprep(Arrates(i).name,'\+$','<sup>+</sup>');
    Arrates(i).name = strrep(Arrates(i).name,'2','<sub>2</sub>');
    Arrates(i).name = strrep(Arrates(i).name,' ','');
end
    
%% Mass database
MassDB = {  'Ar',     39.94400e+0,  0.0000000000e0, 1.5e0,   0  ,  1;
    'Ar+'    39.94345e+0,  3.8068120000e7, 1.5e0,   1  ,  1;
    'C',      12.01100e+0,  5.9211889000e7, 1.5e0,   0  ,  1;
    'C+'     12.01045e+0,  1.4967366000e8, 1.5e0,   1  ,  1;
    'C2',     24.02200e+0,  3.4234785000e7, 2.5e0,   0  ,  2;
    'C2H',    25.03000e+0,  2.2572910000e7, 2.5e0,   0  ,  2;
    'C2H2',   26.03800e+0,  8.7548580000e6, 2.5e0,   0  ,  2;
    'C3',     36.03300e+0,  2.3062193000e7, 2.5e0,   0  ,  2;
    'CF',     31.00940e+0,  9.9617550000e6, 2.5e0,   0  ,  2;
    'CF2',    50.00780e+0, -2.5187600000e6, 3.0e0,   0  ,  2;
    'CF3',    69.00620e+0, -7.1992350000e6, 3.0e0,   0  ,  2;
    'CF4',    88.00460e+0, -1.0258770000e7, 3.0e0,   0  ,  2;
    'CH',     13.01900e+0,  4.5627850000e7, 2.5e0,   0  ,  2;
    'CH2',    14.02700e+0,  2.7803520000e7, 3.0e0,   0  ,  2;
    'CH3',    15.03500e+0,  9.9559030000e6, 3.0e0,   0  ,  2;
    'CH4',    16.04300e+0, -4.1532130000e6, 3.0e0,   0  ,  2;
    'Cl',     35.45300e+0,  3.3740400000e6, 1.5e0,   0  ,  1;
    'Cl2',    70.90600e+0,  0.0000000000e0, 2.5e0,   0  ,  2;
    'CN',     26.01900e+0,  1.6795420000e7, 2.5e0,   0  ,  2;
    'CN+'    26.01845e+0,  6.8835800000e7, 2.5e0,   1  ,  2;
    'CO',     28.01100e+0, -4.0630824000e6, 2.5e0,   0  ,  2;
    'CO+'    28.01045e+0,  4.4200904000e7, 2.5e0,   1  ,  2;
    'CO2',    44.01100e+0, -8.9328800000e6, 2.5e0,   0  ,  2;
    'F',      18.99840e+0,  4.0423050000e6, 1.5e0,   0  ,  1;
    'F2',     37.99680e+0,  0.0000000000e0, 2.5e0,   0  ,  2;
    'H',       1.00800e+0,  2.1432040000e8, 1.5e0,   0  ,  1;
    'H+'      1.00745e+0,  1.5167840000e9, 1.5e0,   1  ,  1;
    'H2',      2.01600e+0,  0.0000000000e0, 2.5e0,   0  ,  2;
    'H2+'     2.01545e+0,  7.3847530000e8, 2.5e0,   1  ,  2;
    'H2O',    18.01600e+0, -1.3261710000e7, 3.0e0,   0  ,  2;
    'HCl',    36.46100e+0, -2.5266800000e6, 2.5e0,   0  ,  2;
    'HCN',    27.02700e+0,  4.8982130000e6, 2.5e0,   0  ,  2;
    'He',      4.00300e+0,  0.0000000000e0, 1.5e0,   0  ,  1;
    'He+'     4.00245e+0,  5.9271800000e8, 1.5e0,   1  ,  1;
    'N',      14.00800e+0,  3.3621610000e7, 1.5e0,   0  ,  1;
    'Ne',     20.17900e+0,  0.0000000000e0, 1.5e0,   0  ,  1;
    'N+'     14.00745e+0,  1.3400000000e8, 1.5e0,   1  ,  1;
    'N2',     28.01600e+0,  0.0000000000e0, 2.5e0,   0  ,  2;
    'N2+'    28.01545e+0,  5.3700000000e7, 2.5e0,   1  ,  2;
    'NCO',    42.01900e+0,  4.2124000000e6, 2.5e0,   0  ,  2;
    'NH',     15.01600e+0,  2.3867900000e7, 2.5e0,   0  ,  2;
    'NH+'    15.01545e+0,  1.1050000000e8, 2.5e0,   1  ,  2;
    'NH2',    16.02400e+0,  1.2036000000e7, 3.0e0,   0  ,  2;
    'NH3',    17.03200e+0, -2.2866370000e6, 3.0e0,   0  ,  2;
    'NO',     30.00800e+0,  2.9961230000e6, 2.5e0,   0  ,  2;
    'NO+'    30.00745e+0,  3.2834800000e7, 2.5e0,   1  ,  2;
    'NO2',    46.00800e+0,  8.0420800000e5, 3.0e0,   0  ,  2;
    'O',      16.00000e+0,  1.5420000000e7, 1.5e0,   0  ,  1;
    'O+'     15.99945e+0,  9.7560000000e7, 1.5e0,   1  ,  1;
    'O2',     32.00000e+0,  0.0000000000e0, 2.5e0,   0  ,  2;
    'O2+'    31.99945e+0,  3.6370000000e7, 2.5e0,   1  ,  2;
    'OH',     17.00800e+0,  2.2995060000e6, 2.5e0,   0  ,  2;
    'Si',     28.08550e+0,  1.5868220000e7, 1.5e0,   0  ,  1;
    'SiO',    44.08550e+0, -2.2683200000e6, 2.5e0,   0  ,  2;
    'e',       0.00055e+0,  0.0000000000e0, 1.5e0,  -1  ,  1;
    };
% [m,n] = size(MassDB);
mass = zeros(1,nsp);
for i = 1:length(sp)
    k = cellfind(MassDB(:,1),sp{i});
    mass(i) = MassDB{k,2};
end


end



