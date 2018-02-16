function [zone,VARlist] = tec2mat(fname)
%% load tecplot ASCII file and write result to zone and VARList
fid = fopen(fname,'r');
IVARlist = 0;
NVAR = 0;
%% read file untill we get a line of variables name
while ~feof(fid)
    line = fgetl(fid);
    line = strtrim(line);
    if ~isempty(line)
        if (~strcmp(line(1),'#'))
            if contains(line,'VARIABLES')
                VARlist = regexp(line,'[''"](.*?)[''"]','tokens');
                if ~isempty(VARlist)
                    for i = 1:length(VARlist)
                        VARlist{i} = VARlist{i}{1};
                    end
                    NVAR = length(VARlist);
                    IVARlist = 1;
                    break
                else
                    error("variable list is not correct")
                end
            end
        end
    end
end

nzone = 0;
zone = repmat(struct('data',[],'passivevarlist',[],'title',[]),10,1);
if IVARlist == 0
    error("there is no variable list")
else
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if contains(line,'ZONE')
            idata = regexp(line,'I\s*=\s*(\d+)','tokens');
            idata = str2num(idata{1}{1});
            jdata = regexp(line,'J\s*=\s*(\d+)','tokens');
            if isempty(jdata)
                jdata = 1;
            else
                jdata = str2num(jdata{1}{1});
            end
            ndata = idata*jdata;
            
            title = regexp(line,'T\s*=\s*["''](.*)[''"]','tokens');
            if ~isempty(title)
                title = title{1}{1};
            %    title = strrep(title,'"','');  title = strrep(title,'''','');
            else
                title = regexp(line,'T\s*=\s*(.*)','tokens');
                if ~isempty(title)
                    title = title{1}{1};
                else
                    title ='no title'
                end
            end
            
            % load zone
            mcolum = NVAR;   
            pos = ftell(fid);  %line above PASSIVEVARLIST or data
            line = fgetl(fid);     line = strtrim(line);     
            PASSIVEVARLIST = [];
            if contains(line,"PASSIVEVARLIST")
                list = regexp(line,'[\d-]','match');
                for i =1:length(list)
                    if isstrprop(list{i},'digit')
                        ii = str2num(list{i});
                        PASSIVEVARLIST = [PASSIVEVARLIST ii];
                        if i < length(list)
                            if strcmp(list{i+1},'-')
                                jj = str2num(list{i+2});
                                for j = ii+1:jj-1
                                    PASSIVEVARLIST = [PASSIVEVARLIST j];
                                end
                            end
                        end
                    end
                end
                mcolum = mcolum - length(PASSIVEVARLIST);
            else
                fseek(fid,pos,-1);
            end
            
            % let's load the data
            data = fscanf(fid,'%f',ndata*mcolum);
            data = reshape(data,mcolum,ndata)';
            
            nzone = nzone + 1;
            zone(nzone).data = data;
            zone(nzone).passivevarlist = PASSIVEVARLIST;
            zone(nzone).title = title;
        end
    end
end            

zone = zone(1:nzone);
                                
fclose(fid);
end

