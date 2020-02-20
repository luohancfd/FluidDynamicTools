%% start comsol server
caseName = '40V';
serverLog = sprintf('server_%s.log', caseName);
portFile = sprintf('port_%s.log', caseName);
COMSOL_LIC = '';
cmdLine = sprintf('nohup comsol server %s -portfile %s </dev/null > %s 2>&1 &',COMSOL_LIC, portFile, serverLog);
[status,cmdout] = system(cmdLine);
duration = 20;
java.lang.Thread.sleep(duration*1000)
fid = fopen(portFile,'r');
line = fgetl(fid);
fclose(fid);
% delete(serverLog);
% fprintf(portFile);
port = sscanf(line,'%d');
fprintf('Connect to port %d\n',port)
addpath(genpath('/apps/cent7/comsol/5.4/build295/mli'));
mphstart(port)
import com.comsol.model.util.*

%% your code here
mphInput = '160MHz_std9.mph';

model = mphopen(mphInput);
modelInfo = mphmodel(model);
modelStudy = mphmodel(model.study);
modelDsets = fieldnames(mphmodel(model.result.dataset));
%% get info of last study
studyLabels = fieldnames(modelStudy);
studyIDs = zeros(1,length(studyLabels));
for i = 1:length(studyLabels)
    studyIDs(i) = sscanf(studyLabels{i}, 'std%d',1);
end
lastStudyID = max(studyIDs);  % ID of last study
lastStudyTimeDefine = sprintf('t_std%d',lastStudyID);
lastStudyLabel = sprintf('std%d',lastStudyID);
lastStudy = model.study(lastStudyLabel);
lastStudySolver = lastStudy.feature('time');
lastStudySol = char(lastStudy.getSolverSequences('all'));
% find dset => solution
dsetStd = ones(length(modelDsets),1)*(-1);
maxDsetId = 0;
for i = 1:length(modelDsets)
    if contains(modelDsets{i},'dset')
        dsetId = sscanf(modelDsets{i},'dset%d',1);
        dsetTag = modelDsets{i};
        solTag = model.result.dataset(dsetTag).getString('solution');
        dsetStd(dsetId) = sscanf(char(solTag),'sol%d',1);
        maxDsetId = dsetId;
    end
end
dsetStd = dsetStd(1:maxDsetId);
%% get info of last average
lastAverageNum = 0;
numericalTags = fieldnames(mphmodel(model.result.numerical));
for i = 1:length(numericalTags)
    name = char(model.result.numerical(numericalTags{i}).getDisplayString());
    if contains(name,'Line Average')
        lastAverageNum = lastAverageNum+1;
    end
end
lastAverage.Id = lastAverageNum;
lastAverage.Tag = sprintf('av%d',lastAverageNum);
lastAverage.Dset = char(model.result.numerical(lastAverage.Tag ).getString('data'));
lastAverage.DsetId = sscanf(lastAverage.Dset, 'dset%d',1);
tableNames = fieldnames(mphmodel(model.result.table));
lastTableID = sscanf(tableNames{end},'tbl%d');
%% create new average
newAverageId = lastAverage.Id;
newAverageTable = lastTableID;
for i = lastAverage.DsetId+1:maxDsetId
    newAverageId = newAverageId + 1;
    newAverageTable = newAverageTable + 1;
    
    newAverageLabel = sprintf('av%d',newAverageId);
    newAverageDset = sprintf('dset%d',i);
    newAverageTableTag = sprintf('tbl%d',newAverageTable)
    
    model.result.numerical.duplicate(newAverageLabel, lastAverage.Tag);
    model.result.numerical(newAverageLabel).set('data', newAverageDset);
    model.result.table.create(newAverageTableTag, 'Table');
    model.result.numerical(newAverageLabel).set('table', newAverageTableTag);
    model.result.numerical(newAverageLabel).set('innerinput', 'all');
    model.result.numerical(newAverageLabel).setResult;
end

mphsave(model,mphInput);

