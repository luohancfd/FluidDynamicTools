function [mphOutput,newStudyID] = addStudy(mphInput,timeLength,timeDiff)
cellFind = @(x,y) find(contains(x,y));
%% set input
% mphInput = 'Vpp_122_std2.mph';
% timeLength = 2500; % unit: period
% timeDiff = 10; % unit: period
timeDifference = sprintf('+%d*period',timeLength);
% append to new time definition, new time end is t_previous+timeDifference

%% load model
model = mphopen(mphInput);

%% get model info
%mphnavigator
modelInfo = mphmodel(model);
modelStudy = mphmodel(model.study);

%% get info of last study
studyLabels = fieldnames(modelStudy);
studyIDs = zeros(1,length(studyLabels));
for i = 1:length(studyLabels)
    studyIDs(i) = sscanf(studyLabels{i}, 'std%d',1);
end
lastStudyID = max(studyIDs);  % ID of last study
lastStudyTimeDefine = sprintf('t_std%d',lastStudyID);

%% setup new parameter
% summarize parameters
paramsGroupName = fieldnames(mphmodel(model.param.group()));
paramsGroup = repmat(struct('id',[],'label',[],'variables',[]),length(paramsGroupName),1);
[paramsGroup.id] =  paramsGroupName{:};
for i = 1:length(paramsGroup)
    group = model.param.group(paramsGroup(i).id);
    paramsGroup(i).label = group.label();
    paramsGroup(i).variables = mphgetexpressions(group);
end

% find where t_stdX is defined
paramsTimeGroupIndex = -1;
for i = 1:length(paramsGroup)
    if ~isempty(paramsGroup(i).variables)
        index = cellFind(paramsGroup(i).variables(:,1),lastStudyTimeDefine);
        if ~isempty(index)
            paramsTimeGroupIndex = i;
            break
        end
    end
end
if paramsTimeGroupIndex == -1
    error('No t_stdX found');
end
paramsTimeGroup = model.param.group(paramsGroup(paramsTimeGroupIndex).id);

% add new parameter to the group for new run
newStudyID = lastStudyID + 1;
newStudyTimeDefine = sprintf('t_std%d',newStudyID);
newStudyTimeExperssion = sprintf('%s%s',lastStudyTimeDefine,timeDifference);
paramsTimeGroup.set(newStudyTimeDefine, newStudyTimeExperssion);
paramsGroup(paramsTimeGroupIndex).variables = mphgetexpressions(model.param.group(paramsGroup(paramsTimeGroupIndex).id));

%% check last study
lastStudyLabel = sprintf('std%d',lastStudyID);
lastStudy = model.study(lastStudyLabel);
lastStudySolver = lastStudy.feature('time');
%lastStudyTlist = lastStudySolver.getDoubleArray('tlist');


%% create new study
newStudyLabel = sprintf('std%d',newStudyID);
model.study().create(newStudyLabel);
newStudy = model.study(newStudyLabel);
newStudy.setGenPlots(false);
newStudy.create('time','Transient');
newStudy.feature('time').set('rtol',0.0001);
newStudy.feature('time').set('initstudy',lastStudyLabel);
newStudy.feature('time').set('initmethod','sol');
newStudy.feature('time').set('initsoluse','current');
newStudy.feature('time').set('solnum','last');
newStudy.feature('time').set('useinitsol','on');
newStudy.feature('time').set('usesol','on');

t_start = mphevaluate(model,lastStudyTimeDefine);
period = mphevaluate(model,'period');
t_end = mphevaluate(model,newStudyTimeDefine);
tlist = t_start:timeDiff*period:t_end;
newStudy.feature('time').set('tlist',tlist);

% generate default solvers
newStudy.createAutoSequences('all');
newSolName = char(newStudy.getSolverSequences('all'));
newSol = model.sol(newSolName);

newTimeSolver = newSol.feature('t1');
fullCoupleSolver = newTimeSolver.feature('fc1');
fullCoupleSolver.set('jtech','onevery');
fullCoupleSolver.set('maxiter',25);

model.sol(newSolName).runFromTo('st1', 'v1');

mphOutput = regexprep(mphInput, '(.*_)std\d+(\.mph)',['$1',newStudyLabel,'$2']);

mphsave(model,mphOutput);

end