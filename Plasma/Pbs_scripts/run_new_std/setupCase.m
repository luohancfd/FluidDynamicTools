%% =====================================
%% PUT YOUR INPUT FILE HERE
% module load comsol/5.4
% module load matlab/R2018a
%
% %https://www.comsol.com/forum/thread/9299/tip-batch-calculations-with-matlab-comsol-4-0a
function setupCase(caseName,timeLength,timeDiff,ifile,oldStdInput, newStd)
% case should all be named as caseName_stdXX.mph
% caseName: name of the case
% timeLength: number of RF cycles for new study
% timeDiff: number of RF cycles between output
% ifile: 0 use oldStdInput, if oldStdInput = 0, use caseName.mph
%        1 find last mph files matching pattern caseName_std*.mph
% oldStdInput: old std id, required for ifile == 0
% newStd: double check

%% start comsol server
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


%% backup base case
status = mkdir('backup');
caseToRun = [caseName,'.mph']; % we keep overwrite base file
caseToBack = sprintf('%s_%s.mph',caseName, datestr(now,'mm-dd-HH-MM'));
[status,msg,msgID] = copyfile(caseToRun, ['backup','/',caseToBack]);

if (nargin < 5)
    oldStdInput = 0;
end

if (ifile == 0)
    if (oldStdInput == 0)
        mphInput = [caseName, '.mph'];
    else
        mphInput = sprintf('%s_std%d.mph', caseName, oldStdInput);
    end
    oldStd = oldStdInput;
else
    [oldStd, mphInput] = findLastFile(caseName);
end


[mphOutput, newStdTemp] = addStudy(mphInput,timeLength,timeDiff);
if (newStdTemp ~= newStd)
    warn('Inconsistent STD');
    returnCode = -1;
else
    returnCode = newStdTemp;
end


fid = fopen('stdList.txt','w');
fprintf(fid, '%d\n', returnCode);
fclose(fid);


ModelUtil.disconnect;
