function [lastStd, lastMph] = findLastFile(caseName)
% find files matching caseName_std*.mph

fileNamePattern = [caseName,'_std*.mph'];
mphFilesTemp = ls(fileNamePattern);
nFile = 0;
if isunix
    mphFilesTemp = strsplit(mphFilesTemp);
    for i = 1:length(mphFilesTemp)
        if ~isempty(mphFilesTemp{i})
            nFile = nFile + 1;
        else
            break
        end
    end
else
    [nFile,~] = size(mphFilesTemp);
end

if (nFile == 0)
    lastStd = -1;
    lastMph = '';
else
    mphFiles = cell(1,nFile);
    if isunix
        for i = 1:nFile
            mphFiles{i} = mphFilesTemp{i};
        end
    else
        for i = 1:nFile
            mphFiles{i} = mphFilesTemp(i,:);
        end
    end

    stdID = zeros(1,length(mphFiles));
    fileNamePattern = [caseName,'_std%d.mph'];
    for i = 1:length(mphFiles)
        z = sscanf(mphFiles{i},fileNamePattern);
        if isempty(z)
            z = 0;
        end
        stdID(i) = z;
    end
    [lastStd,index] = max(stdID);
    lastMph = mphFiles{index};
end