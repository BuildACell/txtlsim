%
function mergedExpFile = mergeWholeExpFile(expFile,wellsToCompare)

% check well consistence
% if only one row of well were given, then all experiment has the same order
if size(wellsToCompare,1) == 1
    wellsToCompare =  repmat(wellsToCompare,size(expFile,1),1);
end
% if only one row of expFile were given, then all wells are in the same
% expFile
if size(expFile,1) == 1
    expFile =  repmat(expFile,size(wellsToCompare,1),1);
end

%! TODO zoltuz 9/7/13 check the num of channels are the same in each cases
numOfWells = size(wellsToCompare,2);
numOfChannels = size(expFile{1}.channels,1);

% rawData
rawData = cellfun(@(x,y) x.Data(:,y,:),expFile,num2cell(wellsToCompare,2),'UniformOutput',false);
[rawData_mean rawData_std] = getStdMean(rawData);
% no background
noBackgroundData = cellfun(@(x,y) x.noBg(:,y,:),expFile,num2cell(wellsToCompare,2),'UniformOutput',false);
[noBackgroundData_mean noBackgroundData_std] = getStdMean(noBackgroundData);
% rate
rates = cellfun(@(x,y) x.rate(y,:),expFile,num2cell(wellsToCompare,2),'UniformOutput',false);
ratesProcessed = cellfun(@(x) reshape(cell2mat(x),1,numOfWells,numOfChannels),rates,'UniformOutput',false);
[rates_mean rates_std] = getStdMean(ratesProcessed);
% endTime
endTimes = cellfun(@(x,y) x.endTime(y,:),expFile,num2cell(wellsToCompare,2),'UniformOutput',false);
endTimesProcessed = cellfun(@(x) reshape(cell2mat(x),1,numOfWells,numOfChannels),rates,'UniformOutput',false);
[endTimes_mean endTimes_std] = getStdMean(endTimesProcessed);
% MgCurve
MGSignal = cellfun(@(x,y) x.MgCurve(:,y),expFile,num2cell(wellsToCompare,2),'UniformOutput',false);
[MGSignal_mean MGSignal_std] = getStdMean(MGSignal);

%% build ouput struct

mergedExpFile.expFiles     = expFile;
mergedExpFile.Data_mean    = rawData_mean;
mergedExpFile.Data_std     = rawData_std;
mergedExpFile.noBg_mean    = noBackgroundData_mean;
mergedExpFile.noBg_std     = noBackgroundData_std;
mergedExpFile.rate_mean    = rates_mean;
mergedExpFile.rate_std     = rates_std;
mergedExpFile.endTime_mean = endTimes_mean;
mergedExpFile.endTime_std  = endTimes_std;
mergedExpFile.MgCurve_mean = MGSignal_mean;
mergedExpFile.MgCurve_std  = MGSignal_std;

end
