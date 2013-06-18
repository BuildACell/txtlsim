function [data_mean,data_std] = getStdMean(varargin)

if nargin == 1
    data      = varargin{1};
    data_mean = mean(data,2);
    data_std  = std(data')';
    
else
    dataDim = cell2mat(cellfun(@size,varargin','UniformOutput',false));
    
    assert(all(dataDim(:,2) == dataDim(1,2)),'the number of valves should be the same in each experiment')
    
    NumDataPoints = min(dataDim(:,1));
    
    for k= 1:dataDim(1,2)
        data           = cell2mat(cellfun(@(x) x(1:NumDataPoints,1),varargin,'UniformOutput',false) );
        data_mean(:,k) = mean(data,2);
        data_std(:,k)  =  std(data')';
    end
    
end
end