function [data_mean,data_std] = getStdMean(data)

data_mean = mean(data,2);
data_std =  std(data')';
end