function [endtime_ind,endtime] = expression_endtime(t_vec,data,varargin)

dataChannels = size(data,2);
% hold on
for k=1:dataChannels
    if sum(data(:,k)) > 0
        indend = find(data(end,k)*0.98 < data(:,k));
        endtime_ind(k) = indend(1);
        endtime(k) = t_vec(endtime_ind(k));
    else 
        endtime(k) = 0;
        endtime_ind(k) =0;
    end
end

end

