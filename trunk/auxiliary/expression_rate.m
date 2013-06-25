function [rates,fitInfo] = expression_rate(t_vec,data,varargin)

dataChannels = size(data,2);
% hold on 
for k=1:dataChannels
    if sum(data(:,k)) > 0
        ind10 = find(data(end,k)*0.10 < data(:,k));
        ind90 = find(data(end,k)*0.65 < data(:,k));
         
    
        f = fittype('poly1');
        [fitInfo(k).ff fitInfo(k).gof] = fit(t_vec(ind10(1):ind90(1)),data(ind10(1):ind90(1),k),f);
%         plot(t_vec(ind10(1):ind90(1)),fitInfo(k).ff.p1*t_vec(ind10(1):ind90(1))+fitInfo(k).ff.p2,'LineWidth',2)
%         plot(t_vec,data(:,k))
        rates(k) = fitInfo(k).ff.p1;
    end
end

end

