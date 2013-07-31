function [rates,fitInfo] = expression_rate(t_vec,data,varargin)

dataChannels = size(data,2);
if nargin > 2 && strcmp(varargin{1},'DEBUG')
    displayMode = 2;
else
    displayMode = 1;
end

if displayMode == 2
 hold on 
end
for k=1:dataChannels
    if sum(data(:,k)) > 0
        ind10 = find(data(end,k)*0.10 < data(:,k));
        ind90 = find(data(end,k)*0.65 < data(:,k));
         
    
        f = fittype('poly1');
        [fitInfo(k).ff fitInfo(k).gof] = fit(t_vec(ind10(1):ind90(1)),data(ind10(1):ind90(1),k),f);
        if displayMode == 2
         plot(t_vec(ind10(1):ind90(1)),fitInfo(k).ff.p1*t_vec(ind10(1):ind90(1))+fitInfo(k).ff.p2,'r','LineWidth',2)
         plot(t_vec,data(:,k))
        end
        rates(k) = fitInfo(k).ff.p1;
    end
end

end

