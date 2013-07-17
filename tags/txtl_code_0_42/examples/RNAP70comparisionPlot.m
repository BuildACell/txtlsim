function RNAP70comparisionPlot(Mobj,t_ode, x_ode, parameterName, parameterValue, colororder)
% RNAP70species = {'RNAP70',...
%  'RNAP70:DNA pJ23119--att-rbs--deGFP',...
%   'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP',...
%   'RNAP70:DNA pJ23119--control',...
%   'NTP:RNAP70:DNA pJ23119--control',...
%   'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att',...
%   'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att'};
% markers = {'--', ':', '.-', '-','+', '--*', '--.','--o'};%,  'x', 's', 'd'
figure('OuterPosition',[100 100 880 620])
numParam = length(Mobj);
p = zeros(4,numParam);
labels = cell(numParam*4,1);
% maxval = 0;
% transientTime = 0.03;% minutes
for i = 1:numParam
    hold on
    idx1 = findspecies(Mobj{i}, 'RNAP70');
    idx2 = findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--att-rbs--deGFP'); 
    idx3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP'); 
    idx4 = findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--control'); 
    idx5 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA pJ23119--control'); 
    idx6 = findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att'); 
    idx7 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att');
    
   [p(1,i)] = plot(t_ode{i}/60, x_ode{i}(:,idx1),'-','Color', colororder(i,:),'LineWidth',2);
   [p(2,i)] = plot(t_ode{i}/60, x_ode{i}(:,idx2)+x_ode{i}(:,idx3), '--','Color', colororder(i,:),'LineWidth',2);
   [p(3,i)] = plot(t_ode{i}/60, x_ode{i}(:,idx4)+x_ode{i}(:,idx5), '-.','Color', colororder(i,:),'LineWidth',2);
   [p(4,i)] = plot(t_ode{i}/60, x_ode{i}(:,idx6)+x_ode{i}(:,idx7), ':','Color', colororder(i,:),'LineWidth',2);
   labels{4*(i-1)+1} = [parameterName ' = ' num2str(parameterValue(i))  ' RNAP70'];
   labels{4*(i-1)+2} = [parameterName ' = ' num2str(parameterValue(i))  ' RNAP70\_GFP'];
   labels{4*(i-1)+3} = [parameterName ' = ' num2str(parameterValue(i))  ' RNAP70\_control'];
   labels{4*(i-1)+4} = [parameterName ' = ' num2str(parameterValue(i))  ' RNAP70\_GFP\_RNAatt'];

%     ,..., 
%         t_ode{i}/60, x_ode{i}(:,idx4)+x_ode{i}(:,idx5),...
%         t_ode{i}/60, x_ode{i}(:,idx6)+x_ode{i}(:,idx7),'Color', colororder(i,:),'LineWidth',2,...
%         t_ode{i}/60,  x_ode{i}(:,idx1)+x_ode{i}(:,idx2)+x_ode{i}(:,idx3)+x_ode{i}(:,idx4)+x_ode{i}(:,idx5)+x_ode{i}(:,idx6)+x_ode{i}(:,idx7),...
%         'Color', colororder(i,:),'LineWidth',2);%,'MarkerSymbol', '*'
%         
%     labels{5*(i-1)+1} = [parameterName ' = ' num2str(parameterValues(i)) '_RNAP70'];
%     labels{5*(i-1)+2} = [parameterName ' = ' num2str(parameterValues(i)) '_RNAPbound_deGFP'];
%     labels{5*(i-1)+3} = [parameterName ' = ' num2str(parameterValues(i)) '_RNAPbound_control'];
%     labels{5*(i-1)+4} = [parameterName ' = ' num2str(parameterValues(i)) '_RNAPbound_att'];
%     labels{5*(i-1)+5} = [parameterName ' = ' num2str(parameterValues(i)) '_totalRNAP70'];
%     
%     
%     timeIdx = find(t_ode{i}/60>transientTime,1);
%     m = max(x_ode{i}(1:timeIdx-1,idx));
%     if maxval<m
%         maxval = m;
%     end
end
% title(figTitle);
% axis([0 transientTime 0 (eps+maxval)*1.1])
% lgh = legend(labels{:}, 'Location', 'Northwest'); 
% legend(lgh, 'boxoff');
% xlabel('time (min)')
% print('-dtiff','-r200',[folderdate fileIdentifier])
% legend(h,M) associates each row of the matrix or cell array of strings M 
% with the corresponding graphics object (patch or line) in the vector of handles h.

h = reshape(p, 4*numParam,1);
legend(h, labels, 'Location', 'NorthEastOutside')
title('RNAP-DNA binding dynamics for various control DNA loading conditions')
xlabel('time/minutes')
ylabel('species conc')



end
