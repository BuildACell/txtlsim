function plotAcrossSims2(Mobj,t_ode, x_ode, speciesName, fileIdentifier, figTitle, parameterName1, ...
    parameterValues1, parameterName2, parameterValues2, colororder)
figure('OuterPosition',[100 100 880 620]) 
numParam = length(Mobj);
p = cell(numParam,1);
labels = cell(numParam,1);
maxval = 0;
transientTime = 0.03;% minutes
for i = 1:numParam
    hold on 
    idx = findspecies(Mobj{i}, speciesName);
p{i} = plot(t_ode{i}/60, x_ode{i}(:,idx), 'Color', colororder(i,:), 'linewidth', 2);
labels{i} = [parameterName1 ' = ' num2str(parameterValues1(i)) ' ' parameterName2 '=' num2str(parameterValues2(i))];
timeIdx = find(t_ode{i}/60>transientTime,1);
m = max(x_ode{i}(1:timeIdx-1,idx));
if maxval<m
    maxval = m;
end
end
title(figTitle);
axis([0 transientTime 0 (eps+maxval)*1.1])
lgh = legend(labels{:}, 'Location', 'Northwest'); 
legend(lgh, 'boxoff');
xlabel('time (min)')
print('-dtiff','-r200',[fileIdentifier])
close


figure('OuterPosition',[100 100 880 620]) 
numParam = length(Mobj);
p = cell(numParam,1);
labels = cell(numParam,1);

for i = 1:numParam
    hold on 
    idx = findspecies(Mobj{i}, speciesName);
p{i} = plot(t_ode{i}/60, x_ode{i}(:,idx), 'Color', colororder(i,:), 'linewidth', 2);
labels{i} = [parameterName1 ' = ' num2str(parameterValues1(i)) ' ' parameterName2 '=' num2str(parameterValues2(i))];

end
title(figTitle);

lgh = legend(labels{:}, 'Location', 'Northwest'); 
legend(lgh, 'boxoff');
xlabel('time (min)')
print('-dtiff','-r200',[fileIdentifier '_full'])
close
end