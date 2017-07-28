% Example created by Zoltan A. Tuza Aug 2014
%
% This example imports Biotek Plate Reader data (exported to excel, then converted to csv file)
% delimiter of the CSV file is also adjusted to ';' in order to prevent data loss

%% import
% the csv file is saved into the "data" folder, reading may take couple of
% seconds
expData = wholeExpfileReader('data/e15_pr_gfp_pr_gfp_mg_grad_082113_data.csv',';');

% expData struct cotains basic information about the experiment such as
% * length of the experiment
% * reading intervals
% * number of reads
% * recorded chanels with name
% * rate of change
% * endTime 
expData


%% plot data
% the Data matrix is organized as follows: 
% number of reads x number of wells x number channels
% For example, for an experiment that runs for 840 min sampled at every 180
% seconds, that gives us ~281 reads. Let's assume we've run 14 different
% well and measuring 3 diffetent emission wavelength. Given that the size of our Data
% matrix is 281x14x3
figure('Name','GFP channel')
plot(expData.t_vec/60,expData.Data(:,3:end,2))
title(sprintf('Data source: %s',expData.FileName),'interpreter','none')
xlabel('Time [min]');
ylabel('GFP A.U.');

%% plot STD with data
% Let's assume that we've run the same experiments multiple time, either
% running the same sample in different well or in actually different
% experiments.
% We can merge these reapetes into one single data structure and get some
% basic statistics with it.

colorCodes = {'r','b','g','c','m','k',[1 1 .5],[.7 .5 .2],[0 1 .2],[.35 .8 .8],[.9 0 .4],[1 .2 .2]};
intensity = 0.3;
% h = stdshade(expData.t_vec/60,expData.Data(:,3:end,2),std,intensity,colorCodes);


%% basic data analysis

% Plot the time when protein production stops (defined as 95% of the end point signal value)
figure('Name','expression end time')
bar(expData.endTime(:,2)/60)
title(sprintf('Data source: %s',expData.FileName),'interpreter','none')
set(gca,'XTickLabel',expData.wellNames)
ylabel('End time [min]')
xlabel('wells on the plate')

% Plot GFP expression rate (defined as the slop of a linear curve fitted between 10% and 65% of the total signal value)
figure('Name','expression rate')
bar(expData.rate(:,2))
title(sprintf('Data source: %s',expData.FileName),'interpreter','none')
set(gca,'XTickLabel',expData.wellNames)
ylabel('Expression rate [A.U./sec]')
xlabel('wells on the plate')



