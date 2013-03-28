% geneexpr_delay.m - basic gene expression reaction
% Vipul Singhal 4 March 2013
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%
close all
clear all

delayedTimeList = [1 2 3 4 5 6]; %[1 2 3 4 5 ];
% t_ode = cell(length(delayedTimeList)+1,1);
% x_ode = cell(length(delayedTimeList)+1,1);
t_ode = cell(length(delayedTimeList)+1,1);
x_ode = cell(length(delayedTimeList)+1,1);
Mobj = cell(length(delayedTimeList)+1,1);
configsetObj = cell(length(delayedTimeList)+1,1);
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');
tube4 = txtl_newtube('gene_expression_delay');
% % Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
 1*10/2.25, ...					% concentration (nM)
  'plasmid');					% type
count = 1;
Mobj{count} = txtl_combine([tube1, tube2, tube3]);%
configsetObj{count} = getconfigset(Mobj{count}, 'active');
simulationTime = 14*60*60;
set(configsetObj{count}, 'SolverType', 'ode23s');
set(configsetObj{count}, 'StopTime', simulationTime);
[t_ode{1},x_ode{1}] = txtl_runsim(Mobj{count},configsetObj{count});

colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--', 'g--'};

%figure(1)
hold on;
idx = findspecies(Mobj{count}, 'protein deGFP*');
%p1 = plot(t_ode{count}/3600, x_ode{count}(:, idx), colors{count});
labels{count} = '0 Hrs delay';
%title('protein deGFP*')
%plots = p1;
count = 2;

dna_deGFP = txtl_add_dna(tube4, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
 0*10/2.25, ...					% concentration (nM)
  'plasmid');					% type

for startTime = delayedTimeList
    Mobj{count} = txtl_combine([tube1,tube2 tube4]);
    configsetObj{count} = getconfigset(Mobj{count}, 'active');
    simulationTime = startTime*60*60;
    set(configsetObj{count}, 'SolverType', 'ode23s');
    set(configsetObj{count}, 'StopTime', simulationTime);
    [t_ode{count},x_ode{count}] = txtl_runsim(Mobj{count},configsetObj{count});
    %Mobj = txtl_combine([Mobj tube3],[7.75 2.25]);
    % this doesnt work since run sim takes the last value in x_ode
    %! TODO vs do something about the addspecies and runsim thing. 
    txtl_addspecies(Mobj{count}, 'DNA p70--rbs--deGFP', 1);
    idx = findspecies(Mobj{count}, 'DNA p70--rbs--deGFP');
    x_ode{count}(end,idx) = 1;
    configsetObj{count} = getconfigset(Mobj{count}, 'active');
    simulationTime = (14-startTime)*60*60;
    set(configsetObj{count}, 'SolverType', 'ode23s');
    set(configsetObj{count}, 'StopTime', simulationTime);
    [t_ode{count},x_ode{count}] = txtl_runsim(Mobj{count},configsetObj{count}, t_ode{count},x_ode{count});%
    count = count +1;
end

% for count = 2:length(delayedTimeList)+1
%   % Plot the time trace  
%   %subplot(2,2,1); 
%   hold on;
%   idx = findspecies(Mobj, 'protein deGFP*');
%   eval(['p' int2str(count) '= plot(t_ode{count}/3600, x_ode{count}(:, idx), colors{count});'])
%   eval(['plots = horzcat(plots, p' int2str(count) ');']);
%   labels{count} = [int2str(delayedTimeList(count-1)) ' Hrs delay'];
%   title('protein deGFP*')
% end
% %title('Time Responses');
% lgh = legend(plots, labels, 'Location', 'Northwest');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [nM]');
% xlabel('Time (Hrs)');
figure(1)
species = 'protein deGFP*';
subplot(2,1,1)
ploteverything
species = 'NTP';
subplot(2,1,2)
ploteverything
% species = 'AA';
% subplot(2,2,3)
% ploteverything
% species = 'Ribo';
% subplot(2,2,4)
% ploteverything  
% figure(2)
%     species = 'DNA p70--rbs--deGFP';
% subplot(2,2,1)
% ploteverything  
%     species = 'RNA rbs--deGFP';
% subplot(2,2,2)
% ploteverything 
%     species = 'Ribo:RNA rbs--deGFP';
% subplot(2,2,3)
% ploteverything 


%    
%     count = 1
% %% plot the result
% % % 
% % DNA and mRNA plot
% dataGroups{1,1} = 'DNA and mRNA';
% dataGroups{1,2} = {'ALL_DNA'}; 
% dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};
% 
% % Gene Expression Plot
% dataGroups{2,1} = 'Gene Expression';
% dataGroups{2,2} = {'protein deGFP*','[protein deGFP]_tot'};
% dataGroups{2,3} = {'g','g--','r-','g--','b-.'};
% 
% % Resource Plot
% dataGroups{3,1} = 'Resource usage';
% 
% txtl_plot(t_ode{count},x_ode{count},Mobj{count},dataGroups);
% figure(5)
%     count = 2
% %% plot the result
% % % 
% % DNA and mRNA plot
% dataGroups{1,1} = 'DNA and mRNA';
% dataGroups{1,2} = {'ALL_DNA'}; 
% dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};
% 
% % Gene Expression Plot
% dataGroups{2,1} = 'Gene Expression';
% dataGroups{2,2} = {'protein deGFP*','[protein deGFP]_tot'};
% dataGroups{2,3} = {'g','g--','r-','g--','b-.'};
% 
% % Resource Plot
% dataGroups{3,1} = 'Resource usage';
% 
% txtl_plot(t_ode{count},x_ode{count},Mobj{count},dataGroups);
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
