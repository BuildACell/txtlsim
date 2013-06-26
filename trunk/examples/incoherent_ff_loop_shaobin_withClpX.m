% IFFL as tested in lab by S. Guo
% VS 2013

%% clean up

clear all
close all
clc


%% no clpX
% Set up the standard TXTL tubes
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');


tube3 = txtl_newtube('circuit_closed_loop_withClpX');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.05*4.5, 'plasmid');
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP-lva(1000)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.2*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 1000);
txtl_addspecies(tube3, 'aTc', 100);
% set up well_a1
well_a1 = txtl_combine([tube1, tube2, tube3]);

%% Run a simulation
 simulationTime = 12*60*60;

% 1st run
simData = txtl_runsim(well_a1,simulationTime);

 t_ode = simData.Time;
 x_ode = simData.Data;
%% plot the result

%
% txtl_plot(simData, well_a1);
%  iaTc = findspecies(well_a1, 'aTc');
%  iarabinose = findspecies(well_a1, 'arabinose');
%  iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
%  itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
%   figure
%  plot(t_ode, x_ode(:,iaTc), 'r', t_ode, x_ode(:,iarabinose), 'g',...
%      t_ode, x_ode(:,iaraAraC), 'b', t_ode, x_ode(:,itetRaTc), 'k')
%  legend('aTc','arabinose','arabinose:AraC','aTc:tetR')
 ideGFP = findspecies(well_a1, 'protein deGFP-lva*')
 figure
 plot(t_ode, x_ode(:,ideGFP))
 %%
tube3 = txtl_newtube('circuit_closed_loop_withClpX');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.2*4.5, 'plasmid');
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP-lva(1000)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.2*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 1000);
txtl_addspecies(tube3, 'aTc', 100);
% set up well_a1
well_a2 = txtl_combine([tube1, tube2, tube3]);






 
%% Run a simulation

 simulationTime = 12*60*60;


% 1st run
 simData2 = txtl_runsim(well_a2,simulationTime);
 t_ode2 = simData2.Time;
 x_ode2 = simData2.Data;
%% plot the result
%
% txtl_plot(simData2,well_a2);
  ideGFP2 = findspecies(well_a2, 'protein deGFP-lva*')
 figure
 plot(t_ode2, x_ode2(:,ideGFP2))
%  iaTc = findspecies(well_a1, 'aTc');
%  iarabinose = findspecies(well_a1, 'arabinose');
%  iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
%  itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
%   figure
%  plot(t_ode, x_ode(:,iaTc), 'r', t_ode, x_ode(:,iarabinose), 'g',...
%      t_ode, x_ode(:,iaraAraC), 'b', t_ode, x_ode(:,itetRaTc), 'k')
%  legend('aTc','arabinose','arabinose:AraC','aTc:tetR')
 %%
 
 tube3 = txtl_newtube('circuit_closed_loop_withClpX');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.8*4.5, 'plasmid');
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP-lva(1000)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.2*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 1000);
txtl_addspecies(tube3, 'aTc', 100);
% set up well_a1
well_a3 = txtl_combine([tube1, tube2, tube3]);


 simulationTime = 12*60*60;
 simData3 = txtl_runsim(well_a3,simulationTime);
 t_ode3 = simData3.Time;
 x_ode3 = simData3.Data;
%% plot the result
%
% txtl_plot(simData3,well_a3);

  ideGFP3 = findspecies(well_a3, 'protein deGFP-lva*')
 figure
 plot(t_ode3, x_ode3(:,ideGFP3))
%  iaTc = findspecies(well_a1, 'aTc');
%  iarabinose = findspecies(well_a1, 'arabinose');
%  iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
%  itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
%   figure
%  plot(t_ode, x_ode(:,iaTc), 'r', t_ode, x_ode(:,iarabinose), 'g',...
%      t_ode, x_ode(:,iaraAraC), 'b', t_ode, x_ode(:,itetRaTc), 'k')
%  legend('aTc','arabinose','arabinose:AraC','aTc:tetR')
 %%
 tube3 = txtl_newtube('circuit_closed_loop_withClpX');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', 1.6*4.5, 'plasmid');
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP-lva(1000)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.2*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 1000);
txtl_addspecies(tube3, 'aTc', 100);
% set up well_a1
well_a4 = txtl_combine([tube1, tube2, tube3]);

 
 simulationTime = 12*60*60;
 simData4 = txtl_runsim(well_a4,simulationTime);
 t_ode4 = simData4.Time;
 x_ode4 = simData4.Data;
%% plot the result

% txtl_plot(simData4,well_a4);
  ideGFP4 = findspecies(well_a4, 'protein deGFP-lva*');
 figure
 plot(t_ode4, x_ode4(:,ideGFP4))
 
 figure
plot(t_ode/60, x_ode(:,33), 'c', t_ode2/60, x_ode2(:,33), 'g', t_ode3/60, x_ode3(:,33), 'b', t_ode4/60, x_ode4(:,33), 'r');
 legend('DNA tetR = 0.225 nM', 'DNA tetR = 0.9 nM', 'DNA tetR = 3.6 nM', 'DNA tetR = 7.2 nM')
 title('IFFL pulse in the presence of protein degradation')
 ylabel('deGFP signal (arbitrary units)');
 xlabel('time (min)');
%  iaTc = findspecies(well_a1, 'aTc');
%  iarabinose = findspecies(well_a1, 'arabinose');
%  iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
%  itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
%   figure
%  plot(t_ode, x_ode(:,iaTc), 'r', t_ode, x_ode(:,iarabinose), 'g',...
%      t_ode, x_ode(:,iaraAraC), 'b', t_ode, x_ode(:,itetRaTc), 'k')
%  legend('aTc','arabinose','arabinose:AraC','aTc:tetR')

% 
% --------------------
% 
% % IFFL as tested in lab by S. Guo
% % VS 2013
% 
% %% clean up
% 
% 
% numSims = 4
% T = cell(numSims,1);
% X = cell(numSims,1);
% tetR_ICs = [0.002 0.01 0.05 0.2 0.8 2];
% %% no clpX
% % Set up the standard TXTL tubes
% tube1 = txtl_extract('E9');
% tube2 = txtl_buffer('E9');
% 
% 
% tube3 = txtl_newtube('circuit_closed_loop_withClpX');
% 
% for i = 1:numSims
% txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', tetR_ICs(i)*4.5, 'plasmid');
% txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
% txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP-lva(1000)',0.5*4.5, 'plasmid');
% txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.2*4.5, 'plasmid');
% txtl_addspecies(tube3, 'arabinose', 1000);
% txtl_addspecies(tube3, 'aTc', 100);
% % set up well_a1
% well_a1 = txtl_combine([tube1, tube2, tube3]);
% 
% %% Run a simulation
% configsetObj = getconfigset(well_a1, 'active');
%  simulationTime = 12*60*60;
% set(configsetObj, 'SolverType', 'ode23s');
% set(configsetObj, 'StopTime', simulationTime);
% 
% % 1st run
%  [t_ode,x_ode] = txtl_runsim(well_a1,configsetObj);
% 
%  ideGFP = findspecies(well_a1, 'protein deGFP-lva*');
% T{i} = t_ode;
% X{i} = x_ode(:,ideGFP);
% 
% end
% 
% plot(T{1}, X{1}, 'c', T{2}, X{2}, 'g', T{3}, X{3}, 'b', T{4}, X{4}, 'r');
% 
% 





