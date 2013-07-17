% IFFL as tested in lab by S. Guo
% VS 2013

%% clean up

clear variables
clc
close all

%% no clpX
% Set up the standard TXTL tubes
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

%% tetR and deGFP
% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_open_loop_noClpX');
txtl_add_dna(tube3, ...
    'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.25*4.5, 'plasmid');

% Mix the contents of the individual tubes
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
 iGFP = findspecies(well_a1, 'protein deGFP*')
T1 = t_ode;
GFP1 = x_ode(:,iGFP);
 
 %% AraC Arabinose tetR deGFP
 % Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_open_loop_noClpX');
txtl_add_dna(tube3, ...
    'p70(50)', 'rbs(20)', 'AraC(600)', 0.1*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.25*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 3000);

% Mix the contents of the individual tubes
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
 figure
 iarabinose = findspecies(well_a1, 'arabinose');
 iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
 plot(t_ode, x_ode(:,iarabinose), 'g',...
     t_ode, x_ode(:,iaraAraC), 'b')
 legend('arabinose','arabinose:AraC')
 iGFP = findspecies(well_a1, 'protein deGFP*')
T2 = t_ode;
GFP2 = x_ode(:,iGFP);
% %% deGFP
%  % Set up a tube that will contain our DNA
% tube3 = txtl_newtube('circuit_open_loop_noClpX');
% txtl_add_dna(tube3, ...
%     'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
% 
% % Mix the contents of the individual tubes
% well_a1 = txtl_combine([tube1, tube2, tube3]);
% 
% % Run a simulation
% configsetObj = getconfigset(well_a1, 'active');
%  simulationTime = 12*60*60;
% set(configsetObj, 'SolverType', 'ode23s');
% set(configsetObj, 'StopTime', simulationTime);
% 
% % 1st run
%  [t_ode,x_ode] = txtl_runsim(well_a1,configsetObj);
%  
% % plot the result
% % DNA and mRNA plot
% dataGroups{1,1} = 'DNA and mRNA';
% dataGroups{1,2} = {'ALL_DNA'};
% dataGroups{1,3} = {'b-','r-','g--','r--','y-','c-','g-','g--'};
% 
% % Gene Expression Plot
% dataGroups{2,1} = 'Gene Expression';
% dataGroups{2,2} = {'protein deGFP*','protein tetRdimer'};
% %dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};
% 
% % Resource Plot
% dataGroups{3,1} = 'Resource usage';
%  txtl_plot(t_ode,x_ode,well_a1,dataGroups);
%  

 %% tetR deGFP aTc
% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_open_loop_noClpX');
txtl_add_dna(tube3, ...
    'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.25*4.5, 'plasmid');
txtl_addspecies(tube3, 'aTc', 1500);

% Mix the contents of the individual tubes
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
 
  iaTc = findspecies(well_a1, 'aTc');
 itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
  figure
  plot(t_ode, x_ode(:,iaTc), 'r',  t_ode, x_ode(:,itetRaTc), 'k')
 legend('aTc','aTc:tetR')
   iGFP = findspecies(well_a1, 'protein deGFP*')
T3 = t_ode;
GFP3 = x_ode(:,iGFP);
 %% AraC arabinose deGFP
 % Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_open_loop_noClpX');
txtl_add_dna(tube3, ...
    'p70(50)', 'rbs(20)', 'AraC(600)', 0.1*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
txtl_addspecies(tube3, 'arabinose', 3000);

% Mix the contents of the individual tubes
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
 iarabinose = findspecies(well_a1, 'arabinose');
 iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
  figure
  plot(t_ode, x_ode(:,iarabinose), 'g',t_ode, x_ode(:,iaraAraC));
 legend('arabinose','arabinose:AraC')
  iGFP = findspecies(well_a1, 'protein deGFP*')
T4 = t_ode;
GFP4 = x_ode(:,iGFP);
%  %% AraC, arabinose, tetR, aTc, deGFP
%  % Set up a tube that will contain our DNA
% tube3 = txtl_newtube('circuit_open_loop_noClpX');
% txtl_add_dna(tube3, ...
%     'p70(50)', 'rbs(20)', 'AraC(600)', 0.1*4.5, 'plasmid');
% txtl_add_dna(tube3, ...
%     'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)', 0.5*4.5, 'plasmid');
% txtl_add_dna(tube3, ...
%     'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.25*4.5, 'plasmid');
% txtl_addspecies(tube3, 'arabinose', 3000);
% txtl_addspecies(tube3, 'aTc', 1500);
% 
% % Mix the contents of the individual tubes
% well_a1 = txtl_combine([tube1, tube2, tube3]);
% 
% % Run a simulation
% configsetObj = getconfigset(well_a1, 'active');
%  simulationTime = 12*60*60;
% set(configsetObj, 'SolverType', 'ode23s');
% set(configsetObj, 'StopTime', simulationTime);
% 
% % 1st run
%  [t_ode,x_ode] = txtl_runsim(well_a1,configsetObj);
% % plot the result
% % DNA and mRNA plot
% dataGroups{1,1} = 'DNA and mRNA';
% dataGroups{1,2} = {'ALL_DNA'};
% dataGroups{1,3} = {'b-','r-','g--','r--','y-','c-','g-','g--'};
% 
% % Gene Expression Plot
% dataGroups{2,1} = 'Gene Expression';
% dataGroups{2,2} = {'protein deGFP*','protein tetRdimer'};
% %dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};
% 
% % Resource Plot
% dataGroups{3,1} = 'Resource usage';
%  txtl_plot(t_ode,x_ode,well_a1,dataGroups);
%  
%  iaTc = findspecies(well_a1, 'aTc');
%  iarabinose = findspecies(well_a1, 'arabinose');
%  iaraAraC = findspecies(well_a1, 'arabinose:protein AraC');
%  itetRaTc = findspecies(well_a1, '2 aTc:protein tetRdimer');
%   figure
%  plot(t_ode, x_ode(:,iaTc), 'r', t_ode, x_ode(:,iarabinose), 'g',...
%      t_ode, x_ode(:,iaraAraC), 'b', t_ode, x_ode(:,itetRaTc), 'k')
%  legend('aTc','arabinose','arabinose:AraC','aTc:tetR')
%  
figure
 plot(T1/60, GFP1, 'm', T2/60, GFP2, 'y', T3/60, GFP3, 'b', T4/60, GFP4, 'g')
 legend('tetR, deGFP',...
     'AraC, arabinose, tetR, deGFP',...
     'tetR, deGFP, aTc',...
     'AraC, arabinose, deGFP');
 ylabel('deGFP signal (arbitrary units)');
 xlabel('time (min)')
 title('deGFP signal as a function of different simulated experimental conditions')
 
 
 %}
 
 
 
 
 
 
 
 
 %%
 


