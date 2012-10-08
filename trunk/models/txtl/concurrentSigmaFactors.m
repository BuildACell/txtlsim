%% clean up

clear variables
close all
clc

%% simple circuit with 2 genes

% Set up the standard TXTL tubes
tube1 = txtl_extract('e1');
tube2 = txtl_buffer('b1');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
%dna_sigma28 = txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'sigma28(600)', 10, 'linear');
dna_deGFP =   txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 10, 'linear');
dna_gamS =    txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3], [22, 8, 3]);

txtl_setup_parameters(well_a1);


%% Run a simulation
configsetObj = getconfigset(well_a1, 'active');
simulationTime = 1.5*60*60;
set(configsetObj, 'StopTime', simulationTime);
set(configsetObj, 'SolverType', 'ode23s');

% 1st run
simData = sbiosimulate(well_a1, configsetObj);
t_ode = simData.Time;
x_ode = simData.Data;
names = simData.DataNames;


%2dn run
% txtl_continue_simulation(simData,well_a1);
% 
% %txtl_addspecies(well_a1, 'RNAP', 10);
% simData_2 = sbiosimulate(well_a1, configsetObj);
% t_ode_2 = simData_2.Time;
% x_ode_2 = simData_2.Data;
% 
% t_ode = [t_ode; t_ode_2+simulationTime]; % don't forget to adjust the time
% x_ode = [x_ode; x_ode_2];


%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'#(^DNA (\w+[-=]*)*)'}
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-'}



% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*'};
dataGroups{2,3} = {'b-','g-','r-','g--','b--','b-.'}



% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups)




