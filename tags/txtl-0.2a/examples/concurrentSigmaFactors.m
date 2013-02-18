%% clean up

clear variables
close all
clc

%% simple circuit with 2 genes

% Set up the standard TXTL tubes
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('E6');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_sigma28 = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'sigma28(600)', 4, 'linear');
dna_deGFP =   txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 4, 'linear');
dna_gamS =    txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

tube4 = txtl_newtube('no_s28_circuit');
txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 4, 'linear');
txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3]);
well_b1 = txtl_combine([tube1, tube2, tube4]);

%% Run a simulation
configsetObj_a1 = getconfigset(well_a1, 'active');
simulationTime = 5.5*60*60;
set(configsetObj_a1, 'StopTime', simulationTime);



% 1st run
[t_ode, x_ode, mObj, simData] = txtl_runsim(well_a1, configsetObj_a1);

configsetObj_b1 = getconfigset(well_b1, 'active');
set(configsetObj_b1, 'StopTime', simulationTime);

[t_ode_nos28, x_ode_nos28, mObj_nos28, simData_nos28] = txtl_runsim(well_b1, configsetObj_b1);


%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','m-','k-'};



% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*'};
dataGroups{2,3} = {'b-','g-','r-','g--','b--','b-.','m-','k+'};



% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups)
txtl_plot(t_ode_nos28,x_ode_nos28,well_b1,dataGroups)






