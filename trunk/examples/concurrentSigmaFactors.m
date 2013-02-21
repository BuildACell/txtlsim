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
dna_sigma28 = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'sigma28(600)', 4*4.2, 'linear');
dna_deGFP =   txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 4*4.2, 'linear');
dna_gamS =    txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1*4.2, 'plasmid');


tube4 = txtl_newtube('no_s28_circuit');
dna_sigma28 = txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'sigma28(600)', 0.4*4.2, 'linear');
dna_deGFP =   txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 4*4.2, 'linear');
dna_gamS =    txtl_add_dna(tube4, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1*4.2, 'plasmid');

tube5 = txtl_newtube('no_s28_circuit');
dna_sigma28 = txtl_add_dna(tube5, 'p70(50)', 'rbs(20)', 'sigma28(600)', 1*4.2, 'linear');
dna_deGFP =   txtl_add_dna(tube5, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 4*4.2, 'linear');
dna_gamS =    txtl_add_dna(tube5, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1*4.2, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3]);
well_b1 = txtl_combine([tube1, tube2, tube4]);
well_c1 = txtl_combine([tube1, tube2, tube5]);

%% Run a simulation
configsetObj_a1 = getconfigset(well_a1, 'active');
simulationTime = 9*60*60;
set(configsetObj_a1, 'SolverType', 'ode23s');
set(configsetObj_a1, 'StopTime', simulationTime);



% 1st run
[t_ode, x_ode, mObj, simData] = txtl_runsim(well_a1, configsetObj_a1);

configsetObj_b1 = getconfigset(well_b1, 'active');
set(configsetObj_b1, 'StopTime', simulationTime);
set(configsetObj_b1, 'SolverType', 'ode23s');

[t_ode_b1, x_ode_b1, mObj_nos28, simData_nos28] = txtl_runsim(well_b1, configsetObj_b1);

configsetObj_c1 = getconfigset(well_c1, 'active');
set(configsetObj_c1, 'StopTime', simulationTime);
set(configsetObj_c1, 'SolverType', 'ode23s');


[t_ode_c1, x_ode_c1, mObj_nos28, simData_nos28] = txtl_runsim(well_c1, configsetObj_c1);


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
txtl_plot(t_ode_b1,x_ode_b1,well_b1,dataGroups)

%%

figure(5)
hold on
plot(t_ode/60,x_ode(:,findspecies(well_b1,'RNAP')),'b*')
plot(t_ode/60,x_ode(:,findspecies(well_b1,'RNAP28')),'r*')
plot(t_ode/60,x_ode(:,findspecies(well_b1,'RNAP70')),'g*')
xlabel('Time [min]');
ylabel('Concentration [nM]');
title('Sigma factors - Simulation')
axis([0 100 0 110])
hold off
legend('free RNAP','RNAP28','RNAP70')

figure(6)
hold on
plot(t_ode/60,x_ode(:,findspecies(well_a1,'protein deGFP*')),'r*')

plot(t_ode_c1/60,x_ode_c1(:,findspecies(well_c1,'protein deGFP*')),'b*')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein deGFP*')),'g*')
legend('4nM \sigma70 - p28','1nM \sigma70 - p28','0.4nM \sigma70 - p28')
xlabel('Time [min]');
ylabel('Concentration [nM]');
title('deGFP expression with presence of competing \sigma28 factor - Simulation')
hold off



