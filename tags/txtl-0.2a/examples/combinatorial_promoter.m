%% clean up

clear variables
clc
close all
%% My test run

% Set up the standard TXTL tubes
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('E6');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_notetR');
txtl_add_dna(tube3, ...
    'p70(50)', 'rbs(20)', 'sigma28(600)', 5, 'linear');
txtl_add_dna(tube3, ...
    'p28_ptet(150)', 'rbs(20)', 'deGFP(1000)-lva-terminator', 5, 'linear');
txtl_add_dna(tube3, ...
'p70(50)', 'rbs(20)', 'gamS(1000)', 2, 'plasmid');


tube4 = txtl_newtube('circuit_with_tetR');
txtl_add_dna(tube4, ...
    'p70(50)', 'rbs(20)', 'tetR(600)', 3, 'linear');
txtl_add_dna(tube4, ...
    'p70(50)', 'rbs(20)', 'sigma28(600)', 5, 'linear');
txtl_add_dna(tube4, ...
    'p28_ptet(150)', 'rbs(20)', 'deGFP(1000)-lva-terminator', 5, 'linear');
txtl_add_dna(tube4, ...
'p70(50)', 'rbs(20)', 'gamS(1000)', 3, 'plasmid');


% Mix the contents of the individual tubes
well_a1 = txtl_combine([tube1, tube2, tube3]);



% set up well_b1
well_b1 = txtl_combine([tube1, tube2, tube4]);




 
%% Run a simulation
configsetObj = getconfigset(well_a1, 'active');
simulationTime = 10*60*60;

set(configsetObj, 'StopTime', simulationTime);

% 1st run
[t_ode,x_ode] = txtl_runsim(well_a1,configsetObj);


configsetObj_b1 = getconfigset(well_b1, 'active');
set(configsetObj_b1, 'StopTime', simulationTime);


[t_ode_b1,x_ode_b1] = txtl_runsim(well_b1,configsetObj_b1);


%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP-lva-terminator*'};%,'protein tetRtetramer'};
%dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups);

txtl_plot(t_ode_b1,x_ode_b1,well_b1,dataGroups);

figure(3)
hold on
plot(t_ode/60,x_ode(:,findspecies(well_a1,'protein deGFP-lva-terminator*')),'b')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein deGFP-lva-terminator*')),'r')
xlabel('Time [min]');
ylabel('Concentration [nM]');
hold off
legend('deGFP_no_tetR','deGFP_tetR')


