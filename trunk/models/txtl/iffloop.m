%% clean up

clear variables
clc
close all
%% My test run

% Set up the standard TXTL tubes
tube1 = txtl_extract('e1');
tube2 = txtl_buffer('b1');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_lacI = txtl_adddna(tube3, ...
    'p70(50)', 'rbs(20)', 'lacI(600)', 3, 'linear');
dna_deGFP = txtl_adddna(tube3, ...
    'placI(50)', 'rbs(20)', 'deGFP(1000)', 3, 'linear');
dna_gamS = txtl_adddna(tube3, ...
    'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3], [22, 15, 1]);


%% Run a simulation
configsetObj = getconfigset(well_a1, 'active');
simulationTime = 4.5*60*60;
set(configsetObj, 'SolverType', 'ode23s');
set(configsetObj, 'StopTime', simulationTime);
% set(configsetObj, 'SolverType', 'ode23s');

% 1st run
[t_ode,x_ode] = txtl_runsim('basic',well_a1,configsetObj,[],[]);

disp('1')
pause(1)
% 2nd run
%txtl_continue_simulation(simData,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70--rbs--lacI', 1);
txtl_addspecies(well_a1, 'DNA placI--rbs--deGFP', 1);
set(configsetObj, 'StopTime', simulationTime/2);
[t_ode_2,x_ode_2] = txtl_runsim('basic',well_a1,configsetObj,t_ode, x_ode);
disp('2')
pause(1)
% t_ode_2 = simData_2.Time;
% x_ode_2 = simData_2.Data;


% 3rd run
%txtl_continue_simulation(simData_2,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70--rbs--lacI', 1);
txtl_addspecies(well_a1, 'DNA placI--rbs--deGFP', 1);
%simData_3 = sbiosimulate(well_a1, configsetObj);
%t_ode_3 = simData_3.Time;
%x_ode_3 = simData_3.Data;
set(configsetObj, 'StopTime', simulationTime*2);
[t_ode_3,x_ode_3] = txtl_runsim('basic',well_a1,configsetObj,t_ode_2,x_ode_2);
disp('3')
pause(1)
% concatante data
t_ode = t_ode_3;%[t_ode; t_ode_2+simulationTime; t_ode_3+simulationTime+simulationTime]; % don't forget to adjust the time
x_ode = x_ode_3;[x_ode; x_ode_2; x_ode_3];

%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
%dataGroups{1,2} = {'#(^DNA (\w+[-=]*)*)'}
dataGroups{1,2} = {'DNA p70--rbs--lacI'};
dataGroups{1,3} = {'b-','r-','b--','r--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','protein gamS','[protein lacI]_tot','[protein gamS]_tot'};
dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups)
