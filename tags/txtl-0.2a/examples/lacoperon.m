% lac operon example file
% model is taken from BFS

%% clean up
clear variables
clc
close all

%% set up the tubes and Species

% Set up the standard TXTL tubes
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('b1');


% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_deGFP = txtl_add_dna(tube3, ...
  'placI(50)', 'rbs(20)', 'betaGal(1000)', 15, 'linear');
dna_gamS = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'gamS(1000)', 10, 'plasmid');


% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3], [20, 8, 2]);


% The concentration of lacI is constant
txtl_addspecies(well_a1, 'protein lacItetramer', 8);
% adding external Lactose
txtl_addspecies(well_a1, 'Lactose_ext', 5);


% Vl = 0.006;
% Kl = 1.25;
% % transporting Lactose_ext into the cell
% Robj1= addreaction(well_a1, 'Lactose_ext -> Lactose');
% Kobj1 = addkineticlaw(Robj1, 'Henri-Michaelis-Menten');
% Pobj1f = addparameter(Kobj1, 'Vl_Lactose_ext',Vl);
% Pobj1r = addparameter(Kobj1, 'Kl_Lactose_ext',Kl);
% set(Kobj1, 'ParameterVariableNames', {'Vl_Lactose_ext', 'Kl_Lactose_ext'});
% set(Kobj1,'SpeciesVariableNames', {'Lactose_ext'});
% 
% get (Robj1, 'ReactionRate')

%% Run the simulation
configsetObj = getconfigset(well_a1, 'active');
simulationTime = 6*60*60;
set(configsetObj, 'StopTime', simulationTime);
set(configsetObj, 'SolverType', 'ode23s');

% 1st run
[t_ode, x_ode, mObj, simData] = txtl_runsim(well_a1, configsetObj,[], []);



%% plot the results


% DNA and mRNA plot
% This should go first to have auto name extraction
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'r-','b-','r--','b--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'Lactose','alloLactose','protein lacItetramer','Glu+Gal'};
dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};



% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,well_a1,dataGroups);
