%% clean up

% clear variables
% clc
% close all
%% My test run

% Set up the standard TXTL tubes
tube1 = txtl_extract('E7');
tube2 = txtl_buffer('E7');

% Set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit_open_loop');
txtl_add_dna(tube3, ...
    'p70(50)', 'rbs(20)', 'sigma28(600)', 0.2*4.5, 'plasmid');
txtl_add_dna(tube3, ...
    'p28_ptet(150)', 'rbs(20)', 'deGFP(1000)', 1*4.5, 'plasmid');



tube4 = txtl_newtube('circuit_closed_loop');
txtl_add_dna(tube4, ...
    'p28(50)', 'rbs(20)', 'tetR(600)', 0.01*4.5, 'plasmid');
txtl_add_dna(tube4, ...
    'p70(50)', 'rbs(20)', 'sigma28(600)',0.2*4.5, 'plasmid');
txtl_add_dna(tube4,'p28_ptet(150)', 'rbs(20)', 'deGFP(1000)',1*4.5, 'plasmid');




% Mix the contents of the individual tubes
 well_a1 = txtl_combine([tube1, tube2, tube3]);

% set up well_b1
well_b1 = txtl_combine([tube1, tube2, tube4]);






 
%% Run a simulation
configsetObj = getconfigset(well_a1, 'active');
 simulationTime = 12*60*60;
set(configsetObj, 'SolverType', 'ode23s');
set(configsetObj, 'StopTime', simulationTime);

% 1st run
 [t_ode,x_ode] = txtl_runsim(well_a1,configsetObj);

configsetObj_b1 = getconfigset(well_b1, 'active');

set(configsetObj_b1, 'SolverType', 'ode23s');
set(configsetObj_b1, 'StopTime', simulationTime);


[t_ode_b1,x_ode_b1] = txtl_runsim(well_b1,configsetObj_b1);


%% plot the result
close all
% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'b-','r-','g--','r--','y-','c-','g-','g--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','protein tetRdimer'};
%dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';
%
 txtl_plot(t_ode,x_ode,well_a1,dataGroups);

 txtl_plot(t_ode_b1,x_ode_b1,well_b1,dataGroups);

% figure(3)
% hold on
% plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein tetR')),'b')
% plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein deGFP*')),'g')
% plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein tetRdimer')),'r')
% xlabel('Time [min]');
% ylabel('Concentration [nM]');
% hold off
% legend('tetR','deGFP*','tetRdimer')

figure(4)
hold on
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'RNAP28:DNA p28_ptet--rbs--deGFP')),'b')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'DNA p28_ptet--rbs--deGFP:protein tetRdimer')),'r')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'RNAP28:DNA p28_ptet--rbs--deGFP:protein tetRdimer')),'g')
xlabel('Time [min]');
ylabel('Concentration [nM]');
hold off
legend('DNA:RNAP28','DNA:tetR','DNA:RNAP28:tetR','Location','Best')


figure(5)
hold on
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'RNAP')),'b')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'RNAP28')),'r')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'RNAP70')),'g')
xlabel('Time [min]');
ylabel('Concentration [nM]');
hold off
legend('RNAP','RNAP28','RNAP70')


figure(6)
hold on
plot(t_ode/60,x_ode(:,findspecies(well_a1,'protein deGFP*')),'b')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein deGFP*')),'r')

xlabel('Time [min]');
ylabel('Concentration [nM]');
hold off
legend('open loop - deGFP*','closed loop - deGFP*','Location','Best')


