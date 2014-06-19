%% clean up

clearvars
close all
%% My test run

% Set up the standard TXTL tubes
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

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
tic
[simData] = txtl_runsim(well_a1,10*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

simData_b1= txtl_runsim(well_b1,10*60*60);
t_ode_b1 = simData_b1.Time;
x_ode_b1 = simData_b1.Data;




%% plot the result


txtl_plot(t_ode,x_ode,well_a1);

txtl_plot(t_ode_b1,x_ode_b1,well_b1);

figure(4)
hold on
plot(t_ode/60,x_ode(:,findspecies(well_a1,'protein deGFP-lva-terminator*')),'b')
plot(t_ode_b1/60,x_ode_b1(:,findspecies(well_b1,'protein deGFP-lva-terminator*')),'r')
xlabel('Time [min]');
ylabel('Concentration [nM]');
hold off
legend('deGFP_no_tetR','deGFP_tetR')


