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
dna_LacI = txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'LacI(600)', 5, 'linear');
dna_deGFP = txtl_dna(tube3, 'placi(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes and add some inducer
well_a1 = txtl_combine([tube1, tube2, tube3], [22, 18, 2]);





%% Run a simulation
configsetObj = getconfigset(well_a1, 'active');
simulationTime = 3*60*60;
set(configsetObj, 'StopTime', simulationTime);
set(configsetObj, 'SolverType', 'ode23s');

% 1st run
simData = sbiosimulate(well_a1, configsetObj);
t_ode = simData.Time;
x_ode = simData.Data;
names = simData.DataNames;

% 2nd run
txtl_continue_simulation(simData,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70=rbs=LacI', 2);
txtl_addspecies(well_a1, 'DNA placi=rbs=deGFP', 2);

simData_2 = sbiosimulate(well_a1, configsetObj);
t_ode_2 = simData_2.Time;
x_ode_2 = simData_2.Data;


% 3rd run
txtl_continue_simulation(simData_2,well_a1);

% add more DNA
txtl_addspecies(well_a1, 'DNA p70=rbs=LacI', 2);
txtl_addspecies(well_a1, 'DNA placi=rbs=deGFP', 2);
simData_3 = sbiosimulate(well_a1, configsetObj);
t_ode_3 = simData_3.Time;
x_ode_3 = simData_3.Data;

% concatnate data
t_ode = [t_ode; t_ode_2+simulationTime; t_ode_3+simulationTime+simulationTime]; % don't forget to adjust the time
x_ode = [x_ode; x_ode_2; x_ode_3];






%% plot the result

figure(1); clf(); subplot(2,1,1);
%txtl_genexpression_plot(simData,well_a1,{'LacI', 'deGFP', 'deGFP*','gamS'})
%title('Gene Expression');

iLacI = findspecies(well_a1, 'protein LacI');
iGFP = findspecies(well_a1, 'protein deGFP');
iGFPs = findspecies(well_a1, 'protein deGFP*');
iGamS = findspecies(well_a1, 'protein gamS');
plot(t_ode/60, x_ode(:, iLacI), 'b-', ...
  t_ode/60, x_ode(:, iGFP) + x_ode(:, iGFPs), 'g--', ...
  t_ode/60, x_ode(:, iGFPs), 'g-', ...
  t_ode/60, x_ode(:, iGamS), 'r-');
% plot(t_ode/60, x_ode(:, iLacI), 'b-',  t_ode/60, x_ode(:, iGamS), 'r-');



title('Gene Expression');
lgh = legend({'LacI', 'GFPt', 'GFP*','GamS'}, 'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

subplot(223)
iDNA_LacI = findspecies(well_a1, 'DNA p70=rbs=LacI');
iDNA_deGFP = findspecies(well_a1, 'DNA placi=rbs=deGFP');
iRNA_LacI = findspecies(well_a1, 'RNA rbs=LacI');
iRNA_deGFP = findspecies(well_a1, 'RNA rbs=deGFP');
plot(t_ode/60, x_ode(:, iDNA_LacI), 'b-', ...
  t_ode/60, x_ode(:, iDNA_deGFP), 'r-', ...
  t_ode/60, x_ode(:, iRNA_LacI), 'b--', ...
  t_ode/60, x_ode(:, iRNA_deGFP), 'r--');

title('DNA and mRNA');
lgh = legend(...
  names([iDNA_LacI, iDNA_deGFP, iRNA_LacI, iRNA_deGFP]), ...
  'Location', 'Northwest');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');



subplot(224)
iNTP = findspecies(well_a1, 'NTP');
iAA  = findspecies(well_a1, 'AA');
iRNAP  = findspecies(well_a1, 'RNAP70');
iRibo  = findspecies(well_a1, 'Ribo');
mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
plot(...
  t_ode/60, x_ode(:, iAA)/x_ode(1, iAA), 'b-', ...
  t_ode/60, x_ode(:, iNTP)/x_ode(1, iNTP), 'r-', ...
  t_ode/60, x_ode(:, iRNAP)/x_ode(1, iRNAP), 'b--', ...
  t_ode/60, x_ode(:, iRibo)/x_ode(1, iRibo), 'r--');

title('Resource usage');
lgh = legend(...
  {'NTP [mM]', 'AA [mM]', 'RNAP70 [nM]', 'Ribo [nM]'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [normalized]');
xlabel('Time [min]');

figure(3)
iDNA_gamS = findspecies(well_a1, 'DNA p70=rbs=gamS');
iRNA_gamS = findspecies(well_a1, 'RNA rbs=gamS');
plot(t_ode/60, x_ode(:, iDNA_gamS), 'r-', ...
    t_ode/60, x_ode(:, iRNA_gamS), 'r--');

title('DNA and mRNA');
lgh = legend(...
  names([iDNA_gamS, iRNA_gamS]), ...
  'Location', 'Northwest');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');


