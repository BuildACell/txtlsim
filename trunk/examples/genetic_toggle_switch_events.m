% Genetic Toggle Switch example
% Vipul Singhal, September 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch. It shows the bistability of the circuit, and the
% switching that can be accomplished by the tetR repressing inducer aTc and
% the lacI repressing inducer IPTG. 
%
close all 
clear all
clc


tetR_initialConc = 10; 
lacI_initialConc = 1; 
simulationDuration = 5*60*60;
inductionTimePoints = [1*60*60 3*60*60];
inducerConc = 9; % setting this to 10 gives strange fluctuations in the AA conc

% Set up the standard TXTL tubes
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('b1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

dna_lacI = txtl_add_dna(tube3,'thio-junk(500)-ptet(50)', 'rbs(20)', 'lacI(647)-lva(40)-terminator(100)', 5, 'linear');
dna_tetR = txtl_add_dna(tube3, 'thio-junk(500)-ptrc2(50)', 'rbs(20)', 'tetR(647)-lva(40)-terminator(100)', 5, 'linear');
dna_deGFP = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_add_dna(tube3,  'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Set initial lacI and tetR concentrations
lacIprotein = sbioselect(tube3, 'Name','protein lacI-lva-terminator');
tetRprotein = sbioselect(tube3, 'Name','protein tetR-lva-terminator');
set(lacIprotein, 'InitialAmount', lacI_initialConc);
set(tetRprotein, 'InitialAmount', tetR_initialConc);

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);
configsetObj = getconfigset(Mobj, 'active');

set(configsetObj, 'StopTime', simulationDuration)
% set up events triggers and events functions as arrays of cells. which in
% turn contain strings
eventTriggers = [{['time >= ' num2str(inductionTimePoints(1))]}, {['time >= ' num2str(inductionTimePoints(2))]}];
eventFunctions = [{['aTc = ' num2str(inducerConc)]}, {['IPTG= ' num2str(inducerConc)]}];

% run the events script

[t_ode, x_ode, simData] = txtl_runsim_events(Mobj, configsetObj, eventTriggers, eventFunctions);


%% Plot results
% Top row: protein and RNA levels
figure(1); clf(); subplot(2,1,1);
ilacI = findspecies(Mobj, 'protein lacI-lva-terminator');
itetR = findspecies(Mobj, 'protein tetR-lva-terminator');
iGamS = findspecies(Mobj, 'protein gamS');
iGFP = findspecies(Mobj, 'protein deGFP');
iGFPs = findspecies(Mobj, 'protein deGFP*');

p = plot(t_ode/60, x_ode(:, itetR),'k-', t_ode/60, x_ode(:, ilacI), 'b-', t_ode/60, x_ode(:, iGamS), 'r-', ...
  t_ode/60, x_ode(:, iGFP) + x_ode(:, iGFPs), 'g--', ...
  t_ode/60, x_ode(:, iGFPs), 'g-');

title('Gene Expression');
lgh = legend({'tetR', 'lacI', 'GamS', 'GFPt', 'GFP*'}, 'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
axis([0 simulationDuration/60 0 10])

% Second row, left: resource limits
subplot(2,2,3);
iNTP = findspecies(Mobj, 'NTP');
iAA  = findspecies(Mobj, 'AA');
iRNAP  = findspecies(Mobj, 'RNAP70');
iRibo  = findspecies(Mobj, 'Ribo');
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

% Second row, right: DNA and mRNA
subplot(2,2,4);
iDNA_tetR = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator');
iDNA_lacI = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator');
iDNA_gamS = findspecies(Mobj, 'DNA p70--rbs--gamS');
iRNA_tetR = findspecies(Mobj, 'RNA rbs--tetR-lva-terminator');
iRNA_lacI = findspecies(Mobj, 'RNA rbs--lacI-lva-terminator');
iRNA_gamS = findspecies(Mobj, 'RNA rbs--gamS');
plot(t_ode/60, x_ode(:, iDNA_tetR), 'k-', ...
    t_ode/60, x_ode(:, iDNA_lacI), 'b-', ...
  t_ode/60, x_ode(:, iDNA_gamS), 'r-', ...
  t_ode/60, x_ode(:, iRNA_tetR), 'k--', ...
    t_ode/60, x_ode(:, iRNA_tetR), 'b--', ...
  t_ode/60, x_ode(:, iRNA_gamS), 'r--');

title('DNA and mRNA');
lgh = legend({'DNA tetR','DNA lacI','DNA gamS', 'RNA tetR', 'RNA lacI', 'RNA gamS'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');


%% diagnostics

% interested in these
ilacI = findspecies(Mobj, 'protein lacI-lva-terminator');
itetR = findspecies(Mobj, 'protein tetR-lva-terminator');
ilacIdimer = findspecies(Mobj, 'protein lacI-lva-terminatordimer');
itetRdimer = findspecies(Mobj, 'protein tetR-lva-terminatordimer');
iDNA_tetR = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator');
iDNA_lacI = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator');
iDNA_lacI_tetRbounddimer1 = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer1');
iDNA_tetR_lacIbounddimer = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator:protein lacI-lva-terminatordimer1');

% not interested in these for now:
iLacItetramer = findspecies(Mobj, 'protein lacI-lva-terminatortetramer');
iDNA_lacI_tetRbounddimer2 = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer2');
iDNA_lacI_2tetRbound = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer');
iDNA_tetR_2lacIbound = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator:protein lacI-lva-terminatortetramer');

figure(4)
subplot(2,1,1)
plot(t_ode/60, x_ode(:, itetR), 'k-', ...
     t_ode/60, x_ode(:, ilacI), 'b-', ...
     t_ode/60, x_ode(:, itetRdimer), 'k--', ...
     t_ode/60, x_ode(:, ilacIdimer), 'b--');

title('Protein Conc');
lgh = legend(...
  {'tetR','lacI','tetR dimer','lacI dimer'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

subplot(2,1,2)
plot(t_ode/60, x_ode(:, iDNA_tetR), 'k-', ...
     t_ode/60, x_ode(:, iDNA_lacI), 'b-', ...
     t_ode/60, x_ode(:, iDNA_tetR_lacIbounddimer), 'k--', ...     
     t_ode/60, x_ode(:, iDNA_lacI_tetRbounddimer1), 'b--')

title('DNA bound to repressors');
lgh = legend({'DNA tetR', 'DNA lacI', 'DNA tetR:lacIdimer','DNA lacI:tetRdimer'},...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

figure(2)
iaTc = findspecies(Mobj, 'aTc');
iIPTG = findspecies(Mobj, 'IPTG');
plot(t_ode/60, x_ode(:, iaTc), 'k-', ...
    t_ode/60, x_ode(:, iIPTG), 'b-')

title('inducers');
lgh = legend({'aTc', 'IPTG'},...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
