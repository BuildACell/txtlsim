% Genetic Toggle Switch example
% Vipul Singhal, September 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch. It shows the bistability of the circuit, and the
% switching that can be accomplished by the tetR repressing inducer aTc and
% the lacI repressing inducer IPTG. 
%

tetR_initialConc = 10; 
lacI_initialConc = 1; 
simulationDuration = 5*60*60;
inductionTimePoints = [1*60*60 3*60*60];
inductionDuration = 5*60; % seconds 
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

lacIprotein = sbioselect(tube3, 'Name','protein lacI-lva-terminator');
tetRprotein = sbioselect(tube3, 'Name','protein tetR-lva-terminator');
set(lacIprotein, 'InitialAmount', lacI_initialConc);
set(tetRprotein, 'InitialAmount', tetR_initialConc);
txtl_addspecies(tube3, 'aTc',0);
txtl_addspecies(tube3, 'IPTG',0);


% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);
sobj_aTc = sbioselect(Mobj, 'Type', 'species', 'Name', 'aTc');
sobj_IPTG = sbioselect(Mobj, 'Type', 'species', 'Name', 'IPTG');
configsetObj = getconfigset(Mobj, 'active');
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end

% sensitivity analysis
configsetObj.SolverOptions.SensitivityAnalysis = true;
sensitivityOpt = configsetObj.SensitivityAnalysisOptions;
lacIprotein = sbioselect(Mobj, 'Name','protein lacI-lva-terminator');
sensitivityOpt.Outputs = [lacIprotein];
params = sbioselect(Mobj,'Type','parameter','Name',{'TXTL_NTP_RNAP_F','TXTL_TX_rate_RNArbs--lacI-lva-terminator_NTP_consumption','TXTL_TL_rate_proteinlacI-lva-terminator_AA_consumption'});
sensitivityOpt.Inputs = params;
sensitivityOpt.Normalization = 'Full';

% first sim
set(configsetObj, 'StopTime', inductionTimePoints(1))
[t_ode,x_ode, simData, t_sen, x_sen, senOutputs, senInputs] = txtl_runsim_sensitivity(Mobj,configsetObj,[],[],[],[]);

% second sim
set(configsetObj, 'StopTime', inductionDuration)
iaTc = findspecies(Mobj, 'aTc');
x_ode(end,iaTc) = inducerConc;
[t_ode1,x_ode1, simData, t_sen1 x_sen1, senOutputs, senInputs] = txtl_runsim_sensitivity(Mobj,configsetObj, t_ode,x_ode, t_sen, x_sen);

% third sim
set(configsetObj, 'StopTime', inductionTimePoints(2)-inductionTimePoints(1)-inductionDuration)
[t_ode2,x_ode2, simData, t_sen2, x_sen2, senOutputs, senInputs] = txtl_runsim_sensitivity(Mobj,configsetObj,t_ode1,x_ode1,t_sen1,x_sen1);

% fourth sim
set(configsetObj, 'StopTime', inductionDuration)
iIPTG = findspecies(Mobj, 'IPTG');
x_ode2(end,iIPTG) = inducerConc;
[t_ode3,x_ode3,simData, t_sen3, x_sen3, senOutputs, senInputs] = txtl_runsim_sensitivity(Mobj,configsetObj,t_ode2,x_ode2,t_sen2,x_sen2);

% fifth sim
set(configsetObj, 'StopTime', simulationDuration-inductionTimePoints(2)-inductionDuration)
[t_ode4,x_ode4, simData, t_sen4, x_sen4, senOutputs, senInputs] = txtl_runsim_sensitivity(Mobj,configsetObj,t_ode3,x_ode3,t_sen3,x_sen3);


%% Plot results
% Top row: protein and RNA levels
figure(1); clf(); subplot(2,1,1);
ilacI = findspecies(Mobj, 'protein lacI-lva-terminator');
itetR = findspecies(Mobj, 'protein tetR-lva-terminator');
iGamS = findspecies(Mobj, 'protein gamS');
iGFP = findspecies(Mobj, 'protein deGFP');
iGFPs = findspecies(Mobj, 'protein deGFP*');

p = plot(t_ode4/60, x_ode4(:, itetR),'k-', t_ode4/60, x_ode4(:, ilacI), 'b-', t_ode4/60, x_ode4(:, iGamS), 'r-', ...
  t_ode4/60, x_ode4(:, iGFP) + x_ode4(:, iGFPs), 'g--', ...
  t_ode4/60, x_ode4(:, iGFPs), 'g-');

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
  t_ode4/60, x_ode4(:, iAA)/x_ode4(1, iAA), 'b-', ...
  t_ode4/60, x_ode4(:, iNTP)/x_ode4(1, iNTP), 'r-', ...
  t_ode4/60, x_ode4(:, iRNAP)/x_ode4(1, iRNAP), 'b--', ...
  t_ode4/60, x_ode4(:, iRibo)/x_ode4(1, iRibo), 'r--');

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
plot(t_ode4/60, x_ode4(:, iDNA_tetR), 'k-', ...
    t_ode4/60, x_ode4(:, iDNA_lacI), 'b-', ...
  t_ode4/60, x_ode4(:, iDNA_gamS), 'r-', ...
  t_ode4/60, x_ode4(:, iRNA_tetR), 'k--', ...
    t_ode4/60, x_ode4(:, iRNA_tetR), 'b--', ...
  t_ode4/60, x_ode4(:, iRNA_gamS), 'r--');

title('DNA and mRNA');
lgh = legend({'DNA tetR','DNA lacI','DNA gamS', 'RNA tetR', 'RNA lacI', 'RNA gamS'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');


%% Plot sensitivity analysis data
 t_sen4 = squeeze(t_sen4);

 figure(3);
 plot(t_sen4/60,x_sen4(:,:,1), t_sen4/60,x_sen4(:,:,2),t_sen4/60,x_sen4(:,:,3));
 title('Normalized sensitivity of lacIprotein with respect to various parameters');
 xlabel('Time (min)');
 ylabel('Sensitivity');
 legend(senInputs, 'Location', 'NorthEastOutside');
 grid on;

%% diagnostics

% interested in these
ilacI = findspecies(Mobj, 'protein lacI-lva-terminator');
itetR = findspecies(Mobj, 'protein tetR-lva-terminator');
ilacIdimer = findspecies(Mobj, 'protein lacI-lva-terminatordimer');
itetRdimer = findspecies(Mobj, 'protein tetR-lva-terminatordimer');
iDNA_tetR = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator');
iDNA_lacI = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator');
iDNA_lacI_tetRbounddimer1 = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer1');
iDNA_tetR_lacIbounddimer = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator:protein lacI-lva-terminatordimer');

% not interested in these for now:
iLacItetramer = findspecies(Mobj, 'protein lacI-lva-terminatortetramer');
iDNA_lacI_tetRbounddimer2 = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer2');
iDNA_lacI_2tetRbound = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lva-terminator:protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer');
iDNA_tetR_2lacIbound = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lva-terminator:protein lacI-lva-terminatortetramer');

figure(4)
subplot(2,1,1)
plot(t_ode4/60, x_ode4(:, itetR), 'k-', ...
     t_ode4/60, x_ode4(:, ilacI), 'b-', ...
     t_ode4/60, x_ode4(:, itetRdimer), 'k--', ...
     t_ode4/60, x_ode4(:, ilacIdimer), 'b--');

title('Protein Conc');
lgh = legend(...
  {'tetR','lacI','tetR dimer','lacI dimer'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

subplot(2,1,2)
plot(t_ode4/60, x_ode4(:, iDNA_tetR), 'k-', ...
     t_ode4/60, x_ode4(:, iDNA_lacI), 'b-', ...
     t_ode4/60, x_ode4(:, iDNA_tetR_lacIbounddimer), 'k--', ...     
     t_ode4/60, x_ode4(:, iDNA_lacI_tetRbounddimer1), 'b--')

title('DNA bound to repressors');
lgh = legend({'DNA tetR', 'DNA lacI', 'DNA tetR:lacIdimer','DNA lacI:tetRdimer'},...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

figure(2)
iaTc = findspecies(Mobj, 'aTc');
iIPTG = findspecies(Mobj, 'IPTG');
plot(t_ode4/60, x_ode4(:, iaTc), 'k-', ...
    t_ode4/60, x_ode4(:, iIPTG), 'b-')

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
