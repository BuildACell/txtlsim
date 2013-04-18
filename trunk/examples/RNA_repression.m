% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9_RNAcascade');
tube2 = txtl_buffer('E9_RNAcascade');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('RNA_repression');

% DNA with attenuator and deGFP
dna_att_deGFP = txtl_add_dna(tube3, ...
  'pJ23119(50)', 'att(1)-rbs(20)', 'deGFP(1000)', ...	% promoter, utr, gene
1*4.2, ...					% concentration (nM)
  'plasmid');					% type


% DNA with antisense RNA and dummy protein
dna_anti_dummyprotein = txtl_add_dna(tube3, ...
  'pJ23119(50)', 'anti(10)', 'dummy', ...	% promoter, utr, gene
0.01*4.2, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
simulationTime = 14*60*60;
set(configsetObj, 'SolverType', 'ode23s');
set(configsetObj, 'StopTime', simulationTime);

[simData] = txtl_runsim(Mobj,configsetObj);
t_ode = simData.Time;
x_ode = simData.Data;

iGFP = findspecies(Mobj, 'protein deGFP*');
%
iRNA = findspecies(Mobj, 'RNA att-rbs--deGFP');
iRIBOBOUND = findspecies(Mobj, 'Ribo:RNA att-rbs--deGFP'); 
iAA_RIBOBOUND = findspecies(Mobj, 'AA:Ribo:RNA att-rbs--deGFP'); 
%
iATT = findspecies(Mobj, 'RNA att');
iCOMPLEX1 = findspecies(Mobj, 'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att');
iNTP_COMPLEX1 = findspecies(Mobj, 'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att');
iCOMPLEX2 =  findspecies(Mobj, 'RNA att:RNA anti');
iRNAPBOUND_COMPLEX2 =  findspecies(Mobj, 'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att:RNA anti');
% 
iANTI = findspecies(Mobj, 'RNA anti');

close all
figure(1)
plot(t_ode, x_ode(:,iGFP), 'r', ...
    t_ode, x_ode(:,iRIBOBOUND)+x_ode(:,iAA_RIBOBOUND)+x_ode(:,iRNA), 'b', ...
    t_ode, x_ode(:,iATT)+x_ode(:,iCOMPLEX1)+x_ode(:,iNTP_COMPLEX1)+x_ode(:,iCOMPLEX2)+x_ode(:,iRNAPBOUND_COMPLEX2), 'g', ...
    t_ode, x_ode(:,iANTI)+x_ode(:,iCOMPLEX2)+x_ode(:,iRNAPBOUND_COMPLEX2), 'k')
legend('protein deGFP*', 'TOTAL RNA att-rbs--deGFP', 'TOTAL RNA att', 'TOTAL RNA anti')

figure(2)
plot(t_ode, x_ode(:,iRIBOBOUND), 'b', ...
    t_ode, x_ode(:,iAA_RIBOBOUND), 'g', ...
    t_ode, x_ode(:,iRNA), 'k')
legend('Ribo:RNA att-rbs--deGFP', 'AA:Ribo:RNA att-rbs--deGFP', 'RNA att-rbs--deGFP')

figure(3)
plot(t_ode, x_ode(:,iATT), 'b', ...
    t_ode, x_ode(:,iCOMPLEX1), 'g', ...
    t_ode, x_ode(:,iNTP_COMPLEX1), 'g.-', ...
    t_ode, x_ode(:,iCOMPLEX2), 'k', ...
    t_ode, x_ode(:,iRNAPBOUND_COMPLEX2), 'k.-')
legend('RNA att', 'RNAPbound:RNA att', 'NTP:RNAPbound:RNA att', 'RNA att:RNA anti', 'RNAPbound:RNA att:RNA anti')

figure(6)
plot(t_ode, x_ode(:,iANTI), 'b', ...
    t_ode, x_ode(:,iCOMPLEX2), 'k', ...
    t_ode, x_ode(:,iRNAPBOUND_COMPLEX2), 'k.-')
legend('RNA anti', 'RNA att:RNA anti', 'RNAPbound:RNA att:RNA anti')



% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
