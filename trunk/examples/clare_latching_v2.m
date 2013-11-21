% clare_latching_v2.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%
close all
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E15');
tube2 = txtl_buffer('E15');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('clare_latching_v2');

% Define the DNA strands (defines TX-TL species + reactions)
dna_AraC = txtl_add_dna(tube3, 'placI(50)', 'anti2(91)-rbs(20)', 'AraC(647)', 1, 'plasmid');
 dna_LuxR = txtl_add_dna(tube3, 'ptet(50)', 'anti1(91)-rbs(20)', 'LuxR(647)', 1, 'plasmid');
 % separate att2-kinase and att2_luxR because thats what Clare is doing experimentally!
dna_att2_LuxR = txtl_add_dna(tube3, 'plux(50)', 'att2(287)-rbs(20)', 'LuxR(647)-lva(20)', 1, 'plasmid'); 
dna_att2_Kinase = txtl_add_dna(tube3, 'plux(50)', 'att2(287)-rbs(20)', 'NRII_L16R(1050)', 1, 'plasmid'); %NRII_L16R is the kinase
% 
 dna_att1_AraC = txtl_add_dna(tube3, 'pBAD(50)', 'att1(287)-rbs(20)', 'AraC(647)-lva(20)', 1, 'plasmid'); 
 dna_att1_NRI = txtl_add_dna(tube3, 'pBAD(50)', 'att1(287)-rbs(20)', 'NRI(1410)', 1, 'plasmid'); 
% 
 dna_deGFP = txtl_add_dna(tube3, 'pGlnA(50)', 'rbs(20)', 'deGFP(1000)', 1, 'plasmid'); 

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

%txtl_plot(simData,Mobj);

dataGroups = txtl_getDefaultPlotDataStruct();
%dataGroups(2).SpeciesToPlot   = {'protein NRI','protein NRI-p'};
txtl_plot(simData,Mobj,dataGroups);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
