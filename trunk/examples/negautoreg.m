% negautoreg.m - negative autoregulation example
% R. M. Murray, 8 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a negatively autoregulated gene.  The constants for this example
% come from the Simbiology toolbox example page:
%
%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('negautoreg');


% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'tetR(1200)', 1, 'plasmid');%
% dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)',...
% 'deGFP(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
 txtl_addspecies(Mobj, 'aTc', 5);


% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation

simulationTime = 14*60*60; % 14 hours


tic
[simData] = txtl_runsim(Mobj,simulationTime);
toc
t_ode = simData.Time;
x_ode = simData.Data;


%% plot the result

dataGroups = txtl_getDefaultPlotDataStruct();
dataGroups(2).SpeciesToPlot   = {'ALL_PROTEIN','[protein tetRdimer]_tot'};


txtl_plot(t_ode,x_ode,Mobj,dataGroups);


 
 
% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
