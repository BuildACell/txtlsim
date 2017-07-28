% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%
close all
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_AraC = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC', 1, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(647)-lva(20)', 1, 'plasmid');
dna_anti1 = txtl_add_dna(tube3, 'pBAD_ptet(50)', 'anti1(91)', 'no_protein', 1, 'plasmid');
dna_att1 = txtl_add_dna(tube3, 'plac(50)', 'att1(287)-anti2(91)', 'no_protein', 1, 'plasmid');
dna_lacI = txtl_add_dna(tube3, ...
  'p70(50)', 'att2(287)-rbs(20)', 'lacI(1000)', ...	% promoter, rbs, gene
   0.5, ...					% concentration (nM)
  'plasmid');					% type
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'att2(287)-rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type
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
[simData] = txtl_runsim(Mobj,2*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

%txtl_plot(simData,Mobj);



% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
