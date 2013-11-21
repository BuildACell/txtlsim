% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression_pr1');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'p70pr1(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
 4.5, ...					% concentration (nM)
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

% Run a simulation
simulationTime = 14*60*60;

% 1st run
[t_ode,x_ode] = txtl_runsim(Mobj,simulationTime);

%% plot the result


txtl_plot(t_ode,x_ode,Mobj);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
