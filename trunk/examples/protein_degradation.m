% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E7');
tube2 = txtl_buffer('E7');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('protein_deg');

% Define the DNA strands (defines TX-TL species + reactions)
dna_clpx = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'ClpX(1269)', ...	% promoter, rbs, gene
 1, ...					% concentration (nM)
  'plasmid');					% type

dna_clpx = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'ClpP(1269)', ...	% promoter, rbs, gene
 0, ...					% concentration (nM)
  'plasmid');					% type

dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP-lva(1000)', ...	% promoter, rbs, gene
 0, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

txtl_addspecies(Mobj,'protein deGFP-lva*',3300, 'Internal');
txtl_addspecies(Mobj,'protein ClpX*',130, 'Internal');

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation

simulationTime = 12*60*60;

% 1st run
[t_ode,x_ode] = txtl_runsim(Mobj,simulationTime);

%% plot the result
close all
% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'}; 
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein ClpX*','protein ClpP*'};
dataGroups{2,3} = {'g','g--','r-','g--','b-.','k'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,Mobj,dataGroups);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
