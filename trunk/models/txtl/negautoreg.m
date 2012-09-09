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
tube1 = txtl_extract('c1');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_dna(tube3, 'ptet(50)', 'rbs(20)', 'tetR(1000)', 100, 'linear');
dna_gamS = txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 10, 'plasmid');

%
% Next we have to set up the reactions that describe how the circuit
% works.  Transcription and translation are already included above, so
% we just need to include protein-protein and protein-DNA interactions.
%
% Note that the commands in this section are standard Simbiology commands,
% so you can put anything you want here.
%

% No additional reactions required for this circuit
% TetR-DNA interactions are automatically included in TetR setup

%
% Describe the actual experiment that we want to run.  This includes 
% combining the various tubes and also adding any additional inducers
% or purified proteins that you want to include in the run.
%

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 1000)
[t_ode, x_ode, names] = sbiosimulate(Mobj, configsetObj);

% Get the name of the species we want to plot
indexlist = findspecies(Mobj, ...
  {'DNA ptet=rbs=tetR', 'RNA rbs=tetR', 'protein tetR', ...
   'RNAP70', 'Ribo', 'protein gamS', 'NTP', 'AA'});
plot(t_ode, x_ode(:,indexlist));
legend(names(indexlist), 'Location', 'NorthEastOutside')
title('Gene Regulation');
xlabel('Time (seconds)');
ylabel('Species Amounts');

%
% Run a set of experiments to explore the effect of inducers
%
%! TODO: write up this section

% Mobj = txtl_combine(tube1, 6, tube2, 2, tube3, 2);
% Add inducer at a given concentration (must be nanomolar)
% txtl_addspecies(Mobj, 'aTc', 50);

