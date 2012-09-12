% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('c1');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
  100, ...					% concentration (nM)
  'plasmid');					% type

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
set(configsetObj, 'StopTime', 6*60*60)
[t_ode, x_ode, names] = sbiosimulate(Mobj, configsetObj);

% Get the name of the species we want to plot
iGFP = findspecies(Mobj, 'protein deGFP*');
iRNA = findspecies(Mobj, 'RNA rbs=deGFP');
iNTP = findspecies(Mobj, 'NTP');
iAA  = findspecies(Mobj, 'AA');

% Top row: protein and RNA levels
figure(1); clf(); subplot(2,1,1);
plot(t_ode/60, x_ode(:, iGFP), 'b-', t_ode/60, x_ode(:, iRNA), 'r-')

title('Gene Expression');
lgh = legend(names([iGFP, iRNA]), 'Location', 'Northwest');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Second row: resource limits + internal variables
subplot(2,2,3);
plot(t_ode/60, x_ode(:, iAA), 'b-', t_ode/60, x_ode(:, iNTP), 'r-')

title('Resource usage');
lgh = legend(names([iAA, iNTP]), 'Location', 'West');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
