% induction.m - gene expression with an inducer
% R. M. Murray, 11 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a negatively autoregulated gene.  The constants for this example
% come from the Simbiology toolbox example page:
%
%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%
clear all
close all
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('E6');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'ptet(50)', 'rbs(20)', 'tetR(1200)', 15, 'linear');
%dna_gamS = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 0, 'plasmid');
% dna_deGFP = txtl_add_dna(tube3, ...
%   'ptet(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');

%
% Next we have to set up the reactions that describe how the circuit
% works.  Transcription and translation are already included above, so
% we just need to include protein-protein and protein-DNA interactions.
%
% Note that the commands in this section are standard Simbiology commands,
% so you can put anything you want here.
%

% No additional reactions required for this circuit
% tetR-DNA interactions are automatically included in tetR setup

%
% Describe the actual experiment that we want to run.  This includes 
% combining the various tubes and also adding any additional inducers
% or purified proteins that you want to include in the run.
%

% Put in gamS to protect linear DNA
gamS = txtl_addspecies(tube3, 'protein gamS', 100);

% Set up the plot
figure(2); clf();
count = 1;

% Do runs at different inducer levels, linearly spaced
levels = [0 2 5 10 20 40 100 300 1000];
maxtetR = zeros(1, length(levels));
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--'};
% Mix the contents of the individual tubes
  Mobj = txtl_combine([tube1, tube2, tube3]);

for atc = levels 
  % Run a simulation
  configsetObj = getconfigset(Mobj, 'active');
  set(configsetObj, 'SolverType', 'ode23s');
  set(configsetObj, 'StopTime', 6*60*60);
   
  [t_ode, x_ode, mObj, simData] = txtl_runsim(Mobj, configsetObj);
  
  % Plot the time trace
  figure(2); hold on;
  itetR = findspecies(Mobj, 'protein tetR');
  plot(t_ode, x_ode(:, itetR), colors{count});
  labels{count} = [int2str(atc) ' nM aTc'];
  
  % Keep track of the max expression for later plotting
  maxtetR(count) = max(x_ode(:, itetR));
  
  % Add additional inducer for the next run
  if count < size(levels,2)
  inducer = txtl_addspecies(Mobj, 'aTc', levels(count+1)-levels(count));
  count = count + 1;
  end
end

% Label the time trace
title('Time Responses');
lgh = legend(labels, 'Location', 'Best');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Plot the characteristic curve
figure(3); clf();
title('Characteristic Curve');
plot(levels, maxtetR);
xlabel('aTc concentration [nM]');
ylabel('max tetR expression [nM]');

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
