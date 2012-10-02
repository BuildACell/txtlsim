% induction.m - gene expression with an inducer
% R. M. Murray, 11 Sep 2012
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
dna_tetR = txtl_dna(tube3, 'p70(50)', 'rbs(20)', 'tetR(647)', 5, 'linear');
dna_deGFP = txtl_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');

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

% Put in gamS to protect linear DNA
gamS = txtl_addspecies(tube3, 'protein gamS', 100);

% Set up the plot
figure(2); clf();
count = 1;

% Do runs at different inducer levels, linearly spaced
levels = [0 2 5 10 20 40 60 80 100];
maxGFP = zeros(1, length(levels));
for atc = levels 
  % Mix the contents of the individual tubes
  Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);

  % Run a simulation
  configsetObj = getconfigset(Mobj, 'active');
  set(configsetObj, 'StopTime', 6*60*60);
    if ~strcmp(version('-release'),'2012a')
        set(configsetObj, 'SolverType', 'ode23s');
    end
  [t_ode, x_ode, names] = sbiosimulate(Mobj, configsetObj);

  % Plot the time trace
  figure(2); hold on;
  iGFP = findspecies(Mobj, 'protein deGFP*');
  plot(t_ode, x_ode(:, iGFP));
  labels{count} = [int2str(atc) ' nM aTc'];
  
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode(:, iGFP));
  
  % Add additional inducer for the next run
  inducer = txtl_addspecies(tube3, 'aTc', 25);
  count = count + 1;
end

% Label the time trace
title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Plot the characteristic curve
figure(3); clf();
title('Characteristic Curve');
plot(levels, maxGFP);
xlabel('aTc concentration [nM]');
ylabel('max GFP expression [nM]');

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
