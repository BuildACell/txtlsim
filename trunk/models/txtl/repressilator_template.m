% repressilator.m - repressilator example
%
% This file contains a simple example of setting up a TXTL simulation
% for a negatively autoregulated gene.  The constants for this example
% come from the Simbiology toolbox example page:
%
%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%

%% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations


% Now set up a tube that will contain our DNA




%% Define the DNA strands (defines TX-TL species + reactions)


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

%% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 1, 1.5]);

%
%% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 5*60*60)

% first rum
[t_ode, x_ode, mObj, simData] = txtl_runsim(Mobj, configsetObj,[], []);



%% plot the result
close all
% DNA and mRNA plot

% Gene Expression Plot

% Resource Plot







% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
