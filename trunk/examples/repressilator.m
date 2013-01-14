% repressilator.m - repressilator example
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
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');



% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_adddna(tube3, 'thio-junk(500)-placI(50)', 'rbs(20)', 'tetR(647)-lva(40)-terminator(100)', 5, 'linear');%
dna_lacI = txtl_adddna(tube3, 'thio-junk(500)-plambda(50)', 'rbs(20)', 'lacI(647)-lva(40)-terminator(100)', 5, 'linear');%
dna_lambda = txtl_adddna(tube3, 'thio-junk(500)-ptet(50)', 'rbs(20)', 'lambda(647)-lva(40)-terminator(100)', 5, 'linear');%
dna_deGFP = txtl_adddna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_adddna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');


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

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 1, 1.5]);

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 5*60*60)
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end

[t_ode, x_ode, mObj, simData] = txtl_runsim(Mobj, configsetObj,[], []);
names = simData.DataNames;

[t_ode1, x_ode1, mObj, simData] = txtl_runsim(Mobj, configsetObj,t_ode, x_ode);

%% plot the result
close all
% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'#(^DNA (\w+[-=]*)*)'};
%dataGroups{1,2} = {'DNA p70--rbs--sigma28'};
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--','m-','k-','y--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
%dataGroups{2,2} = {'protein deGFP-lva-terminator*'};
dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.','c-','y--','m-','k-','r-'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';
%
 txtl_plot(t_ode,x_ode,mObj,dataGroups);




% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
