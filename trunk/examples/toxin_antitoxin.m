%
%clear all
clc
%close all

% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('b1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_adddna(tube3, ...
  'ppard(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
  16, ...					% concentration (nM)
  'plasmid');					% type


dna_parD = txtl_adddna(tube3, ...
  'p70(50)', 'rbs(20)', 'parD(1000)', ...	% promoter, rbs, gene
  16, ...					% concentration (nM)
  'plasmid');					% type

dna_parE = txtl_adddna(tube3, ...
  'p70(50)', 'rbs(20)', 'parE(1000)', ...	% promoter, rbs, gene
  10, ...					% concentration (nM)
  'plasmid');					% type


% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 8, 2]);

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
simulationTime = 14*60*60;
set(configsetObj, 'SolverType', 'ode23s');
set(configsetObj, 'StopTime', simulationTime);

% 1st run
[t_ode,x_ode] = txtl_runsim(Mobj,configsetObj);

%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'#(^DNA (\w+[-=]*)*)'};
dataGroups{1,3} = {'r','g','b','k','m'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','protein parEdimerparDdimer','protein parEdimer','protein parDdimer'};
dataGroups{2,3} = {'g','g--','b-','b--','r-.','k','m','y'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot_gui(t_ode,x_ode,Mobj,dataGroups);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
