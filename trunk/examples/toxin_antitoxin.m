%
% clear all
clc
%close all

% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E8');
tube2 = txtl_buffer('E8');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'ppard(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
  5, ...					% concentration (nM)
  'plasmid');					% type


dna_parD = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'parD(1000)', ...	% promoter, rbs, gene
  0, ...					% concentration (nM)
  'plasmid');					% type

dna_parE = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'parE(1000)', ...	% promoter, rbs, gene
  0, ...					% concentration (nM)
  'plasmid');					% type


% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

txtl_addspecies(Mobj,'protein parD',80000/50); % parD Anti 80uM 
txtl_addspecies(Mobj,'protein parE',0); % parE Toxin 68 uM

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
simulationTime = 10*60*60;
set(configsetObj, 'StopTime', simulationTime);

% 1st run
[t_ode,x_ode,~,simData] = txtl_runsim(Mobj,configsetObj);

%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'r','g','b','k','m'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','protein parEdimerparDdimer','protein parEdimer','protein parDdimer'};
dataGroups{2,3} = {'g','g--','b-','b--','r-.','k','m','y'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,Mobj,dataGroups);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End: