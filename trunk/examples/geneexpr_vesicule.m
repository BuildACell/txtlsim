% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E15');
tube2 = txtl_buffer('E15');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
Mobj.UserData.Vesicule = 1;



%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

tic
txtl_runsim(Mobj);
toc


Outside = addcompartment(Mobj, 'outside');

parameters = {'TXTL_membrane transport_F',0.01};

txtl_addspecies(Mobj,'outside.NTP',10); 
txtl_addspecies(Mobj,'outside.ATP',10); 
txtl_addspecies(Mobj,'outside.AA',10); 
txtl_addreaction(Mobj,'outside.NTP -> contents.NTP','MassAction',parameters,'membrane_transport');
txtl_addreaction(Mobj,'outside.ATP -> contents.ATP','MassAction',parameters,'membrane_transport');
txtl_addreaction(Mobj,'outside.AA -> contents.AA','MassAction',parameters,'membrane_transport');




tic
[simData] = txtl_runsim(Mobj,36*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

txtl_plot(simData,Mobj);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
