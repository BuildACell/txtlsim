% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
%

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
1*4.2, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

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
txtl_addspecies(Mobj, 'protein tetR', 10)
% 1st run


 
% txtl_runsim(Mobj);
% set(Mobj.Reactions(6).KineticLaw.Parameters(2),'Value',54261*log(2)/0.1) 
 
[simData] = txtl_runsim(Mobj,configsetObj);
t_ode = simData.Time;
x_ode = simData.Data;



% configsetObj2 = getconfigset(Mobj, 'active');
% simulationTime = 12*60*60;
% set(configsetObj2, 'SolverType', 'ode23s');
% set(configsetObj2, 'StopTime', simulationTime);
% 
% x_ode_2(end,16) = 4.2*2;
% [t_ode,x_ode] = txtl_runsim(Mobj,configsetObj2,t_ode_2,x_ode_2);

% t_ode = simData2.Time;
% x_ode = simData2.Data;

% simDataResp = resample(simData2,[0:10:t_ode(end)])
% dntp = diff(simDataResp.Data(:,10))
% %%
% 
% figure(4)
% plotyy(simDataResp.Time(1:end-1)/60,dntp,simDataResp.Time/60,simDataResp.Data(:,10))


%% plot the result

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'}; 
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
dataGroups{2,2} = {'protein deGFP*','[protein deGFP]_tot'};
dataGroups{2,3} = {'g','g--','r-','g--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,Mobj,dataGroups);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
