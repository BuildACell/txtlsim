%


close all

tube1 = txtl_extract('E15');
tube2 = txtl_buffer('E15');

tube3 = txtl_newtube('gene_expression');

% sigma factor
txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'sigma54(1434)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'linear');					% type

% NRII_L16R KINASE
txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'NRII_L16R(1050)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'linear');					% type

% NRII_H139N
txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'NRII_H139N(1050)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'linear');					% type

% NRI
txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'NRI(1410)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'linear');					% type

txtl_add_dna(tube3, ...
  'pGlnA(308)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'linear');					% type



Mobj = txtl_combine([tube1, tube2, tube3]);


tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;


%%

dataGroups = txtl_getDefaultPlotDataStruct();
dataGroups(2).SpeciesToPlot   = {'protein NRI','protein NRI-p'};
txtl_plot(simData,Mobj,dataGroups);


