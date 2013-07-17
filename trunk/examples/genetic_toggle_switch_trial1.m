% Genetic Toggle Switch example
% Vipul Singhal, September 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch. It shows the bistability of the circuit, and the
% switching that can be accomplished by the tetR repressing inducer aTc and
% the lacI repressing inducer IPTG. 

% the dimerization and tetramerization rates of lacI and tetR here are
% high. need to check what the actual values are, since that is intimitaly
% related to the working of this switch. 
% another thing to be played with (and hence should be told to anyone who
% uses this simulation, is that TX and TL rates might be mo

%
close all
% it seems that 'standard parts' may be hard, and there may need to exist a
% whole range of strengths for each part. like tetR dimerization rates. 
ptet_DNA = 1;
placI_DNA = 0.5;
initial_tetR = [25];%:10:40
initial_lacI = [15];%:100:400
[a,b] = meshgrid(initial_tetR,initial_lacI);
c = [reshape(a, numel(a), 1) reshape(b, numel(b), 1)]
% Set up the standard TXTL tubes
tube1 = txtl_extract('E9_toggle_switch'); % increased transcription and TL rate, compared to E9, so that we can fit in some purturbations within the 6 hour window. 
tube2 = txtl_buffer('E9_toggle_switch');
Mobj = cell(numel(a),1);
simData = cell(numel(a),1);
x_ode = cell(numel(a),1);
t_ode = cell(numel(a),1);

for i = 1:numel(a)
% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_lacI = txtl_add_dna(tube3,'ptet2(50)', 'rbs(20)', 'lacI2-lva(647)', ptet_DNA, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'placI2(50)', 'rbs(20)', 'tetR2-lva(647)', placI_DNA, 'plasmid');
% dna_deGFP = txtl_add_dna(tube3, 'ptet2(50)', 'rbs(20)', 'deGFP(1000)', ptet_DNA, 'plasmid');
% dna_deCFP = txtl_add_dna(tube3, 'placI2(50)', 'rbs(20)', 'deCFP(1000)', placI_DNA, 'plasmid');
dna_ClpX = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'ClpX', 0, 'plasmid');
% dna_lacI = txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'lacI-lva(647)', 0, 'plasmid');
% dna_tetR = txtl_add_dna(tube3, '70(50)', 'rbs(20)', 'tetR-lva(647)', 0, 'plasmid');

% Mix the contents of the individual tubes
Mobj{i} = txtl_combine([tube1, tube2, tube3]);

%add various species
txtl_addspecies(Mobj{i}, 'protein tetR2-lva',c(i,1));
txtl_addspecies(Mobj{i}, 'protein lacI2-lva',c(i,2));
txtl_addspecies(Mobj{i}, 'aTc',0);
txtl_addspecies(Mobj{i}, 'IPTG',0);
txtl_addspecies(Mobj{i}, 'protein ClpX*',500);

simulationTime = 14*60*60;
tic
simData{i} = txtl_runsim(Mobj{i},simulationTime);
toc


t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;

end

plotGeneticToggle