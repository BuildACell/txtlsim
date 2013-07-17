% Genetic Toggle Switch example
% Vipul Singhal, September 2012
%
close all
clear all
%
SimNumber = 19;
close all 
ptet_DNA = 0.5;
placI_DNA = 0.5;
initial_tetR = [0:2:4];
initial_lacI = [0:2:4];
[a,b] = meshgrid(initial_tetR,initial_lacI);
c = [reshape(a, numel(a), 1) reshape(b, numel(b), 1)]

% Set up the standard TXTL tubes
tube1 = txtl_extract('E9');%_toggle_switch_2
tube2 = txtl_buffer('E9');%_toggle_switch_2
Mobj = cell(numel(a),1);
simData = cell(numel(a),1);
x_ode = cell(numel(a),1);
t_ode = cell(numel(a),1);

for i = 1:numel(a)
% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_lacI3 = txtl_add_dna(tube3,'ptet3(50)', 'rbs(20)', 'lacI3-lva(647)', ptet_DNA, 'plasmid');
dna_tetR3 = txtl_add_dna(tube3, 'placI3(50)', 'rbs(20)', 'tetR3-lva(647)', placI_DNA, 'plasmid');
dna_GFP = txtl_add_dna(tube3,'ptet3(50)', 'rbs(20)', 'deGFP', 0.5, 'plasmid'); 
% if tetR wins, GFP is low
dna_RFP = txtl_add_dna(tube3, 'placI3(50)', 'rbs(20)', 'RFP', 0.5, 'plasmid');
% if lacI wins, RFP is low. 
dna_ClpX = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'ClpX', 0, 'plasmid');
dna_lacI = txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'lacI4-lva(647)', 0, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'tetR4-lva(647)', 0, 'plasmid');

% Mix the contents of the individual tubes
Mobj{i} = txtl_combine([tube1, tube2, tube3]);

%add species
txtl_addspecies(Mobj{i}, 'protein tetR3-lva',c(i,1));
txtl_addspecies(Mobj{i}, 'protein lacI3-lva',c(i,2));
txtl_addspecies(Mobj{i}, 'protein tetR4-lva',max(initial_tetR) - c(i,1));
txtl_addspecies(Mobj{i}, 'protein lacI4-lva',max(initial_lacI) - c(i,2));

txtl_addspecies(Mobj{i}, 'aTc',0);
txtl_addspecies(Mobj{i}, 'IPTG',0);
txtl_addspecies(Mobj{i}, 'protein ClpX*',60);

% simulate
simulationTime = 40*60*60;
tic
simData{i} = txtl_runsim(Mobj{i},simulationTime);
toc


t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;

end

plotGT2_troubleshoot