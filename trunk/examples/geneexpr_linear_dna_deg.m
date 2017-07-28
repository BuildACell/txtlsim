% Characterize linear DNA based on Z suns data. 
% VS Aug 2013
close all
clear all

tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
initialDNA = [1 2 4 8 16];
Mobj = cell(2*length(initialDNA),1);
simData = cell(2*length(initialDNA),1);
t_ode = cell(2*length(initialDNA),1);
x_ode = cell(2*length(initialDNA),1);

%% Plasmid DNA

for i = 1:length(initialDNA)
tube3 = txtl_newtube('gene_expression');
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', initialDNA(i),'plasmid');
Mobj{i} = txtl_combine([tube1, tube2, tube3]);
[simData{i}] = txtl_runsim(Mobj{i},14*60*60);
t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;
end

cellOfSpecies = {'protein deGFP*', 'protein deGFP'
                 'AA:AGTP:Ribo:RNA rbs--deGFP','Ribo:RNA rbs--deGFP'
                 'CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP', 'RNAP70:DNA p70--rbs--deGFP'
                 'CUTP', 'AGTP'};
plotCustomSpecies2(Mobj(1:length(initialDNA)), x_ode(1:length(initialDNA)),...
    t_ode(1:length(initialDNA)), cellOfSpecies)


%% Linear DNA
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
initialDNA = [1 2 4 8 16];
for i = length(initialDNA)+1:2*length(initialDNA)
tube3 = txtl_newtube('gene_expression');
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', initialDNA(i-length(initialDNA)),'linear');
Mobj{i} = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj{i},'protein gamS',25);
[simData{i}] = txtl_runsim(Mobj{i},10*60*60);
t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;
end
cellOfSpecies = {'protein deGFP*', 'protein deGFP', 'RecBCD:gamS';
                 'AA:AGTP:Ribo:RNA rbs--deGFP','Ribo:RNA rbs--deGFP','protein gamS';
                 'CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP', 'RNAP70:DNA p70--rbs--deGFP','DNA p70--rbs--deGFP:RecBCD';
                 'CUTP','AGTP','RecBCD'};
plotCustomSpecies2(Mobj(length(initialDNA)+1:2*length(initialDNA)),...
    x_ode(length(initialDNA)+1:2*length(initialDNA)), t_ode(length(initialDNA)+1:2*length(initialDNA)), cellOfSpecies)

plasmidGFPendpoints = zeros(1,length(initialDNA));
linearGFPendpoints = zeros(1,length(initialDNA));
for i = 1:length(initialDNA)
    iGFP = findspecies(Mobj{i}, 'protein deGFP*');
    plasmidGFPendpoints(i) = x_ode{i}(end,iGFP);
    iGFP = findspecies(Mobj{length(initialDNA)+i}, 'protein deGFP*');
    linearGFPendpoints(i) = x_ode{length(initialDNA)+i}(end,iGFP);
end

figure
p = plot([0 initialDNA], [0 plasmidGFPendpoints],'bo--',[0 initialDNA], [0 linearGFPendpoints],'rs-.');
legend('plasmid', 'linear GamS','Location','NorthEastOutside')
set(p,'LineWidth',2)
axis([0 initialDNA(end) 0 plasmidGFPendpoints(end)*1.5])
% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
