% protein_degradation Characterization
% VS Jun 2014

% This file recreates the ClpX characterization that Zach did. It is from
% his powerpont presentaton 10/31/13. There is also on the biocircuits wiki
% Parameter estimation page:
% https://www.cds.caltech.edu/biocircuits/index.php/Modeling_Parameter_Estimation_Summary#ClpX_mediated_protein_Degradation


%% Basic simulation: spike in GFP and ClpX (Zach slide 15 RHS)
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
tube3 = txtl_newtube('protein_deg');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(1269)',0, 'plasmid');	
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)-lva', 0, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
simulationTime = 12*60*60;
txtl_addspecies(Mobj,'protein ClpX*', 200);
txtl_addspecies(Mobj,'protein deGFP-lva*', 4500);
[t,x] = txtl_runsim(Mobj,simulationTime);

% plot the result
close all

dataGroups = txtl_getDefaultPlotDataStruct();
dataGroups(2).SpeciesToPlot   = {'protein deGFP-lva*','protein ClpX*'};

txtl_plot(t,x,Mobj,dataGroups);


%% sweep for spike GFP and ClpX
ct = 1;
clpXconc = [0 12.5 25 50 100 200 400];
t = cell(length(clpXconc),1);
h = zeros(length(clpXconc),1);
x = t;
figure
colororder1 = lines;
for cc = clpXconc
Mobj = txtl_combine([tube1, tube2, tube3]);
simulationTime = 12*60*60;
txtl_addspecies(Mobj,'protein ClpX*', cc);
txtl_addspecies(Mobj,'protein deGFP-lva*', 4500);
[t{ct},x{ct}] = txtl_runsim(Mobj,simulationTime);
iGFP = findspecies(Mobj,'protein deGFP-lva*');
h(ct) = plot(t{ct}/60, x{ct}(:,iGFP));
hold on
set(h(ct), 'color', colororder1(ct,:), 'linewidth', 2);
ct = ct + 1;
end

legendList = {'0 nM', '12.5 nM', '25 nm', '50 nM', '100 nM', '200 nM', '400 nM'};
lgh = legend(h, legendList, 'Location','NorthEastOutside');
title('ClpX Mediated degradation: protein GFP @ 4500nM, protein ClpX vaying','fontsize', 14)
xlabel('time/min','fontsize', 14)
ylabel('deGFP conc/nM','fontsize', 14)
set(gca,'fontsize', 14)


%% Gene expression in the presence of ClpX
clear all
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
tube3 = txtl_newtube('protein_deg');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(1269)',0.001, 'plasmid');	
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)-lva', 30, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
simulationTime = 8*60*60;
% txtl_addspecies(Mobj,'protein ClpX*', 200);
% txtl_addspecies(Mobj,'protein deGFP-lva*', 4500);
[t,x] = txtl_runsim(Mobj,simulationTime);

% plot the result
close all

dataGroups = txtl_getDefaultPlotDataStruct();
dataGroups(2).SpeciesToPlot   = {'protein deGFP-lva*','protein ClpX*'};

txtl_plot(t,x,Mobj,dataGroups);
cellofspecies = {'protein deGFP-lva***','protein ClpX*','protein ClpX';
    'protein deGFP-lva','protein deGFP-lva*:protein ClpX*','protein deGFP-lva*'};
plotCustomSpecies2({Mobj}, {x}, {t}, cellofspecies)

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
