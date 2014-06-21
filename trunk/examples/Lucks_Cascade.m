
close all
simNumber = 6;
figure
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

% import data
rawdata_CSHL

%% No Repression
% Corresponds to Experiment 1.015 carried out at CSHL Aug 2013. We titrated
% Att1-GFP at 0.25, 0.5, 1, 2 nM. There was no control DNA.
conc_array= [0.25 0.5 1 2];

Mobj = cell(length(conc_array), 1);
simData = cell(length(conc_array), 1);
x_ode = cell(length(conc_array), 1);
t_ode = cell(length(conc_array), 1);

for i = 1:length(conc_array)
    tube4 = txtl_newtube('RNA_repression');
    txtl_add_dna(tube4, ...
        'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', conc_array(i)*4.2, 'plasmid');
    
    % Mix the contents of the individual tubes
    Mobj{i} = txtl_combine([tube1, tube2, tube4]);
    
    % Run a simulation
    simulationTime = 2*60*60;
    [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
    t_ode{i} = simData{i}.Time;
    x_ode{i} = simData{i}.Data;
end

% Plot 0: This is a plot that compares the experimental data to the
% simulation data. Dashed line is experimental data. Solid line is
% simulation. The 2nM experimental data has been scaled to match the
% simulation levels. all experimental data is scaled by the same factor.
% The purpose is to get the relative repressions to match up.



% Scale the simulation data (assuming the shapes are the same)
iGFP = findspecies(Mobj{4}, 'protein sfGFP*');
GFPsim_interp = interp1(t_ode{4}', x_ode{4}(:,iGFP )', 0:300:7200);
plotthis = mean(exp015(:,3*(4-1)+1:3*(4-1)+3),2)-exp015(:,13);
scalefactor = mean(GFPsim_interp./plotthis');



p = zeros(8,1);
hold on
for i = 1:4
    plotthis = mean(exp015(:,3*(i-1)+1:3*(i-1)+3),2)-exp015(:,13);
    hold on
    p(2*(i-1)+1) = plot(0:300:7200,plotthis);
    set(p(2*(i-1)+1), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
    hold on
    iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
    p(2*(i-1)+2) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
    set(p(2*(i-1)+2), 'Color', colororder3(i,:), 'LineWidth', 1.5)
    
end
legend(p,'Att1GFP 0.25 Exp','Att1GFP 0.25 Sim','Att1GFP 0.5 Exp',...
    'Att1GFP 0.5 Sim','Att1GFP 1 Exp','Att1GFP 1 Sim',...
    'Att1GFP 2 Exp','Att1GFP 2 Sim','Location', 'Northwest')

title('Att1-GFP expression, Dashed = Experiment, Solid = Simulation')
xlabel('time/s')
ylabel('RFU')
hold off



%% Single repression: compared to control RNA
att = [1 1 2 2];
anti_array= [0 8 0 8];
ctrl_array= [8 0 8 0];

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations


%% Single repression, att1  (and controls)
Mobj = cell(length(anti_array), 1);
simData = cell(length(anti_array), 1);
x_ode = cell(length(anti_array), 1);
t_ode = cell(length(anti_array), 1);
iGFP = cell(length(anti_array), 1);
for i = 1:length(anti_array)
    tube4 = txtl_newtube('RNA_repression');
    
    txtl_add_dna(tube4, ...
        'pJ23119(35)', ['att' num2str(att(i)) '(287)-rbs(14)'], 'sfGFP(714)', 0.5*4.2, 'plasmid');
    txtl_add_dna(tube4, ...
        'pJ23119(35)', ['anti' num2str(att(i)) '(91)'], 'no_protein', anti_array(i)*4.2,'plasmid');
    txtl_add_dna(tube4, ...
        'pJ23119(35)', 'control(91)', 'no_protein', ctrl_array(i)*4.2,'plasmid');
    
    % Mix the contents of the individual tubes
    Mobj{i} = txtl_combine([tube1, tube2, tube4]);
    
    % Run a simulation
    simulationTime = 2*60*60;
    [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
    t_ode{i} = simData{i}.Time;
    x_ode{i} = simData{i}.Data;
end

% Plot 0: Compare GFP in experiment and simulation
p = zeros(8,1);
figure
for i = 1:4
    plotthis = mean(exp014(:,3*(i-1)+1:3*(i-1)+3),2)-exp014(:,13);
    hold on
    p(2*(i-1)+1) = plot(0:300:7200,plotthis);
    set(p(2*(i-1)+1), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
    
    iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
    p(2*(i-1)+2) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
    set(p(2*(i-1)+2), 'Color', colororder3(i,:), 'LineWidth', 1.5)
    
end
legend(p,'Att1 Ctrl E','Att1 Ctrl S',...
    'Att1 Anti1 E','Att1 Anti1 S',...
    'Att2 Ctrl E','Att2 Ctrl S',...
    'Att2 Anti2 E','Att2 Anti2 S', 'Location', 'Northwest')
title('Rep: Att1-GFP @0.5nM, Exp and Sim, Anti/Ctrl @8nM')
xlabel('time/s')
ylabel('RFU')
hold off



%% double cascade repression.
anti1max = 2; anti2max = 18; att1GFP = 0.5;
anti1_array = [0        0               anti1max        anti1max        anti1max        anti1max];
anti2_array = [0        0               anti2max        10              4               0];
ctrl_array =  [0    anti1max+anti2max      0           anti2max-6    anti2max-1      anti2max];

Mobj = cell(length(anti1_array), 1);
simData = cell(length(anti1_array), 1);
x_ode = cell(length(anti1_array), 1);
t_ode = cell(length(anti1_array), 1);

for i = 1:length(anti1_array)
    tube5 = txtl_newtube('RNA_repression_3level');
    % DNA with antisense RNA and dummy protein
    dna_anti1 = txtl_add_dna(tube5, ...
        'pJ23119(35)', 'anti2(91)', 'no_protein', anti2_array(i)*4.2,'plasmid');
    
    dna_att1_anti2 = txtl_add_dna(tube5, ...
        'pJ23119(35)', 'att2(287)-anti1(91)-anti1(91)', 'no_protein', anti1_array(i)*4.2, 'plasmid');
    
    dna_att2_deGFP = txtl_add_dna(tube5, ...
        'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', att1GFP*4.2, 'plasmid');
    
    dna_anti = txtl_add_dna(tube5, ...
        'pJ23119(35)', 'control(91)', 'no_protein', ctrl_array(i)*4.2,'plasmid');
    
    Mobj{i} = txtl_combine([tube1, tube2, tube5]);
    
    % Run a simulation
    simulationTime = 2*60*60;
    [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
    t_ode{i} = simData{i}.Time;
    x_ode{i} = simData{i}.Data;
end



% plotEverything
cellOfSpecies = {'RNA anti2','RNA anti2:RNA att2',           'RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2:RNA anti2'
    'RNA att2-anti1',        'RNA att2-anti1:RNA att1',      'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1' 
    'RNA att2-anti1-anti1',  'RNA att2-anti1-anti1:RNA att1','RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1-anti1'
    'RNA att1-rbs--sfGFP',   'RNA att2',                     'protein sfGFP*'
    'RNA control',           'RNA att1' ,                    'RNase'};
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies)

cellOfSpecies = {'RNA anti2:RNase',              'RNA att2:RNase',   'Ribo:RNA att1-rbs--sfGFP:RNase'          
                     'RNA att2-anti1:RNase',         'RNA att1:RNase',   'RNA att1-rbs--sfGFP:RNase' 
                     'RNA att2-anti1-anti1:RNase',   'RNA control:RNase','AA:AGTP:Ribo:RNA att1-rbs--sfGFP:RNase'};                 
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies)

cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1', 'CUTP:AGTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
                 'RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2','CUTP:AGTP:RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2';
                 'AGTP'                                           'CUTP'};
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies) %,[], 'save', ['Sim' num2str(simNumber) 'RNA'],'nascentDNA' 

figure
% Extract out the relevant experimental data: plots 1 - 5 and 9.
expid = [1:5 9]; % experiment ID
for i = 1:6
    plotthis = mean(exp013(:,3*(expid(i)-1)+1:3*(expid(i)-1)+3),2)-exp013(:,31);
    hold on
    p(i) = plot(0:300:7200,plotthis);
    set(p(i), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
    iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
    p2(i) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
    set(p2(i), 'Color', colororder3(i,:), 'LineWidth', 1.5)
end

legend(p2,['Att1GFP=' num2str(att1GFP)],...
    ['Att1GFP=' num2str(att1GFP)                                                                        ' C=' num2str(ctrl_array(2))],...
    ['Att1GFP=' num2str(att1GFP) ' Att2AS1=' num2str(anti1_array(3)) ' AS2=' num2str(anti2_array(3))],...
    ['Att1GFP=' num2str(att1GFP) ' Att2AS1=' num2str(anti1_array(4)) ' AS2=' num2str(anti2_array(4))    ' C=' num2str(ctrl_array(4))],...
    ['Att1GFP=' num2str(att1GFP) ' Att2AS1=' num2str(anti1_array(5)) ' AS2=' num2str(anti2_array(5))    ' C=' num2str(ctrl_array(5))],...
    ['Att1GFP=' num2str(att1GFP) ' Att2AS1=' num2str(anti1_array(6))                                    ' C=' num2str(ctrl_array(6))],...
    'Location', 'Northwest')

title('Cascade, Dashed = Exp, Solid = Sim')
xlabel('time/s')
ylabel('RFU')
clear p p2