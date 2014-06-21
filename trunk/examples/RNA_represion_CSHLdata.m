
close all
simNumber = 4;
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



% p = zeros(8,1);
% 
% for i = 1:4
%     plotthis = mean(exp015(:,3*(i-1)+1:3*(i-1)+3),2)-exp015(:,13);
%     hold on
%     p(2*(i-1)+1) = plot(0:300:7200,plotthis);
%     set(p(2*(i-1)+1), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
%     
%     iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
%     p(2*(i-1)+2) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
%     set(p(2*(i-1)+2), 'Color', colororder3(i,:), 'LineWidth', 1.5)
%     
% end
% legend(p,'Att1GFP 0.25 Exp','Att1GFP 0.25 Sim','Att1GFP 0.5 Exp',...
%     'Att1GFP 0.5 Sim','Att1GFP 1 Exp','Att1GFP 1 Sim',...
%     'Att1GFP 2 Exp','Att1GFP 2 Sim','Location', 'Northwest')
% 
% title('Att1-GFP expression, Dashed = Experiment, Solid = Simulation')
% xlabel('time/s')
% ylabel('RFU')
% hold off
%{
% Plots 1 - 4 are plots of hidden species.
% Plot 1
% cellOfSpecies = {'protein sfGFP*',[],[];
%
%     'RNA att1-rbs--sfGFP', 'Ribo:RNA att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-rbs--sfGFP';
%
%     'RNase',[],[];};
%
% titleString = {'GFP, titrate GFP DNA';
%     'RNA att-GFP, titrate GFP DNA';
%     'RNase, titrate GFP DNA'} ;
%
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, conc_array)
%
% % Plot 2
% cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
%
%     'RNA att1', [];};
%
% titleString = {'RNA att1-active transcription, titrate GFP DNA';
%     'RNA att1-free, titrate GFP DNA'} ;
%
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, conc_array)
%
% % Plot 3
% cellOfSpecies = {'RNAP70',[],[],[];
%
%     'RNAP70:DNA pJ23119--att1-rbs--sfGFP','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP','RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
%
%     'NTP',[],[],[];
%
%     'ATP',[],[],[];
%
%     'AA',[],[],[];};
%
% titleString = {'RNAP70-free, titrate GFP DNA';
%     'RNAP70-bound, titrate GFP DNA';
%     'NTP-free, titrate GFP DNA';
%     'ATP-free, titrate GFP DNA';
%     'AA-free, titrate GFP DNA'} ;
%
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, conc_array)
%
% % Plot 4
% cellOfSpecies = {'RNAP';
%
%     'Ribo';
%
%     'protein sfGFP';
%
%     'protein sfGFP*';};
%
% titleString = {'RNAP, titrate GFP DNA';
%     'Ribo, titrate GFP DNA';
%     'GFP-unmature, titrate GFP DNA';
%     'GFP-mature, titrate GFP DNA'} ;
%
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, conc_array)
%
%}


%% Single repression: compared to control RNA
% att = [1 1 2 2];
% anti_array= [0 8 0 8];
% ctrl_array= [8 0 8 0];

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations


%% Single repression, att1  (and controls)
% Mobj = cell(length(anti_array), 1);
% simData = cell(length(anti_array), 1);
% x_ode = cell(length(anti_array), 1);
% t_ode = cell(length(anti_array), 1);
% iGFP = cell(length(anti_array), 1);
% for i = 1:length(anti_array)
%     tube4 = txtl_newtube('RNA_repression');
%     
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', ['att' num2str(att(i)) '(287)-rbs(14)'], 'sfGFP(714)', 0.5*4.2, 'plasmid');
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', ['anti' num2str(att(i)) '(91)'], 'no_protein', anti_array(i)*4.2,'plasmid');
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', 'control(91)', 'no_protein', ctrl_array(i)*4.2,'plasmid');
%     
%     % Mix the contents of the individual tubes
%     Mobj{i} = txtl_combine([tube1, tube2, tube4]);
%     
%     % Run a simulation
%     simulationTime = 2*60*60;
%     [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
%     t_ode{i} = simData{i}.Time;
%     x_ode{i} = simData{i}.Data;
% end
% 
% % Plot 0: Compare GFP in experiment and simulation
% p = zeros(8,1);
% figure
% for i = 1:4
%     plotthis = mean(exp014(:,3*(i-1)+1:3*(i-1)+3),2)-exp014(:,13);
%     hold on
%     p(2*(i-1)+1) = plot(0:300:7200,plotthis);
%     set(p(2*(i-1)+1), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
%     
%     iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
%     p(2*(i-1)+2) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
%     set(p(2*(i-1)+2), 'Color', colororder3(i,:), 'LineWidth', 1.5)
%     
% end
% legend(p,'Att1 Ctrl E','Att1 Ctrl S',...
%     'Att1 Anti1 E','Att1 Anti1 S',...
%     'Att2 Ctrl E','Att2 Ctrl S',...
%     'Att2 Anti2 E','Att2 Anti2 S', 'Location', 'Northwest')
% title('Rep: Att1-GFP @0.5nM, Exp and Sim, Anti/Ctrl @8nM')
% xlabel('time/s')
% ylabel('RFU')
% hold off
%{
% %{
% % % Plot 1a
% % cellOfSpecies = {'protein sfGFP*',[],[];
% %
% %     'RNA att1-rbs--sfGFP', 'Ribo:RNA att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-rbs--sfGFP';
% %
% %     'RNase',[],[];};
% %
% % titleString = {'GFP';
% %     'RNA att-GFP';
% %     'RNase'} ;
% %
% % plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS0 Ctrl8','Att1 AS8 Ctrl0'})
% %
% % % Plot 1b
% % cellOfSpecies = {'protein sfGFP*',[],[];
% %
% %     'RNA att2-rbs--sfGFP', 'Ribo:RNA att2-rbs--sfGFP','AA:ATP:Ribo:RNA att2-rbs--sfGFP';
% %
% %     'RNase',[],[];};
% %
% % titleString = {'GFP';
% %     'RNA att-GFP';
% %     'RNase'} ;
% %
% % plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS0 Ctrl8','Att2 AS8 Ctrl0'})
% %
% % % Plot 2a
% % cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
% %
% %     'RNA att1', [];};
% %
% % titleString = {'RNA att1-active transcription';
% %     'RNA att1-free'} ;
% %
% % plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS0 Ctrl8','Att1 AS8 Ctrl0'})
% %
% % % Plot 2b
% % cellOfSpecies = {'RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2';
% %
% %     'RNA att2', [];};
% %
% % titleString = {'RNA att2-active transcription';
% %     'RNA att1-free'} ;
% %
% % plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS0 Ctrl8','Att2 AS8 Ctrl0'})
% %
% % % Plot 3a
% % cellOfSpecies = {'RNAP70',[],[],[];
% %
% %     'RNAP70:DNA pJ23119--att1-rbs--sfGFP','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP','RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
% %
% %     'NTP',[],[],[];
% %
% %     'ATP',[],[],[];
% %
% %     'AA',[],[],[];};
% %
% % titleString = {'RNAP70-free';
% %     'RNAP70-bound';
% %     'NTP-free';
% %     'ATP-free';
% %     'AA-free'} ;
% %
% % plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS0 Ctrl8','Att1 AS8 Ctrl0'})
% %
% % % Plot 3b
% % cellOfSpecies = {'RNAP70',[],[],[];
% %
% %     'RNAP70:DNA pJ23119--att2-rbs--sfGFP','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP','RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2';
% %
% %     'NTP',[],[],[];
% %
% %     'ATP',[],[],[];
% %
% %     'AA',[],[],[];};
% %
% % titleString = {'RNAP70-free';
% %     'RNAP70-bound';
% %     'NTP-free';
% %     'ATP-free';
% %     'AA-free'} ;
% %
% % plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS0 Ctrl8','Att2 AS8 Ctrl0'})
% %
% % % Plot 4a
% % cellOfSpecies = {'RNAP';
% %
% %     'Ribo';
% %
% %     'protein sfGFP';
% %
% %     'protein sfGFP*';};
% %
% % titleString = {'RNAP';
% %     'Ribo';
% %     'GFP-unmature';
% %     'GFP-mature'} ;
% %
% % plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS0 Ctrl8','Att1 AS8 Ctrl0'})
% %
% % % Plot 4b
% % cellOfSpecies = {'RNAP';
% %
% %     'Ribo';
% %
% %     'protein sfGFP';
% %
% %     'protein sfGFP*';};
% %
% % titleString = {'RNAP';
% %     'Ribo';
% %     'GFP-unmature';
% %     'GFP-mature'} ;
% %
% % plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS0 Ctrl8','Att2 AS8 Ctrl0'})
% %
% %}
% 
% %% Cross repression
% att = [1 1 2 2];
% anti_array= [0 8 0 8];
% ctrl_array= [8 0 8 0];
% 
% 
% Mobj = cell(length(anti_array), 1);
% simData = cell(length(anti_array), 1);
% x_ode = cell(length(anti_array), 1);
% t_ode = cell(length(anti_array), 1);
% iGFP = cell(length(anti_array), 1);
% for i = 1:length(anti_array)
%     
%     if att(i) == 1
%         anti = 2;
%     elseif att(i) == 2
%         anti = 1;
%     end
%     
%     tube4 = txtl_newtube('RNA_repression');
%     
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', ['att' num2str(att(i)) '(287)-rbs(14)'], 'sfGFP(714)', 0.5*4.2, 'plasmid');
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', ['anti' num2str(anti) '(91)'], 'no_protein', anti_array(i)*4.2,'plasmid');
%     txtl_add_dna(tube4, ...
%         'pJ23119(35)', 'control(91)', 'no_protein', ctrl_array(i)*4.2,'plasmid');
%     
%     % Mix the contents of the individual tubes
%     Mobj{i} = txtl_combine([tube1, tube2, tube4]);
%     
%     % Run a simulation
%     simulationTime = 2*60*60;
%     [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
%     t_ode{i} = simData{i}.Time;
%     x_ode{i} = simData{i}.Data;
% end
% 
% % Plot 0: Compare GFP in experiment and simulation
% gainOfOldReaderoverNewReader = mean(exp012(floor(end/2):end,2+1)./(mean(exp014(floor(end/2):end,1:3),2)-exp014(floor(end/2):end,13)));
% 
% p = zeros(8,1);
% figure
% for i = 1:4
%     hold on
%     p(2*(i-1)+1) = plot(0:300:7200,exp012(:,i+1)/gainOfOldReaderoverNewReader);
%     set(p(2*(i-1)+1), 'Color', colororder3(i,:), 'LineWidth', 1.5, 'LineStyle', '--')
%     
%     iGFP = findspecies(Mobj{i}, 'protein sfGFP*');
%     p(2*(i-1)+2) = plot(t_ode{i}, x_ode{i}(:,iGFP )/scalefactor);
%     set(p(2*(i-1)+2), 'Color', colororder3(i,:), 'LineWidth', 1.5)
%     
% end
% legend(p,'Att1 Ctrl E','Att1 Ctrl S',...
%     'Att2 Anti1 E','Att2 Anti1 S',...
%     'Att1 Anti2 E','Att1 Anti2 S',...
%     'Att2 Ctrl E','Att2 Ctrl S', 'Location', 'Northwest')
% title('Xtalk: Att1-GFP @0.5nM, Exp and Sim, Anti/Ctrl @8nM')
% xlabel('time/s')
% ylabel('RFU')
% hold off
%}
%{
%
% % Plot 1a
% cellOfSpecies = {'protein sfGFP*',[],[];
%
%     'RNA att1-rbs--sfGFP', 'Ribo:RNA att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-rbs--sfGFP';
%
%     'RNase',[],[];};
%
% titleString = {'GFP';
%     'RNA att-GFP';
%     'RNase'} ;
%
% plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS2-0 Ctrl-8','Att1 AS2-8 Ctrl-0'})
%
% % Plot 1b
% cellOfSpecies = {'protein sfGFP*',[],[];
%
%     'RNA att2-rbs--sfGFP', 'Ribo:RNA att2-rbs--sfGFP','AA:ATP:Ribo:RNA att2-rbs--sfGFP';
%
%     'RNase',[],[];};
%
% titleString = {'GFP';
%     'RNA att-GFP';
%     'RNase'} ;
%
% plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS1-0 Ctrl-8','Att2 AS1-8 Ctrl-0'})
%
% % Plot 2a
% cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
%
%     'RNA att1', [];};
%
% titleString = {'RNA att1-active transcription';
%     'RNA att1-free'} ;
%
% plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS2-0 Ctrl-8','Att1 AS2-8 Ctrl-0'})
%
% % Plot 2b
% cellOfSpecies = {'RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2';
%
%     'RNA att2', [];};
%
% titleString = {'RNA att2-active transcription';
%     'RNA att1-free'} ;
%
% plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS1-0 Ctrl-8','Att2 AS1-8 Ctrl-0'})
%
% % Plot 3a
% cellOfSpecies = {'RNAP70',[],[],[];
%
%     'RNAP70:DNA pJ23119--att1-rbs--sfGFP','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP','RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
%
%     'NTP',[],[],[];
%
%     'ATP',[],[],[];
%
%     'AA',[],[],[];};
%
% titleString = {'RNAP70-free';
%     'RNAP70-bound';
%     'NTP-free';
%     'ATP-free';
%     'AA-free'} ;
%
% plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS2-0 Ctrl-8','Att1 AS2-8 Ctrl-0'})
%
% % Plot 3b
% cellOfSpecies = {'RNAP70',[],[],[];
%
%     'RNAP70:DNA pJ23119--att2-rbs--sfGFP','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP','RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2','NTP:RNAP70:DNA pJ23119--att2-rbs--sfGFP:RNA att2';
%
%     'NTP',[],[],[];
%
%     'ATP',[],[],[];
%
%     'AA',[],[],[];};
%
% titleString = {'RNAP70-free';
%     'RNAP70-bound';
%     'NTP-free';
%     'ATP-free';
%     'AA-free'} ;
%
% plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS1-0 Ctrl-8','Att2 AS1-8 Ctrl-0'})
%
% % Plot 4a
% cellOfSpecies = {'RNAP';
%
%     'Ribo';
%
%     'protein sfGFP';
%
%     'protein sfGFP*';};
%
% titleString = {'RNAP';
%     'Ribo';
%     'GFP-unmature';
%     'GFP-mature'} ;
%
% plotCustomSpecies(Mobj(1:2), x_ode(1:2), t_ode(1:2), cellOfSpecies, titleString, {'Att1 AS2-0 Ctrl-8','Att1 AS2-8 Ctrl-0'})
%
% % Plot 4b
% cellOfSpecies = {'RNAP';
%
%     'Ribo';
%
%     'protein sfGFP';
%
%     'protein sfGFP*';};
%
% titleString = {'RNAP';
%     'Ribo';
%     'GFP-unmature';
%     'GFP-mature'} ;
%
% plotCustomSpecies(Mobj(3:4), x_ode(3:4), t_ode(3:4), cellOfSpecies, titleString, {'Att2 AS1-0 Ctrl-8','Att2 AS1-8 Ctrl-0'})
%
%}
% %

%% double cascade repression.

anti1_array = [0 0 6 6 6 6];
anti2_array = [0 0 18 10 4 0];
ctrl_array = [0 24 0 4 8 18];

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
        'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', 0.5*4.2, 'plasmid');
    
    dna_anti = txtl_add_dna(tube5, ...
        'pJ23119(35)', 'control(91)', 'no_protein', ctrl_array(i)*4.2,'plasmid');
    
    Mobj{i} = txtl_combine([tube1, tube2, tube5]);
    
    % Run a simulation
    simulationTime = 2*60*60;
    [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
    t_ode{i} = simData{i}.Time;
    x_ode{i} = simData{i}.Data;
end
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

legend(p2,'Att1GFP=0.5',...
    'Att1GFP=0.5 C=24',...
    'Att1GFP=0.5 Att2AS1=6 AS2=18',...
    'Att1GFP=0.5 Att2AS1=6 AS2=14 C=4',...
    'Att1GFP=0.5 Att2AS1=6 AS2=10 C=8',...
    'Att1GFP=0.5 Att2AS1=6 C=18',...
    'Location', 'Northwest')

title('Cascade, Dashed = Exp, Solid = Sim')
xlabel('time/s')
ylabel('RFU')
clear p p2
%{ 
% % plot unbound resources
% LegendList = {'GFP=0.5',...
%     'GFP=0.5 C=24',...
%     'GFP=0.5 Att2AS1=6 AS2=18',...
%     'GFP=0.5 Att2AS1=6 AS2=14 C=4',...
%     'GFP=0.5 Att2AS1=6 AS2=10 C=8',...
%     'GFP=0.5 Att2AS1=6 C=18'}
% anti1_array = [0 0 6 6 6 6];
% anti2_array = [0 0 18 14 10 0];
% ctrl_array = [0 24 0 4 8 18];
% 
% % free resources
% cellOfSpecies = {'RNAP70',[],[];
%     'RNase', [],[];
%     'NTP', [],[];
%     'ATP', [],[];};
% titleString = {'free RNAP70';
%     'Free RNase';
%     'Free NTP';
%     'free ATP'} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% % plot RNAP70 in its various bound states
% cellOfSpecies = {'RNAP70:DNA pJ23119--anti2',   'NTP:RNAP70:DNA pJ23119--anti2';
%     'RNAP70:DNA pJ23119--control', 'NTP:RNAP70:DNA pJ23119--control';};
% titleString = {'RNAPbound (DNA anti2)';
%     'RNAPbound (DNA control)';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)

% 
% cellOfSpecies = {  'RNAP70:DNA pJ23119--att2-anti1',  'NTP:RNAP70:DNA pJ23119--att2-anti1';
%     'RNAP70:DNA pJ23119--att2-anti1:RNA att2',  'NTP:RNAP70:DNA pJ23119--att2-anti1:RNA att2';
%     'RNAP70:DNA pJ23119--att2-anti1:RNA att2:RNA anti2',[];}
% titleString = {'RNAPbound (DNA att2-anti1)';
%     'RNAPbound (DNA att2-anti1):att2';
%     'RNAPbound (DNA att2-anti1):att2:(RNA anti2)';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% cellOfSpecies = {  'RNAP70:DNA pJ23119--att1-rbs--sfGFP',  'NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP';
%     'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1', 'NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1'
%     'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1',[];}
% titleString = { 'RNAPbound (DNA att1-GFP)';
%     'RNAPbound (DNA att1-GFP):att1';
%     'RNAPbound (DNA att1-GFP):att1:(RNA anti1)';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% 
% % plot RNase in its various bound states
% cellOfSpecies = {'RNA att2-anti1:RNase',[],[];
%     'RNA att2:RNase',[],[];  };
% titleString = {'RNasebound att2-AS1';
%     'RNasebound att2';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% cellOfSpecies = { 'RNA att1:RNase',[],[];
%     'Ribo:RNA att1-rbs--sfGFP:RNase','AA:ATP:Ribo:RNA att1-rbs--sfGFP:RNase','RNA att1-rbs--sfGFP:RNase';};
% titleString = {'RNasebound att1';
%     'RNasebound att1-GFP';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% cellOfSpecies = {'RNA control:RNase',[],[];
%     'RNA anti2:RNase',[],[]; };
% titleString = {'RNasebound control';
%     'RNasebound AS2';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% 
% % Plot the RNAs in their various states (functional grouplings, can discuss
% % how best to look at things. this is informed by the level of coarse graining we intend to carry out in the model)
% cellOfSpecies = {'RNA anti2',[],[];
%     'RNA att2-anti1',[],[];
%     'RNA att1-rbs--sfGFP','Ribo:RNA att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-rbs--sfGFP'
%     'RNA control',[],[];};
% titleString = {'RNA AS2';
%     'RNA att2-AS1';
%     'RNA att1-GFP';
%     'RNA ctrl';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% titleString = { ' free RNA att1';
%     'free RNA att2';};
% cellOfSpecies =  {'RNA att1';
%     'RNA att2';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
%     'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1','RNA att2-anti1:RNA att1';
%     'RNAP70:DNA pJ23119--att2-anti1:RNA att2','NTP:RNAP70:DNA pJ23119--att2-anti1:RNA att2';
%     'RNAP70:DNA pJ23119--att2-anti1:RNA att2:RNA anti2','RNA anti2:RNA att2';};
% titleString = {'RNA att1 nascent';
%     'RNA att1 repressed';
%     'RNA att2 nascent';
%     'RNA att2 repressed';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
% 
% 
% cellOfSpecies = {'RNAP';
%     'protein sigma70';};
% titleString = {'RNAP';
%     'sig70';} ;
% plotCustomSpecies(Mobj, x_ode, t_ode, cellOfSpecies, titleString, LegendList)
%} 


% plotEverything
cellOfSpecies = {'RNA anti2','RNA anti2:RNA att2',           'RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2:RNA anti2'
    'RNA att2-anti1',        'RNA att2-anti1:RNA att1',      'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1' 
    'RNA att2-anti1-anti1',  'RNA att2-anti1-anti1:RNA att1','RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1:RNA att2-anti1-anti1'
    'RNA att1-rbs--sfGFP',   'RNA att2',                     'protein sfGFP*'
    'RNA control',           'RNA att1' ,                    'RNase'};
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies) %,[], 'save', ['Sim' num2str(simNumber) 'RNA'],'VariousRNA' 

cellOfSpecies = {'RNA anti2:RNase',              'RNA att2:RNase',   'Ribo:RNA att1-rbs--sfGFP:RNase'          
                     'RNA att2-anti1:RNase',         'RNA att1:RNase',   'RNA att1-rbs--sfGFP:RNase' 
                     'RNA att2-anti1-anti1:RNase',   'RNA control:RNase','AA:AGTP:Ribo:RNA att1-rbs--sfGFP:RNase'};                 
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies)%,[], 'save', ['Sim' num2str(simNumber) 'RNA'],'RNase' 

cellOfSpecies = {'RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1', 'CUTP:AGTP:RNAP70:DNA pJ23119--att1-rbs--sfGFP:RNA att1';
                 'RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2','CUTP:AGTP:RNAP70:DNA pJ23119--att2-anti1-anti1:RNA att2';
                 'AGTP'                                           'CUTP'};
plotCustomSpecies2(Mobj, x_ode, t_ode, cellOfSpecies)% ,[], 'save', ['Sim' num2str(simNumber) 'RNA'],'nascentDNA' 