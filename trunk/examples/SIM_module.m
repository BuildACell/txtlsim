close all

% % 
% % 
% % tube1 = txtl_extract('E22_RNAcascade_6Aug');
% % tube2 = txtl_buffer('E22_RNAcascade_6Aug');
% % 
% % conc004 = 0;
% % conc339 = 4;
% % conc124 = 1;
% % conc007 = 1;
% % conc419 = 0;
% % conc002 = 0;
% % array = conc004;
% % Mobj = cell(length(array), 1);
% % simData = cell(length(array), 1);
% % x_ode = cell(length(array), 1);
% % t_ode = cell(length(array), 1);
% % for i = 1:length(array)
% % tube4 = txtl_newtube('RNA_repression');
% % % txtl_add_dna(tube4, ...
% % %   'pJ23119(35)', 'anti2(91)', 'no_protein', conc004*4.2,'plasmid');
% % % txtl_add_dna(tube4, ...
% % %   'pJ23119(35)', 'att2(287)-anti1(91)', 'no_protein', conc339*4.2,'plasmid');
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'att1(287)-att1(287)-rbs(14)', 'sfGFP(714)', conc124*4.2, 'linear');	
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', conc007*4.2, 'plasmid');
% % % txtl_add_dna(tube4, ...
% % %   'pJ23119(35)', 'att1(287)-rbs(14)', 'RFP(714)', conc419*4.2, 'plasmid');
% % % txtl_add_dna(tube4, ...
% % %   'pJ23119(35)', 'control(91)', 'no_protein', conc002*4.2,'plasmid');
% % 
% % % Mix the contents of the individual tubes
% % Mobj{i} = txtl_combine([tube1, tube2, tube4]);
% % 
% % % Run a simulation
% % simulationTime = 2*60*60;
% % [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
% % t_ode{i} = simData{i}.Time;
% % x_ode{i} = simData{i}.Data;
% % end
% % 
% % Mobj{1}.reactions
% % Mobj{1}.species
% % txtl_plot(simData{1},Mobj{1} )
% % 
% % 
% % cellOfSpecies = {'RNAP70',[];  
% %     'NTP',[]; 
% %     'ATP',[];
% %    'RNA att1-rbs--sfGFP',[];
% %    'AA:ATP:Ribo:RNA att1-rbs--sfGFP','Ribo:RNA att1-rbs--sfGFP';};
% %    
% % 
% % titleString = {'RNAP70'; 'free NTP'; 'free ATP';
% %    'free RNA'; 'Ribobound RNA'};
% % plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, 2)
% % 
% % cellOfSpecies = {'RNase',[];
% %     'RNA att1-rbs--sfGFP:RNase',[];
% %     'RNA att1:RNase',[];
% %    'AA:ATP:Ribo:RNA att1-rbs--sfGFP:RNase','Ribo:RNA att1-rbs--sfGFP:RNase';};
% %    
% % 
% % titleString = {'free RNase'; 'RNA rbs--deGFP:RNase';'RNA att1:RNase'; '(AA:ATP:)Ribo:RNA rbs--deGFP:RNase';};
% % plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, 2)
% 
% 
%% bottom 2 levels
% tube1 = txtl_extract('E22_RNAcascade_6Aug');
% tube2 = txtl_buffer('E22_RNAcascade_6Aug');
% 
% conc004 = [0 0 0 0 ];
% conc339 = [0 20];
% conc124 = 0.5*ones(size(conc339));
% conc007 = 0;
% conc419 = 0.5*ones(size(conc339));
% conc002 = 0;
% array = conc339;
% Mobj = cell(length(array), 1);
% simData = cell(length(array), 1);
% x_ode = cell(length(array), 1);
% t_ode = cell(length(array), 1);
% for i = 1:length(array)
% tube4 = txtl_newtube('RNA_repression');
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'anti2(91)', 'no_protein', conc004(i)*4.2,'plasmid');
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att2(287)-anti1(91)', 'no_protein', conc339(i)*4.2,'plasmid'); % 339
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att1(287)-att1(287)-rbs(14)', 'sfGFP(714)', conc124(i)*4.2, 'plasmid'); % 124	
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', conc007(i)*4.2, 'plasmid');
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att1(287)-rbs(14)', 'RFP(714)', conc419(i)*4.2, 'plasmid'); % 419
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'control(91)', 'no_protein', conc002(i)*4.2,'plasmid');
% 
% % Mix the contents of the individual tubes
% Mobj{i} = txtl_combine([tube1, tube2, tube4]);
% 
% % Run a simulation
% simulationTime = 2*60*60;
% [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
% t_ode{i} = simData{i}.Time;
% x_ode{i} = simData{i}.Data;
% end
% 
% Mobj{1}.reactions
% Mobj{1}.species
% 
% cellOfSpecies ={'protein sfGFP*';
% 'protein RFP*'};
% 
% titleString = {'protein sfGFP*';
% 'protein RFP*'};
% 
% legends = {[num2str(array(1)) 'nM Anti'],[num2str(array(2)) 'nM Anti']};%,num2str(array(5)),num2str(array(3)),num2str(array(4))
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% titleString = {'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att1-att1-rbs--sfGFP','Ribo:RNA att1-att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-att1-rbs--sfGFP'; 
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1_SIM',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1-att1',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att2-anti1:RNA att1-att1',[],[]};
% 
% titleString = {'full RNA att1 att1 GFP';
% '1st nascent';
% '1st att rep';
% '1st att term';
% '2nd nascent';
% 'two att 1 anti rep';
% 'two att 1 anti term';
% 'two att 2 anti rep';
% 'two att 2 anti term';};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% 
% cellOfSpecies ={'protein sfGFP*';
% 'protein RFP*'};
% 
% titleString = {'protein sfGFP*';
% 'protein RFP*'};
% 
% legends = {[num2str(array(1)) 'nM Anti'],[num2str(array(2)) 'nM Anti']}%,num2str(array(5)){num2str(array(1)),num2str(array(2)),num2str(array(3)),num2str(array(4))};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% titleString = {'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];
%     'RNAP70',[],[];
%     'DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'RNA att1_SIM',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];};
% 
% 
% 
% 
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-att1-rbs--sfGFP';
%     'RNA att1_SIM';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM'};
%     
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];
%         'RNAP70',[],[];
%     'DNA pJ23119--att1-rbs--RFP',[],[];
%     'RNA att1',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];};
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-rbs--RFP';
%     'RNA att1';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1'};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% iGFP = findspecies(Mobj{1}, 'protein sfGFP*');
% iRFP = findspecies(Mobj{1}, 'protein RFP*');
% figure
% plot(t_ode{1}, x_ode{1}(:,iGFP),'b', t_ode{1}, x_ode{1}(:,iRFP),'r')
% legend('GFP', 'RFP')
% % 
%% all 3 levels
tube1 = txtl_extract('E22_RNAcascade_6Aug');
tube2 = txtl_buffer('E22_RNAcascade_6Aug');

conc004 = [0  2  8];
conc339 = 4*ones(size(conc004));
conc124 = 0.5*ones(size(conc004));
conc007 = 0;
conc419 = 0.5*ones(size(conc004));
conc002 = 0;
array = conc004;
Mobj = cell(length(array), 1);
simData = cell(length(array), 1);
x_ode = cell(length(array), 1);
t_ode = cell(length(array), 1);
for i = 1:length(array)
tube4 = txtl_newtube('RNA_repression');
txtl_add_dna(tube4, ...
  'pJ23119(35)', 'anti2(91)', 'no_protein', conc004(i)*4.2,'plasmid');
txtl_add_dna(tube4, ...
  'pJ23119(35)', 'att2(287)-anti1(91)', 'no_protein', conc339(i)*4.2,'plasmid'); % 339
txtl_add_dna(tube4, ...
  'pJ23119(35)', 'att1(287)-att1(287)-rbs(14)', 'sfGFP(714)', conc124(i)*4.2, 'plasmid'); % 124	
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', conc007(i)*4.2, 'plasmid');
txtl_add_dna(tube4, ...
  'pJ23119(35)', 'att1(287)-rbs(14)', 'RFP(714)', conc419(i)*4.2, 'plasmid'); % 419
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'control(91)', 'no_protein', conc002(i)*4.2,'plasmid');

% Mix the contents of the individual tubes
Mobj{i} = txtl_combine([tube1, tube2, tube4]);

% Run a simulation
simulationTime = 2*60*60;
[simData{i}] = txtl_runsim(Mobj{i},simulationTime);
t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;
end


% 
% cellOfSpecies ={'protein sfGFP*';
% 'protein RFP*'};
% 
% titleString = {'protein sfGFP*';
% 'protein RFP*'};
% 
% legends = {[num2str(array(1)) 'nM Anti'],[num2str(array(2)) 'nM Anti'],[num2str(array(3)) 'nM Anti']};%,[num2str(array(4)) 'nM Anti'],[num2str(array(5)) 'nM Anti'],num2str(array(5))
% outputCurves = plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% titleString = {'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att1-att1-rbs--sfGFP','Ribo:RNA att1-att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-att1-rbs--sfGFP'; 
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1_SIM',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1-att1',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att2-anti1:RNA att1-att1',[],[]};
% 
% titleString = {'full RNA att1 att1 GFP';
% '1st nascent';
% '1st att rep';
% '1st att term';
% '2nd nascent';
% 'two att 1 anti rep';
% 'two att 1 anti term';
% 'two att 2 anti rep';
% 'two att 2 anti term';};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)





% cellOfSpecies ={'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% titleString = {'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)

% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];
%     'RNAP70',[],[];
%     'DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'RNA att1_SIM',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];};
% 
% 
% 
% 
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-att1-rbs--sfGFP';
%     'RNA att1_SIM';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM'};
%     
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];
%         'RNAP70',[],[];
%     'DNA pJ23119--att1-rbs--RFP',[],[];
%     'RNA att1',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];};
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-rbs--RFP';
%     'RNA att1';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1'};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)

% iGFP = findspecies(Mobj{4}, 'protein sfGFP*');
% iRFP = findspecies(Mobj{4}, 'protein RFP*');
% figure
% plot(t_ode{4}, x_ode{4}(:,iGFP),'b', t_ode{4}, x_ode{4}(:,iRFP),'r')
% legend('GFP', 'RFP')

% %% spike in experiment
% 
% tube1 = txtl_extract('E22_RNAcascade_6Aug');
% tube2 = txtl_buffer('E22_RNAcascade_6Aug');
% 
% conc004 = 0*ones(size(conc004));
% conc339 = 2*ones(size(conc004));
% conc124 = 0.25*ones(size(conc004));
% conc007 = 0;
% conc419 = 0.25*ones(size(conc004));
% conc002 = 0;
% array = conc339;
% Mobj = cell(length(array), 1);
% simData = cell(length(array), 1);
% x_ode = cell(length(array), 1);
% t_ode = cell(length(array), 1);
% for i = 1:length(array)
% tube4 = txtl_newtube('RNA_repression');
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'anti2(91)', 'no_protein', conc004(i)*4.2,'plasmid');
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att2(287)-anti1(91)', 'no_protein', conc339(i)*4.2,'plasmid'); % 339
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att1(287)-att1(287)-rbs(14)', 'sfGFP(714)', conc124(i)*4.2, 'plasmid'); % 124	
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'att1(287)-rbs(14)', 'sfGFP(714)', conc007(i)*4.2, 'plasmid');
% txtl_add_dna(tube4, ...
%   'pJ23119(35)', 'att1(287)-rbs(14)', 'RFP(714)', conc419(i)*4.2, 'plasmid'); % 419
% % txtl_add_dna(tube4, ...
% %   'pJ23119(35)', 'control(91)', 'no_protein', conc002(i)*4.2,'plasmid');
% 
% % Mix the contents of the individual tubes
% Mobj{i} = txtl_combine([tube1, tube2, tube4]);
% 
% % Run a simulation
% simulationTime = 1*60*60;
% [simData{i}] = txtl_runsim(Mobj{i},simulationTime);
% t_ode{i} = simData{i}.Time;
% x_ode{i} = simData{i}.Data;
% 
% iDNA = findspecies(Mobj{i}, 'DNA pJ23119--anti2');
% x_ode{i}(end,iDNA) = x_ode{i}(end,iDNA) + 20;
% simulationTime = 1*60*60;
% [ t_ode{i},x_ode{i}] = txtl_runsim(Mobj{i},simulationTime, t_ode{i},x_ode{i});
% 
% 
% 
% end
% 
% Mobj{1}.reactions
% Mobj{1}.species
% 
% cellOfSpecies ={'protein sfGFP*';
% 'protein RFP*'};
% 
% titleString = {'protein sfGFP*';
% 'protein RFP*'};
% 
% legends = {num2str(array(1)),num2str(array(2)),num2str(array(3)),num2str(array(4))};%,num2str(array(5))
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% titleString = {'RNA att2-anti1';
% 'RNA att1-att1-rbs--sfGFP';
% 'RNA att1-rbs--RFP'};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNA att1-att1-rbs--sfGFP','Ribo:RNA att1-att1-rbs--sfGFP','AA:ATP:Ribo:RNA att1-att1-rbs--sfGFP'; 
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1_SIM',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1','NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1',[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att1-att1',[],[];
%     'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1-att1:RNA att2-anti1:RNA att2-anti1',[],[];
%     'RNA att2-anti1:RNA att2-anti1:RNA att1-att1',[],[]};
% 
% titleString = {'full RNA att1 att1 GFP';
% '1st nascent';
% '1st att rep';
% '1st att term';
% '2nd nascent';
% 'two att 1 anti rep';
% 'two att 1 anti term';
% 'two att 2 anti rep';
% 'two att 2 anti term';};
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNA att1-rbs--RFP','Ribo:RNA att1-rbs--RFP','AA:ATP:Ribo:RNA att1-rbs--RFP'; 
%     'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1','NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[];
% 'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1:RNA att2-anti1',[],[];
% 'RNA att2-anti1:RNA att1',[],[];};
% 
% titleString = {'full RNA att1 RFP';
% '1st nascent'
% '1st att rep';
% '1st att term';};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% 
% cellOfSpecies ={'protein sfGFP*';
% 'protein RFP*'};
% 
% titleString = {'protein sfGFP*';
% 'protein RFP*'};
% 
% legends = {num2str(array(1)),num2str(array(2)),num2str(array(3)),num2str(array(4))};%,num2str(array(5))
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% % cellOfSpecies ={'RNA att2-anti1';
% % 'RNA att1-att1-rbs--sfGFP';
% % 'RNA att1-rbs--RFP'};
% % 
% % titleString = {'RNA att2-anti1';
% % 'RNA att1-att1-rbs--sfGFP';
% % 'RNA att1-rbs--RFP'};
% % 
% % plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];
%     'RNAP70',[],[];
%     'DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'RNA att1_SIM',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM',[],[];};
% 
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-att1-rbs--sfGFP';
%     'RNA att1_SIM';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP';
%     'NTP:RNAP70:DNA pJ23119--att1-att1-rbs--sfGFP:RNA att1_SIM'};
%     
% 
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% cellOfSpecies ={'RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];
%         'RNAP70',[],[];
%     'DNA pJ23119--att1-rbs--RFP',[],[];
%     'RNA att1',[],[];
%     'RNA att2-anti1',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP',[],[];
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1',[],[];};
% 
% titleString = {'1st nascent';
%     'RNAP70';
%         'DNA pJ23119--att1-rbs--RFP';
%     'RNA att1';
%     'RNA att2-anti1';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP';
%     'NTP:RNAP70:DNA pJ23119--att1-rbs--RFP:RNA att1'};
% plotCustomSpecies(Mobj,x_ode, t_ode, cellOfSpecies, titleString, legends)
% 
% iGFP = findspecies(Mobj{4}, 'protein sfGFP*');
% iRFP = findspecies(Mobj{4}, 'protein RFP*');
% figure
% plot(t_ode{4}, x_ode{4}(:,iGFP),'b', t_ode{4}, x_ode{4}(:,iRFP),'r')
% legend('GFP', 'RFP')
% 
% 
% iDNA = findspecies(Mobj{4}, 'DNA pJ23119--anti2');
% iRNA = findspecies(Mobj{4}, 'RNA anti2');
% figure
% subplot(2,1,1)
% plot(t_ode{4}, x_ode{4}(:,iDNA),'b')
% legend('DNA anti2')
% 
% subplot(2,1,2)
% plot(t_ode{4}, x_ode{4}(:,iRNA),'r')
% legend('RNA anti2')



