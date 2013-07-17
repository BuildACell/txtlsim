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
initial_tetRdimer = 0;%:10:40
initial_lacItetramer = 0;%:100:400
[a,b] = meshgrid(initial_tetRdimer,initial_lacItetramer);
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


txtl_addspecies(tube3, 'aTc',0);
txtl_addspecies(tube3, 'IPTG',0);


% Mix the contents of the individual tubes
Mobj{i} = txtl_combine([tube1, tube2, tube3]);
% sobj_aTc = sbioselect(Mobj, 'Type', 'species', 'Name', 'aTc');
% sobj_IPTG = sbioselect(Mobj, 'Type', 'species', 'Name', 'IPTG');
txtl_addspecies(Mobj{i}, 'protein tetR2-lva',c(i,1));
txtl_addspecies(Mobj{i}, 'protein lacI2-lva',c(i,2));
txtl_addspecies(Mobj{i}, 'protein ClpX*',500);
%txtl_addspecies(Mobj{i}, 'protein ClpX',10);
simulationTime = 1*60*60;
tic
simData{i} = txtl_runsim(Mobj{i},simulationTime);
toc


t_ode{i} = simData{i}.Time;
x_ode{i} = simData{i}.Data;

iATC = findspecies(Mobj{i}, 'aTc');
iIPTG = findspecies(Mobj{i}, 'IPTG');

aTc_to_add = 1000;
IPTG_to_add = 0;

x_ode{i}(end,iATC)
x_ode{i}(end,iIPTG)
x_ode{i}(end,iATC) = x_ode{i}(end,iATC) + aTc_to_add;
x_ode{i}(end,iIPTG) = x_ode{i}(end,iIPTG) + IPTG_to_add;
simulationTime = 1*60*60;
[t_ode{i}, x_ode{i}] = txtl_runsim(Mobj{i},simulationTime, t_ode{i}, x_ode{i});


aTc_to_add = 0;
IPTG_to_add = 1000;

x_ode{i}(end,iATC)
x_ode{i}(end,iIPTG)
x_ode{i}(end,iATC) = x_ode{i}(end,iATC) + aTc_to_add;
x_ode{i}(end,iIPTG) = x_ode{i}(end,iIPTG) + IPTG_to_add;
simulationTime = 1*60*60;
[t_ode{i}, x_ode{i}] = txtl_runsim(Mobj{i},simulationTime, t_ode{i}, x_ode{i});

aTc_to_add = 5000;
IPTG_to_add = 0;

x_ode{i}(end,iATC)
x_ode{i}(end,iIPTG)
x_ode{i}(end,iATC) = x_ode{i}(end,iATC) + aTc_to_add;
x_ode{i}(end,iIPTG) = x_ode{i}(end,iIPTG) + IPTG_to_add;
simulationTime = 3*60*60;
[t_ode{i}, x_ode{i}] = txtl_runsim(Mobj{i},simulationTime, t_ode{i}, x_ode{i});

end

%% Plot results
% figure
% hold on
% for i = 1:numel(a)
%     iCFP = findspecies(Mobj{i}, 'protein deCFP*');
%     iGFP = findspecies(Mobj{i}, 'protein deGFP*');
%     plot(x_ode{i}(:,iCFP), x_ode{i}(:,iGFP), 'r', x_ode{i}(end,iCFP), x_ode{i}(end,iGFP), 'b*')
% end
% title('phase plane for CFP and GFP')
% figure
% hold on
% for i = 1:numel(a)
%     itetR = findspecies(Mobj{i}, 'protein tetR2-lva');
%     ilacI = findspecies(Mobj{i}, 'protein lacI2-lva');
%     plot(x_ode{i}(:,itetR), x_ode{i}(:,ilacI), 'r', x_ode{i}(end,itetR), x_ode{i}(end,ilacI), 'b*',x_ode{i}(1,itetR), x_ode{i}(1,ilacI), 'bo')
% end
% title('phase plane for tetR (x axis) and lacI (y axis)')
% 

% figure
%  plot(t_ode, x_ode(:,iCFP), 'c', t_ode, x_ode(:,iGFP), 'g')
%  title('CFP (tetR) and GFP (LacI)')
% figure
% plot(x_ode(:,iCFP), x_ode(:,iGFP), 'r')
% title('CFP (tetR) vs GFP (LacI)')

for i = 1:numel(a)
% plot protein conc
figure
ilacI = findspecies(Mobj{i}, 'protein lacI2-lva');
ilacIdimer = findspecies(Mobj{i}, 'protein lacI2-lvadimer');
ilacItetramer = findspecies(Mobj{i}, 'protein lacI2-lvatetramer');
plot(t_ode{i}/3600, x_ode{i}(:,ilacI), 'k', t_ode{i}/3600, x_ode{i}(:,ilacIdimer), 'k--', t_ode{i}/3600, x_ode{i}(:,ilacItetramer), 'k.-');
legend('lacI', 'lacIdimer', 'lacItetramer')
figure
itetR = findspecies(Mobj{i}, 'protein tetR2-lva');
itetRdimer = findspecies(Mobj{i}, 'protein tetR2-lvadimer');
plot(t_ode{i}/3600, x_ode{i}(:,itetR), 'k', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'k--')
legend('tetR', 'tetRdimer')

%plot resources
figure
iATP = findspecies(Mobj{i}, 'ATP');
iNTP = findspecies(Mobj{i}, 'NTP');
iAA = findspecies(Mobj{i}, 'AA');
iClpX = findspecies(Mobj{i}, 'protein ClpX*');
plot(t_ode{i}/3600, x_ode{i}(:,iATP), 'k', t_ode{i}/3600, x_ode{i}(:,iNTP), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iAA), 'k.-', t_ode{i}/3600, x_ode{i}(:,iClpX), 'r')
legend('ATP', 'NTP', 'AA', 'ClpX')
title('Resource Usage')
figure

itetR_DNA1 = findspecies(Mobj{i}, 'DNA placI2--rbs--tetR2-lva');
itetR_DNA2 = findspecies(Mobj{i}, 'RNAP70:DNA placI2--rbs--tetR2-lva');
itetR_DNA3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA placI2--rbs--tetR2-lva');
itetR_DNA_repressed = findspecies(Mobj{i}, 'DNA placI2--rbs--tetR2-lva:protein lacI2-lvatetramer');

ilacI_DNA1 = findspecies(Mobj{i}, 'DNA ptet2--rbs--lacI2-lva');
ilacI_DNA2 = findspecies(Mobj{i}, 'RNAP70:DNA ptet2--rbs--lacI2-lva');
ilacI_DNA3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA ptet2--rbs--lacI2-lva');
ilacI_DNA_repressed = findspecies(Mobj{i}, 'DNA ptet2--rbs--lacI2-lva:protein tetR2-lvadimer');

plot(t_ode{i}/3600, x_ode{i}(:,itetR_DNA1)+x_ode{i}(:,itetR_DNA2)+x_ode{i}(:,itetR_DNA3), 'k', t_ode{i}/3600, x_ode{i}(:,itetR_DNA_repressed), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,ilacI_DNA1)+x_ode{i}(:,ilacI_DNA2)+x_ode{i}(:,ilacI_DNA3), 'g', t_ode{i}/3600, x_ode{i}(:,ilacI_DNA_repressed), 'g--')
legend('DNA tetR', 'DNA tetR repressed','DNA lacI', 'DNA lacI repressed')
title('DNAs');

figure
iClpX_complex1 = findspecies(Mobj{i}, 'protein tetR2-lva:protein ClpX*');
iClpX_complex2 = findspecies(Mobj{i}, 'protein lacI2-lva:protein ClpX*');
iClpX_unmature = findspecies(Mobj{i}, 'protein ClpX');



plot(t_ode{i}/3600, x_ode{i}(:,iClpX_complex1)+x_ode{i}(:,iClpX_complex2)+x_ode{i}(:,iClpX_unmature)+x_ode{i}(:,iClpX), 'k',...
    t_ode{i}/3600, x_ode{i}(:,iClpX_complex1),'b',t_ode{i}/3600, x_ode{i}(:,iClpX_complex2),'r',...
    t_ode{i}/3600, x_ode{i}(:,iClpX),'g',t_ode{i}/3600, x_ode{i}(:,iClpX_unmature),'c')
legend('all ClpX protein', 'ClpX tetR', 'ClpX lacI', 'ClpX*', 'ClpX')
title('ClpX');
figure
iRNA_lacI1 = findspecies(Mobj{i}, 'RNA rbs--lacI2-lva');
iRNA_lacI2 = findspecies(Mobj{i}, 'Ribo:RNA rbs--lacI2-lva');
iRNA_tetR1 = findspecies(Mobj{i}, 'RNA rbs--tetR2-lva');
iRNA_tetR2 = findspecies(Mobj{i}, 'Ribo:RNA rbs--tetR2-lva');
plot(t_ode{i}/3600, x_ode{i}(:,iRNA_lacI1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNA_lacI2), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iRNA_tetR1), 'r',t_ode{i}/3600, x_ode{i}(:,iRNA_tetR2), 'r--')
legend('lacI RNA', 'ribobound lacI', 'tetR RNA', 'ribobound tetR')
title('RNA');

figure
iRNAP = findspecies(Mobj{i}, 'RNAP');
iRNAP70 = findspecies(Mobj{i}, 'RNAP70');
iRNAP_lacI1 = findspecies(Mobj{i}, 'RNAP70:DNA ptet2--rbs--lacI2-lva');
iRNAP_lacI2 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA ptet2--rbs--lacI2-lva');
iRNAP_tetR1 = findspecies(Mobj{i}, 'RNAP70:DNA placI2--rbs--tetR2-lva');
iRNAP_tetR2 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA placI2--rbs--tetR2-lva');
iRNAP_ClpX1 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA p70--rbs--ClpX');
iRNAP_ClpX2 = findspecies(Mobj{i}, 'RNAP70:DNA p70--rbs--ClpX');
plot(t_ode{i}/3600, x_ode{i}(:,iRNAP), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNAP70), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_lacI1)+x_ode{i}(:,iRNAP_lacI2), 'g',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR1)+x_ode{i}(:,iRNAP_tetR2), 'b',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_ClpX1)+x_ode{i}(:,iRNAP_ClpX2), 'm')
    
legend('RNAP', 'RNAP70', 'RNAP\_lacI', 'RNAP\_tetR', 'RNAP\_ClpX', 'Location', 'NorthEastOutside')
title('RNAP');
% figure
% plot(t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR2), 'k--')
figure
iaTc = findspecies(Mobj{i}, 'aTc');
iaTcbound = findspecies(Mobj{i}, '2 aTc:protein tetR2-lvadimer');
plot(t_ode{i}/3600, x_ode{i}(:,iaTc), 'k-',t_ode{i}/3600, x_ode{i}(:,iaTcbound), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,itetR), 'r', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'r--')
legend('aTc', 'aTc tetRdimer', 'tetR', 'tetRdimer', 'Location', 'NorthEastOutside')
title('aTc');
end

