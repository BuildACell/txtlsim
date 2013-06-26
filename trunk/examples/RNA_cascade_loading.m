function [Mobj, t_ode, x_ode] = RNA_cascade_loading

attconc = 1; anticonc = 0; controlconc = [0 0.5 1 1.5 2 2.5 4 6];%
numParam = length(controlconc);
%colours = {'b', 'r', 'g', 'c', 'm', 'y'
colororder1 = lines;
colororder2 = [0 0 1;
    1 0 0;
    0 1 0;
    160/255 32/255 240/255;
    0 1 1;
    1 69/255 0;
    112/255 138/255 144/255;
    188/255 143/255 143/255];

t_ode = cell(numParam, 1); x_ode = cell(numParam, 1); Mobj = cell(numParam, 1);
j = 1; k = 1;
for i = 1:numParam
% Run the sim
[t_ode{i}, x_ode{i}, Mobj{i}] = RNA_repression(attconc(k), anticonc(j), controlconc(i));
end
folderdate = datestr(now,'yyyymmmmdd_HHMMSS');
mkdir([pwd '\examples\RNA_circuit_modeling\' folderdate '\detailed data'])
dirstr = pwd;
cd([pwd '\examples\RNA_circuit_modeling\' folderdate '\detailed data'])

for i = 1:numParam

%% Species indices
%proteins
iPROT_GFP_MATURE = findspecies(Mobj{i}, 'protein deGFP*');
iPROT_GFP = findspecies(Mobj{i}, 'protein deGFP');
iPROT_SIGMA70  = findspecies(Mobj{i}, 'protein sigma70');
iPROT_GAMS  = findspecies(Mobj{i}, 'protein gamS');

%RNAs, and things bound
iRNA_GFP = findspecies(Mobj{i}, 'RNA att-rbs--deGFP');
iRIBOBOUND_GFP = findspecies(Mobj{i}, 'Ribo:RNA att-rbs--deGFP'); 
iAA_RIBOBOUND_GFP = findspecies(Mobj{i}, 'AA:ATP:Ribo:RNA att-rbs--deGFP'); 
iRNA_ATT = findspecies(Mobj{i}, 'RNA att');
iCOMPLEX1 = findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att');
iNTP_COMPLEX1 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att');
iCOMPLEX2 =  findspecies(Mobj{i}, 'RNA att:RNA anti');
iRNAPBOUND_COMPLEX2 =  findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att:RNA anti');
iRNA_ANTI = findspecies(Mobj{i}, 'RNA anti');
iRNA_CONTROL = findspecies(Mobj{i}, 'RNA control');

% DNAs, and things bound
iRNAP_DNA_CONTROL = findspecies(Mobj{i}, 'RNAP70:DNA pJ23119--control');
iNTP_RNAP_DNA_CONTROL = findspecies(Mobj{i}, 'NTP:RNAP70:DNA pJ23119--control');
iDNA_GFP = findspecies(Mobj{i}, 'DNA pJ23119--att-rbs--deGFP');
iDNA_ANTI = findspecies(Mobj{i}, 'DNA pJ23119--anti');
iDNA_CONTROL = findspecies(Mobj{i}, 'DNA pJ23119--control');

% resources
iNTP  = findspecies(Mobj{i}, 'NTP');
iAA  = findspecies(Mobj{i}, 'AA');
iRIBO  = findspecies(Mobj{i}, 'Ribo');
iRNAP  = findspecies(Mobj{i}, 'RNAP');
iRNAP70  = findspecies(Mobj{i}, 'RNAP70');
iRECBCD  = findspecies(Mobj{i}, 'RecBCD');
% iRNASE  = findspecies(Mobj{i}, 'RNase');

close all

%% Plot figures
figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iPROT_GFP_MATURE), 'r', ...
    t_ode{i}/60, x_ode{i}(:,iRIBOBOUND_GFP)+x_ode{i}(:,iAA_RIBOBOUND_GFP)+x_ode{i}(:,iRNA_GFP), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iRNA_ATT)+x_ode{i}(:,iCOMPLEX1)+x_ode{i}(:,iNTP_COMPLEX1)+x_ode{i}(:,iCOMPLEX2)+x_ode{i}(:,iRNAPBOUND_COMPLEX2), 'g', ...
    t_ode{i}/60, x_ode{i}(:,iRNA_ANTI)+x_ode{i}(:,iCOMPLEX2)+x_ode{i}(:,iRNAPBOUND_COMPLEX2), 'k')
legend('protein deGFP*', 'TOTAL RNA att-rbs--deGFP', 'TOTAL RNA att', 'TOTAL RNA anti','Location','NorthEastOutside')
xlabel('time (min)')
title(['protein expression, and total RNAs att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r100',[num2str(i) '_' folderdate '_protein_expression_and_total_RNAs'])

figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iRIBOBOUND_GFP), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iAA_RIBOBOUND_GFP), 'g', ...
    t_ode{i}/60, x_ode{i}(:,iRNA_GFP), 'k')
legend('Ribo:RNA att-rbs--deGFP', 'AA:Ribo:RNA att-rbs--deGFP', 'RNA att-rbs--deGFP','Location','NorthEastOutside')
xlabel('time (min)')
title(['GFP translation att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r100',[num2str(i) '_' folderdate 'GFPtranslation'] )

figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iRNA_ATT), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iCOMPLEX1), 'g', ...
    t_ode{i}/60, x_ode{i}(:,iNTP_COMPLEX1), 'g.-', ...
    t_ode{i}/60, x_ode{i}(:,iCOMPLEX2), 'k', ...
    t_ode{i}/60, x_ode{i}(:,iRNAPBOUND_COMPLEX2), 'k.-')
legend('RNA att', 'RNAPbound:RNA att', 'NTP:RNAPbound:RNA att', 'RNA att:RNA anti', 'RNAPbound:RNA att:RNA anti','Location','NorthEastOutside')
xlabel('time (min)')
title(['Attenuator RNA dynamics att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r200',[num2str(i) '_' folderdate 'Attenuator RNA dynamics'])


figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iRNA_ANTI), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iCOMPLEX2), 'k', ...
    t_ode{i}/60, x_ode{i}(:,iRNAPBOUND_COMPLEX2), 'k.-')
legend('RNA anti', 'RNA att:RNA anti', 'RNAPbound:RNA att:RNA anti','Location','NorthEastOutside')
xlabel('time (min)')
title(['antisense RNA dynamics att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r200',[num2str(i) '_' folderdate ' antisense RNA dynamics'])

figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iNTP)/x_ode{i}(1,iNTP), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iAA)/x_ode{i}(1,iAA), 'k', ...
    t_ode{i}/60, x_ode{i}(:,iRNAP70)/x_ode{i}(1,iRNAP70), 'g',...
    t_ode{i}/60, x_ode{i}(:,iRECBCD)/x_ode{i}(1,iRECBCD), 'm', ...
    t_ode{i}/60, x_ode{i}(:,iRIBO)/x_ode{i}(1,iRIBO), 'c')
legend('NTP', 'AA', 'RNAP70', 'RecBCD', 'Ribo', 'RNase','Location','NorthEastOutside')
xlabel('time (min)')
title(['resource usage att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r200',[num2str(i) '_' folderdate ' resource usage'])

figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iDNA_CONTROL), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iRNA_CONTROL), 'k', ...
    t_ode{i}/60, x_ode{i}(:,iRNAP_DNA_CONTROL), 'g',...
    t_ode{i}/60, x_ode{i}(:,iNTP_RNAP_DNA_CONTROL), 'r')
legend('DNA pJ23119--control', 'RNA control', 'RNAP70:DNA pJ23119--control', 'NTP:RNAP70:DNA pJ23119--control','Location','NorthEastOutside')
xlabel('time (min)')
title(['control RNA and DNA att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
 print('-dtiff','-r200',[num2str(i) '_' folderdate 'control RNA and DNA'])

figure('OuterPosition',[100 100 860 600])
plot(t_ode{i}/60, x_ode{i}(:,iDNA_GFP), 'b', ...
    t_ode{i}/60, x_ode{i}(:,iDNA_ANTI), 'k', ...
    t_ode{i}/60, x_ode{i}(:,iDNA_CONTROL), 'g')
legend('DNA pJ23119--att-rbs--deGFP', 'DNA pJ23119--anti', 'DNA pJ23119--control','Location','NorthEastOutside')
xlabel('time (min)')
title(['DNAs att=' num2str(attconc(k)) ' anti=' num2str(anticonc(j)) ' control=' num2str(controlconc(i))])
print('-dtiff','-r200',[num2str(i) '_' folderdate '_DNAs'])

end

cd([dirstr '\examples\RNA_circuit_modeling\' folderdate])
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: protein deGFP*'];
plotAcrossSims(Mobj, t_ode, x_ode, 'protein deGFP*', 'protein_deGFP', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RNA att-rbs--deGFP'];
plotAcrossSims(Mobj, t_ode, x_ode, 'RNA att-rbs--deGFP', 'RNA_att_rbs_deGFP', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RNA control'];
plotAcrossSims(Mobj, t_ode, x_ode, 'RNA control', 'RNA_control', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: DNA pJ23119--control '];
plotAcrossSims(Mobj, t_ode, x_ode, 'DNA pJ23119--control', 'DNA_pJ23119_control', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: NTP'];
plotAcrossSims(Mobj, t_ode, x_ode, 'NTP', 'NTP', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RecBCD'];
plotAcrossSims(Mobj, t_ode, x_ode, 'RecBCD', 'RecBCD', figTitle, folderdate, 'control', controlconc, colororder2)
% figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: protein gamS'];
% plotAcrossSims(Mobj, t_ode, x_ode, 'protein gamS', 'protein_gamS', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RNAP'];
plotAcrossSims(Mobj, t_ode, x_ode, 'RNAP', 'RNAP', figTitle, folderdate, 'control', controlconc, colororder2)
figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RNAP70'];
plotAcrossSims(Mobj, t_ode, x_ode, 'RNAP70', 'RNAP70', figTitle, folderdate, 'control', controlconc, colororder2)
% figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: RNAP'];
% plotAcrossSims(Mobj, t_ode, x_ode, 'RNAP', 'RNAP', figTitle, folderdate, 'control', controlconc, colororder2)
speciesToPlot = {'RNAP70:DNA pJ23119--att-rbs--deGFP',...
'NTP:RNAP70:DNA pJ23119--att-rbs--deGFP',...
'RNAP70:DNA pJ23119--control',...
'NTP:RNAP70:DNA pJ23119--control',...
'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att'};
%'RNAP70:DNA pJ23119--att-rbs--deGFP:RNA att:RNA anti',...
%'RNAP70:DNA pJ23119--anti',...
%'NTP:RNAP70:DNA pJ23119--anti',...
for ii = 1:length(speciesToPlot)
    figTitle =  ['att = ' num2str(attconc(k)) '; anti = ' num2str(anticonc(j)) '; species: ' speciesToPlot{ii}];
    plotAcrossSims(Mobj, t_ode, x_ode, speciesToPlot{ii}, ['_fig_' num2str(ii)], figTitle, folderdate, 'control', controlconc, colororder2)
end
% 
% %RNAP70comparisionPlot(Mobj,t_ode, x_ode, 'control', controlconc, colororder2)
% %%
 cd(dirstr)
 end


 



  
  
  