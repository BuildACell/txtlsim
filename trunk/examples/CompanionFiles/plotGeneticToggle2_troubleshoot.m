
close all
%set(0,'DefaultFigureVisible','off');
set(0,'DefaultFigureVisible','on');
scrsz = get(0,'ScreenSize');
mainDir = pwd;
mkdir([mainDir '\examples\Vipul\Genetic Toggle\2013_manual_parameter_study\Sim' num2str(SimNumber)])
figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
hold on
for i = 1:numel(a)
    itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
    ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
    
%     if norm([x_ode{i}(end,itetR),x_ode{i}(end,ilacI)]-tetR_SS) < norm([x_ode{i}(end,itetR),x_ode{i}(end,ilacI)]-lacI_SS)
    plot(x_ode{i}(:,itetR), x_ode{i}(:,ilacI), 'r', x_ode{i}(end,itetR), x_ode{i}(end,ilacI), 'r*',x_ode{i}(1,itetR), x_ode{i}(1,ilacI), 'ro')
%     else
%      plot(x_ode{i}(:,itetR), x_ode{i}(:,ilacI), 'b', x_ode{i}(end,itetR), x_ode{i}(end,ilacI), 'b*',x_ode{i}(1,itetR), x_ode{i}(1,ilacI), 'bo')   
%     end
end
title('tetR (x axis) and lacI (y axis)')
cd([mainDir '\examples\Vipul\Genetic Toggle\2013_manual_parameter_study'])
print('-dtiff','-r130',['Sim_' num2str(SimNumber) 'PP'])
saveas(gcf, ['Sim_' num2str(SimNumber) '_PP.fig'])
cd(mainDir)
close 

for i = 1:numel(a)
    
    figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])
    subplot(2, 1, 1)
iATP = findspecies(Mobj{i}, 'ATP');
plot(t_ode{i}/3600, x_ode{i}(:,iATP))
legend('ATP', 'Location', 'NorthEastOutside')

subplot(2, 1, 2)
iRNAP = findspecies(Mobj{i}, 'RNAP');
iRNAP70 = findspecies(Mobj{i}, 'RNAP70');
iRNAP_lacI1 = findspecies(Mobj{i}, 'RNAP70:DNA ptet3--rbs--lacI3-lva');
iRNAP_lacI3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA ptet3--rbs--lacI3-lva');
iRNAP_tetR1 = findspecies(Mobj{i}, 'RNAP70:DNA placI3--rbs--tetR3-lva');
iRNAP_tetR3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA placI3--rbs--tetR3-lva');
iRNAP_ClpX1 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA p70--rbs--ClpX');
iRNAP_ClpX2 = findspecies(Mobj{i}, 'RNAP70:DNA p70--rbs--ClpX');
plot(t_ode{i}/3600, x_ode{i}(:,iRNAP), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNAP70), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_lacI1)+x_ode{i}(:,iRNAP_lacI3), 'g',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR1)+x_ode{i}(:,iRNAP_tetR3), 'b',...
    t_ode{i}/3600, x_ode{i}(:,iRNAP_ClpX1)+x_ode{i}(:,iRNAP_ClpX2), 'm')
    
legend('RNAP', 'RNAP70', 'RNAP\_lacI', 'RNAP\_tetR', 'RNAP\_ClpX', 'Location', 'NorthEastOutside')%, 'Location', 'NorthEastOutside'
h = suptitle(['init tetR monomer = ' num2str(c(i,1)) ', init lacI monomer = ' num2str(c(i,2))]);
cd([mainDir '\examples\Vipul\Genetic Toggle\2013_manual_parameter_study\Sim' num2str(SimNumber)])
print('-dtiff','-r130',['Resources_Sim_' num2str(SimNumber) '_' num2str(i)])
saveas(gcf, ['Resources_Sim_' num2str(SimNumber) '_' num2str(i) '.fig'])
cd(mainDir)
close

figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])
subplot(4, 1, 4)
iNTP = findspecies(Mobj{i}, 'NTP');
iAA = findspecies(Mobj{i}, 'AA');
iGFP = findspecies(Mobj{i}, 'protein deGFP*');
iRFP = findspecies(Mobj{i}, 'protein RFP*');
plot(t_ode{i}/3600, x_ode{i}(:,iGFP), t_ode{i}/3600, x_ode{i}(:,iRFP))
% plot(t_ode{i}/3600, x_ode{i}(:,iATP))
% legend('ATP', 'Location', 'NorthEastOutside')
legend('GFP', 'RFP', 'Location', 'NorthEastOutside')

subplot(4, 1, 1)
itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
itetRdimer = findspecies(Mobj{i}, 'protein tetR3-lvadimer');
ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
ilacIdimer = findspecies(Mobj{i}, 'protein lacI3-lvadimer');
ilacItetramer = findspecies(Mobj{i}, 'protein lacI3-lvatetramer');
plot(t_ode{i}/3600, x_ode{i}(:,itetR), 'r', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'r--',...
    t_ode{i}/3600, x_ode{i}(:,ilacI), 'b', t_ode{i}/3600, x_ode{i}(:,ilacIdimer), 'b--',...
    t_ode{i}/3600, x_ode{i}(:,ilacItetramer), 'b-.');
legend('tetR', 'tetRdimer','lacI', 'lacIdimer', 'lacItetramer', 'Location', 'NorthEastOutside')

subplot(4, 1, 2)
itetR_DNA1 = findspecies(Mobj{i}, 'DNA placI3--rbs--tetR3-lva');
itetR_DNA2 = findspecies(Mobj{i}, 'RNAP70:DNA placI3--rbs--tetR3-lva');
itetR_DNA3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA placI3--rbs--tetR3-lva');
itetR_DNA_repressed = findspecies(Mobj{i}, 'DNA placI3--rbs--tetR3-lva:protein lacI3-lvatetramer');
ilacI_DNA1 = findspecies(Mobj{i}, 'DNA ptet3--rbs--lacI3-lva');
ilacI_DNA2 = findspecies(Mobj{i}, 'RNAP70:DNA ptet3--rbs--lacI3-lva');
ilacI_DNA3 = findspecies(Mobj{i}, 'NTP:RNAP70:DNA ptet3--rbs--lacI3-lva');
ilacI_DNA_repressed = findspecies(Mobj{i}, 'DNA ptet3--rbs--lacI3-lva:protein tetR3-lvadimer');

plot(t_ode{i}/3600, x_ode{i}(:,itetR_DNA1)+x_ode{i}(:,itetR_DNA2)+x_ode{i}(:,itetR_DNA3), 'k', t_ode{i}/3600, x_ode{i}(:,itetR_DNA_repressed), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,ilacI_DNA1)+x_ode{i}(:,ilacI_DNA2)+x_ode{i}(:,ilacI_DNA3), 'g', t_ode{i}/3600, x_ode{i}(:,ilacI_DNA_repressed), 'g--')
legend('DNAtetR', 'DNAtetR rep','DNA lacI', 'DNAlacI rep', 'Location', 'NorthEastOutside')


subplot(4, 1, 3)
iClpX_tetR3 = findspecies(Mobj{i}, 'protein tetR3-lva:protein ClpX*');
iClpX_lacI3 = findspecies(Mobj{i}, 'protein lacI3-lva:protein ClpX*');
iClpX_tetR4 = findspecies(Mobj{i}, 'protein tetR4-lva:protein ClpX*');
iClpX_lacI4 = findspecies(Mobj{i}, 'protein lacI4-lva:protein ClpX*');
iClpX_unmature = findspecies(Mobj{i}, 'protein ClpX');
iClpX_mature = findspecies(Mobj{i}, 'protein ClpX*');
ilacI3_denatured = findspecies(Mobj{i}, 'protein lacI3-lva*');
itetR3_denatured = findspecies(Mobj{i}, 'protein tetR3-lva*');
ilacI4_denatured = findspecies(Mobj{i}, 'protein lacI4-lva*');
itetR4_denatured = findspecies(Mobj{i}, 'protein tetR4-lva*');

plot(t_ode{i}/3600, x_ode{i}(:,iClpX_tetR3),'r', t_ode{i}/3600,x_ode{i}(:,iClpX_lacI3),'b',...
    t_ode{i}/3600,x_ode{i}(:,iClpX_unmature),'k--', t_ode{i}/3600,x_ode{i}(:,iClpX_mature), 'k',...
    t_ode{i}/3600, x_ode{i}(:,itetR3_denatured),'m',t_ode{i}/3600, x_ode{i}(:,ilacI3_denatured),'c')
%{
,...
    t_ode{i}/3600, x_ode{i}(:,iClpX_tetR4),'r-.', t_ode{i}/3600,x_ode{i}(:,iClpX_lacI4),'b-.',...
    t_ode{i}/3600, x_ode{i}(:,itetR4_denatured),'m-.',t_ode{i}/3600, x_ode{i}(:,ilacI4_denatured),'c-.'
 'tetR4 denature','lacI4 denature',
'ClpX:tetR4','ClpX:lacI4',
    %}
    legend('ClpX:tetR3','ClpX:lacI3',...
    'ClpX unmature', 'ClpX mature', ...
    'tetR3 denature','lacI3 denature',...
    'Location', 'NorthEastOutside')
h = suptitle(['init tetR monomer = ' num2str(c(i,1)) ', init lacI monomer = ' num2str(c(i,2))]);
cd([mainDir '\examples\Vipul\Genetic Toggle\2013_manual_parameter_study\Sim' num2str(SimNumber)])
print('-dtiff','-r130',['Sim_' num2str(SimNumber) '_' num2str(i)])
saveas(gcf, ['Sim_' num2str(SimNumber) '_' num2str(i) '.fig'])
cd(mainDir)
% , 'k', t_ode{i}/3600, x_ode{i}(:,iNTP), 'k--',...
%     t_ode{i}/3600, x_ode{i}(:,iAA), 'k-.'
% , 'NTP', 'AA',

% 
% 
% %later
% subplot(4, 2, 6)
% iRNA_lacI1 = findspecies(Mobj{i}, 'RNA rbs--lacI3-lva');
% iRNA_lacI3 = findspecies(Mobj{i}, 'Ribo:RNA rbs--lacI3-lva');
% iRNA_tetR1 = findspecies(Mobj{i}, 'RNA rbs--tetR3-lva');
% iRNA_tetR3 = findspecies(Mobj{i}, 'Ribo:RNA rbs--tetR3-lva');
% plot(t_ode{i}/3600, x_ode{i}(:,iRNA_lacI1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNA_lacI3), 'k--',...
%     t_ode{i}/3600, x_ode{i}(:,iRNA_tetR1), 'r',t_ode{i}/3600, x_ode{i}(:,iRNA_tetR3), 'r--')
% legend('lacI RNA', 'ribo:lacI RNA', 'tetR RNA', 'ribo:tetR RNA', 'Location', 'NorthEastOutside')
% 

% % figure
% % plot(t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR3), 'k--')
% subplot(4, 2, 8)
% iaTc = findspecies(Mobj{i}, 'aTc');
% iaTcbound = findspecies(Mobj{i}, '2 aTc:protein tetR3-lvadimer');
% plot(t_ode{i}/3600, x_ode{i}(:,iaTc), 'k-',t_ode{i}/3600, x_ode{i}(:,iaTcbound), 'k--',...
%     t_ode{i}/3600, x_ode{i}(:,itetR), 'r', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'r--')
% legend('aTc', 'aTc:tetRdimer', 'tetR', 'tetRdimer', 'Location', 'NorthEastOutside')%
% dirstr = pwd;
% cd([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
% print('-dtiff','-r200',['TimeTraces__tetR_' num2str(c(i,1)) '__lacI_' num2str(c(i,2)) '___' folderdate])
% saveas(gcf, ['TimeTraces__tetR_' num2str(c(i,1)) '__lacI_' num2str(c(i,2)) '___' folderdate '.fig'])
% cd(dirstr)
close
end
% set(0,'DefaultFigureVisible','off');
% colororder1 = lines;
% colororder2 = [0.8 0 0;
%     0 0.8 0;
%     160/255 32/255 240/255;
%     0 1 1;
%     1 69/255 0;
%     112/255 138/255 144/255;
%     188/255 143/255 143/255
%     0 0 0.8;];
% colororder3 = [colororder2;colororder1];
% figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
% hold on
% for i = 1:numel(a)
%     itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
%     ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
%     h= plot(t_ode{i}/3600, x_ode{i}(:,itetR), '-',...
%         t_ode{i}/3600, x_ode{i}(:,ilacI),'--');
%     set(h(1), 'Color', colororder3(mod(i, size(colororder3,1)),:), 'LineWidth', 1.5)
%     set(h(2), 'Color', colororder3(mod(i, size(colororder3,1)),:), 'LineWidth', 1.5)
%    
%     
% end
% dirstr = pwd;
% cd([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
% print('-dtiff','-r200',['TimeTraces_tetR_lacI_' folderdate])
% saveas(gcf, ['TimeTraces_tetR_lacI_' folderdate '.fig'])
% cd(dirstr)
% close all
