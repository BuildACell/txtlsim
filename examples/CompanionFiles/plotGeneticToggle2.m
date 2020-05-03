
close all
set(0,'DefaultFigureVisible','off');
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
hold on
tetR_SS = [24 0.02];
lacI_SS = [0.2 1.1];
for i = 1:numel(a)
    itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
    ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
    
    if norm([x_ode{i}(end,itetR),x_ode{i}(end,ilacI)]-tetR_SS) < norm([x_ode{i}(end,itetR),x_ode{i}(end,ilacI)]-lacI_SS)
    plot(x_ode{i}(:,itetR), x_ode{i}(:,ilacI), 'r', x_ode{i}(end,itetR), x_ode{i}(end,ilacI), 'r*',x_ode{i}(1,itetR), x_ode{i}(1,ilacI), 'ro')
    else
     plot(x_ode{i}(:,itetR), x_ode{i}(:,ilacI), 'b', x_ode{i}(end,itetR), x_ode{i}(end,ilacI), 'b*',x_ode{i}(1,itetR), x_ode{i}(1,ilacI), 'bo')   
    end
end
title('tetR (x axis) and lacI (y axis)')
folderdate = datestr(now,'yyyymmmmdd_HHMMSS');
mkdir([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
dirstr = pwd;

cd([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
print('-dtiff','-r200',['allTrajectories_' folderdate])
saveas(gcf, ['allTrajectories_' folderdate '.fig'])
cd(dirstr)
close 
set(0,'DefaultFigureVisible','off');
for i = 1:numel(a)
    
figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
% plot protein conc
subplot(4, 2, 1)
ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
ilacIdimer = findspecies(Mobj{i}, 'protein lacI3-lvadimer');
ilacItetramer = findspecies(Mobj{i}, 'protein lacI3-lvatetramer');
plot(t_ode{i}/3600, x_ode{i}(:,ilacI), 'k', t_ode{i}/3600, x_ode{i}(:,ilacIdimer), 'k--', t_ode{i}/3600, x_ode{i}(:,ilacItetramer), 'k.-');
legend('lacI', 'lacIdimer', 'lacItetramer', 'Location', 'NorthEastOutside')

subplot(4, 2, 2)
itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
itetRdimer = findspecies(Mobj{i}, 'protein tetR3-lvadimer');
plot(t_ode{i}/3600, x_ode{i}(:,itetR), 'k', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'k--')
legend('tetR', 'tetRdimer', 'Location', 'NorthEastOutside')

%plot resources
subplot(4, 2, 3)
iATP = findspecies(Mobj{i}, 'ATP');
iNTP = findspecies(Mobj{i}, 'NTP');
iAA = findspecies(Mobj{i}, 'AA');
iClpX = findspecies(Mobj{i}, 'protein ClpX*');
plot(t_ode{i}/3600, x_ode{i}(:,iATP), 'k', t_ode{i}/3600, x_ode{i}(:,iNTP), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iAA), 'k.-', t_ode{i}/3600, x_ode{i}(:,iClpX), 'r')
legend('ATP', 'NTP', 'AA', 'ClpX', 'Location', 'NorthEastOutside')
subplot(4, 2, 4)

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


subplot(4, 2, 5)
iClpX_complex1 = findspecies(Mobj{i}, 'protein tetR3-lva:protein ClpX*');
iClpX_complex2 = findspecies(Mobj{i}, 'protein lacI3-lva:protein ClpX*');
iClpX_unmature = findspecies(Mobj{i}, 'protein ClpX');



plot(t_ode{i}/3600, x_ode{i}(:,iClpX_complex1)+x_ode{i}(:,iClpX_complex2)+x_ode{i}(:,iClpX_unmature)+x_ode{i}(:,iClpX), 'k',...
    t_ode{i}/3600, x_ode{i}(:,iClpX_complex1),'b',t_ode{i}/3600, x_ode{i}(:,iClpX_complex2),'r',...
    t_ode{i}/3600, x_ode{i}(:,iClpX),'g',t_ode{i}/3600, x_ode{i}(:,iClpX_unmature),'c')
legend('all ClpX', 'ClpX tetR', 'ClpX lacI', 'ClpX*', 'ClpX', 'Location', 'NorthEastOutside')
%title('ClpX');
subplot(4, 2, 6)
iRNA_lacI1 = findspecies(Mobj{i}, 'RNA rbs--lacI3-lva');
iRNA_lacI3 = findspecies(Mobj{i}, 'Ribo:RNA rbs--lacI3-lva');
iRNA_tetR1 = findspecies(Mobj{i}, 'RNA rbs--tetR3-lva');
iRNA_tetR3 = findspecies(Mobj{i}, 'Ribo:RNA rbs--tetR3-lva');
plot(t_ode{i}/3600, x_ode{i}(:,iRNA_lacI1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNA_lacI3), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,iRNA_tetR1), 'r',t_ode{i}/3600, x_ode{i}(:,iRNA_tetR3), 'r--')
legend('lacI RNA', 'ribo:lacI RNA', 'tetR RNA', 'ribo:tetR RNA', 'Location', 'NorthEastOutside')

subplot(4, 2, 7)
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
% figure
% plot(t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR1), 'k-',t_ode{i}/3600, x_ode{i}(:,iRNAP_tetR3), 'k--')
subplot(4, 2, 8)
iaTc = findspecies(Mobj{i}, 'aTc');
iaTcbound = findspecies(Mobj{i}, '2 aTc:protein tetR3-lvadimer');
plot(t_ode{i}/3600, x_ode{i}(:,iaTc), 'k-',t_ode{i}/3600, x_ode{i}(:,iaTcbound), 'k--',...
    t_ode{i}/3600, x_ode{i}(:,itetR), 'r', t_ode{i}/3600, x_ode{i}(:,itetRdimer), 'r--')
legend('aTc', 'aTc:tetRdimer', 'tetR', 'tetRdimer', 'Location', 'NorthEastOutside')%
h = suptitle(['init tetR monomer = ' num2str(c(i,1)) ', init lacI monomer = ' num2str(c(i,2))]);
dirstr = pwd;
cd([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
print('-dtiff','-r200',['TimeTraces__tetR_' num2str(c(i,1)) '__lacI_' num2str(c(i,2)) '___' folderdate])
saveas(gcf, ['TimeTraces__tetR_' num2str(c(i,1)) '__lacI_' num2str(c(i,2)) '___' folderdate '.fig'])
cd(dirstr)
close
end
set(0,'DefaultFigureVisible','off');
colororder1 = lines;
colororder2 = [0.8 0 0;
    0 0.8 0;
    160/255 32/255 240/255;
    0 1 1;
    1 69/255 0;
    112/255 138/255 144/255;
    188/255 143/255 143/255
    0 0 0.8;];
colororder3 = [colororder2;colororder1];
figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.3])%,'Name','Simulation Plot Window','NumberTitle','off')
hold on
for i = 1:numel(a)
    itetR = findspecies(Mobj{i}, 'protein tetR3-lva');
    ilacI = findspecies(Mobj{i}, 'protein lacI3-lva');
    h= plot(t_ode{i}/3600, x_ode{i}(:,itetR), '-',...
        t_ode{i}/3600, x_ode{i}(:,ilacI),'--');
    set(h(1), 'Color', colororder3(mod(i, size(colororder3,1)),:), 'LineWidth', 1.5)
    set(h(2), 'Color', colororder3(mod(i, size(colororder3,1)),:), 'LineWidth', 1.5)
   
    
end
dirstr = pwd;
cd([pwd '\examples\Vipul\Genetic Toggle\' folderdate])
print('-dtiff','-r200',['TimeTraces_tetR_lacI_' folderdate])
saveas(gcf, ['TimeTraces_tetR_lacI_' folderdate '.fig'])
cd(dirstr)
close all
