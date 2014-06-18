close all; clear all

% 30 nM initial DNA conc, match figure 1c in PRL paper. 
tube1 = txtl_extract('E31VNPRL');
tube2 = txtl_buffer('E31VNPRL');
tube3 = txtl_newtube('gene_expression');
dna_deGFP = txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'deGFP(1000)', 30,  'plasmid');				
Mobj = txtl_combine([tube1, tube2, tube3]);
[simData] = txtl_runsim(Mobj,8*60*60);
txtl_plot(simData,Mobj);
t_ode = simData.Time;
    x_ode = simData.Data;
RNA1 = findspecies(Mobj, 'RNA rbs--deGFP');
RNA2 = findspecies(Mobj, 'Ribo:RNA rbs--deGFP');
RNA3 = findspecies(Mobj, 'AA:AGTP:Ribo:RNA rbs--deGFP') ;
RNA4 = findspecies(Mobj, 'RNA rbs--deGFP:RNase')       ;                   
RNA5 = findspecies(Mobj, 'AA:AGTP:Ribo:RNA rbs--deGFP:RNase')  ;                
RNA6 = findspecies(Mobj, 'Ribo:RNA rbs--deGFP:RNase')   ; 
figure
plot(t_ode/60, x_ode(:,RNA1)+x_ode(:,RNA2)+x_ode(:,RNA3)+x_ode(:,RNA4)+x_ode(:,RNA5)+x_ode(:,RNA6))
title('RNA levels, initial DNA conc = 30nM')

count = 0; Mobj = cell(6,1); simData = cell(6,1); t_ode = cell(6,1); x_ode = cell(6,1); 
for DNAinitial = [0.5 5 10 20 30 40]
    count = count+1;
    tube1 = txtl_extract('E31VNPRL');
    tube2 = txtl_buffer('E31VNPRL');
    tube3 = txtl_newtube('gene_expression');
    dna_deGFP = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', DNAinitial, 'plasmid');
    Mobj{count} = txtl_combine([tube1, tube2, tube3]);
    [simData{count}] = txtl_runsim(Mobj{count},8*60*60);
    t_ode{count} = simData{count}.Time;
    x_ode{count} = simData{count}.Data;
end
figure
colororder = lines;
for i = 1:6
    iGFP = findspecies(Mobj{i}, 'protein deGFP*')
    h(i) = plot(t_ode{i}/60, x_ode{i}(:,iGFP));
    hold on
    set(h(i), 'Color', colororder(i,:), 'LineWidth', 1.5);
    hold on
end        
           axis([0 55 0 1350])
        legend(h, {'0.5', '5', '10', '20', '30', '40'}, 'Location', 'NorthEastOutside')
    

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
