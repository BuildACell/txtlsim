% ATP degradation work

close all; clear all; clc
tube1 = txtl_extract('E30');
tube2 = txtl_buffer('E30');
tube3 = txtl_newtube('gene_expression');
dna_deGFP = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 0, 'plasmid');

Mobj = txtl_combine([tube1, tube2, tube3]);




[simData] = txtl_runsim(Mobj,14*60*60);
t_ode = simData.Time;
x_ode = simData.Data;
ATP1 = findspecies(Mobj, 'ATP');
ATP2 = findspecies(Mobj, 'AA:ATP:Ribo:RNA rbs--deGFP:RNase');
ATP3 = findspecies(Mobj, 'AA:ATP:Ribo:RNA rbs--deGFP');

figure
subplot(2,1,1)
semilogy(t_ode/60, x_ode(:,ATP1)+x_ode(:,RNA2)+x_ode(:,RNA3)+x_ode(:,RNA4)+x_ode(:,RNA5)+x_ode(:,RNA6))
axis([0 240 0.01 300])
title('logy, RNA degradation, initial RNA conc = 200nM')

subplot(2,1,2)
plot(t_ode/60, x_ode(:,ATP1)+x_ode(:,RNA2)+x_ode(:,RNA3)+x_ode(:,RNA4)+x_ode(:,RNA5)+x_ode(:,RNA6))
axis([0 240 0.01 300])
title('RNA degradation, initial RNA conc = 200nM')