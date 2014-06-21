% Accessing individual parameters: Mobj.UserData.ReactionConfig.RNase_F for example. 
% Accessing initial concentrations. Mobj.species(2).initialAmount (for the
% second species. can use findspecies to get tthe index). 

% We do a parameter study of the RNA degradation. 

% Nominal Parameters
close all; clear all; clc
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');
tube3 = txtl_newtube('gene_expression');
dna_deGFP = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 0, 'plasmid');

Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj, 'RNA rbs--deGFP', 200);



[simData] = txtl_runsim(Mobj,14*60*60);
t_ode = simData.Time;
x_ode = simData.Data;
RNA1 = findspecies(Mobj, 'RNA rbs--deGFP');
RNA2 = findspecies(Mobj, 'Ribo:RNA rbs--deGFP');
RNA3 = findspecies(Mobj, 'AA:AGTP:Ribo:RNA rbs--deGFP') ;
RNA4 = findspecies(Mobj, 'RNA rbs--deGFP:RNase')       ;                   
RNA5 = findspecies(Mobj, 'AA:AGTP:Ribo:RNA rbs--deGFP:RNase')  ;                
RNA6 = findspecies(Mobj, 'Ribo:RNA rbs--deGFP:RNase')   ;   

figure
subplot(2,1,1)
semilogy(t_ode/60, x_ode(:,RNA1)+x_ode(:,RNA2)+x_ode(:,RNA3)+x_ode(:,RNA4)+x_ode(:,RNA5)+x_ode(:,RNA6))
axis([0.01 240 0.01 300])
title('logy, RNA degradation, initial RNA conc = 200nM')

subplot(2,1,2)
plot(t_ode/60, x_ode(:,RNA1)+x_ode(:,RNA2)+x_ode(:,RNA3)+x_ode(:,RNA4)+x_ode(:,RNA5)+x_ode(:,RNA6))
axis([0 240 0.01 300])
title('RNA degradation, initial RNA conc = 200nM')

% The graph has an initial half life of 15 min, about what we want, and this half life stays at this value. 
% This is pretty good. What we would like is to get the number closer to
% 12 min. 
% Lets see what happens when we change initial RNase conc, the k_deg, and the 
% k_on and k_off. 

% Actually lets do this later. For now, this is good enough. We can try to get it close to 12 nM in the next iteration. 

