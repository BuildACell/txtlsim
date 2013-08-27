% IFFL as tested in lab by S. Guo
% VS 2013

%% clean up

clear all
close all
clc


%% no clpX
% Set up the standard TXTL tubes


tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

tube3 = txtl_newtube('circuit_closed_loop_withClpX');

txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'AraC(600)',0.04*4.5, 'plasmid');
txtl_add_dna(tube3, 'pBAD(50)', 'rbs(20)', 'tetR(600)', 0.05*4.5, 'plasmid');
txtl_add_dna(tube3,'pBAD_ptet(150)', 'rbs(20)', 'deGFP(1000)-lva(20)',0.5*4.5, 'plasmid');
txtl_add_dna(tube3,'p70(50)', 'rbs(20)', 'ClpX(600)',0.1*4.5, 'plasmid');

txtl_addspecies(tube3, 'arabinose', 1000);
txtl_addspecies(tube3, 'aTc', 1000);
txtl_addspecies(tube3, 'protein ClpX', 100);
% set up well_a1
Mobj = txtl_combine([tube1, tube2, tube3]);

%% Run a simulation
 simulationTime = 12*60*60;

simData = txtl_runsim(Mobj,simulationTime);

t_ode = simData.Time;
x_ode = simData.Data;

figure
cellOfSpecies = {'RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:protein tetRdimer', 'RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:protein tetRdimer:arabinose:protein AraC'
                 'NTP:RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:protein tetRdimer','NTP:RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:protein tetRdimer:arabinose:protein AraC'
                 'RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:arabinose:protein AraC', 'RNAP70:DNA pBAD--rbs--tetR:arabinose:protein AraC'
                 'NTP:RNAP70:DNA pBAD_ptet--rbs--deGFP-lva:arabinose:protein AraC', 'NTP:RNAP70:DNA pBAD--rbs--tetR:arabinose:protein AraC'};
plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)
figure
cellOfSpecies = {'protein deGFP-lva*', 'protein tetR','protein ClpX*'
                 'protein AraC','protein tetRdimer','protein deGFP-lva*:protein ClpX*'
                 'arabinose', 'aTc','protein deGFP-lva***'
                 'arabinose:protein AraC', '2 aTc:protein tetRdimer','protein ClpX'};
plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)
figure
cellOfSpecies = {'RNAP', 'protein sigma70','Ribo'
                 'RNAP70','RNase','NTP'
                 'ATP', 'AA','protein deGFP-lva'};
plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies)

