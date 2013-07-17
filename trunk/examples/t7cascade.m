% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E15');
tube2 = txtl_buffer('E15');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_t7rnap = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 't7RNAP(1000)', ...	% promoter, rbs, gene
   0.1, ...					% concentration (nM)
  'plasmid');					% type
dna_deCFP = txtl_add_dna(tube3, ...
  't7(50)', 'rbs(20)', 'deCFP(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% txtl_addreaction(Mobj,'protein t7RNAP -> null','MassAction',{'t7RNAP_deg',0.0011});
% txtl_addreaction(Mobj,'NTP:protein t7RNAP:DNA t7--rbs--deCFP -> NTP + DNA t7--rbs--deCFP','MassAction',{'t7RNAP_deg',0.0011});
% txtl_addreaction(Mobj,'protein t7RNAP:DNA t7--rbs--deCFP -> DNA t7--rbs--deCFP','MassAction',{'t7RNAP_deg',0.0011});


%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;

%%
figure(2)
hold on
plot(t_ode/60,sum(x_ode(:,[22 23 26]),2),'LineWidth',2)
plot(t_ode/60,sum(x_ode(:,[15 16 18]),2),'r','LineWidth',2)
xlabel('Time [min]')
ylabel('mRNA [nM]')
title('TXTL simulation')
legend('1 nM T7-mg-cfp mRNA',' 0.1nM pr-T7RNAP')


%txtl_plot(simData,Mobj);