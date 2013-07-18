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
  't7_lacI(50)', 'rbs(20)', 'lacI(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type

txtl_add_dna(tube3, ...
  't7_lacI(50)', 'rbs(20)', 'alsS(1715)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type
txtl_add_dna(tube3, ...
  't7_lacI(50)', 'rbs(20)', 'budC(768)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type
txtl_add_dna(tube3, ...
  't7_lacI(50)', 'rbs(20)', 'alsD(768)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type





% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);


%% pathway reactions

txtl_addspecies(Mobj,'pyruvate',60000)
txtl_addspecies(Mobj,'IPTG',500);
txtl_addspecies(Mobj,'NADH',1000);



% 1st Stage
txtl_addreaction(Mobj,'protein alsS + pyruvate <-> protein alsS:pyruvate','MassAction',{'F',0.1;'R',0.01});
txtl_addreaction(Mobj,'protein alsS:pyruvate -> protein alsS + ACTL ','MassAction',{'F',0.1});

% 2nd Stage
txtl_addreaction(Mobj,'protein alsD + ACTL <-> protein alsD:ACTL','MassAction',{'F',0.1;'R',0.01});
txtl_addreaction(Mobj,'protein alsD:ACTL -> protein alsD + ACT ','MassAction',{'F',0.1});

% 3rd Stage
txtl_addreaction(Mobj,'protein budC + ACT + NADH <-> protein budC:ACT:NADH','MassAction',{'F',0.1;'R',0.01});
txtl_addreaction(Mobj,'protein budC:ACT:NADH -> protein budC + BDO + NADPlus','MassAction',{'F',0.1});




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
figure(1)
hold on
plot(t_ode/60,sum(x_ode(:,[50 52]),2),'b')
plot(t_ode/60,sum(x_ode(:,[53 54]),2),'r')
plot(t_ode/60,sum(x_ode(:,[55 56]),2),'g')
plot(t_ode/60,x_ode(:,57),'c')
hold off
legend('pyruvate','ACTL','ACT','BDO')

