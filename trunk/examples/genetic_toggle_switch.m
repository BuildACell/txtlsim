% Genetic Toggle Switch example
% Vipul Singhal, September 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch. It shows the bistability of the circuit, and the
% switching that can be accomplished by the tetR repressing inducer aTc and
% the lacI repressing inducer IPTG. 
%


% Set up the standard TXTL tubes
tube1 = txtl_extract('E9');
tube2 = txtl_buffer('E9');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
dna_lacI = txtl_add_dna(tube3,'ptet(50)', 'rbs(20)', 'lacI(647)', 5, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'ptrc2(50)', 'rbs(20)', 'tetR(647)', 5, 'plasmid');
dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP(1000)', 5, 'plasmid');
dna_deCFP = txtl_add_dna(tube3, 'ptrc2(50)', 'rbs(20)', 'deCFP(1000)', 5, 'plasmid');


txtl_addspecies(tube3, 'aTc',0);
txtl_addspecies(tube3, 'IPTG',0);


% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
% sobj_aTc = sbioselect(Mobj, 'Type', 'species', 'Name', 'aTc');
% sobj_IPTG = sbioselect(Mobj, 'Type', 'species', 'Name', 'IPTG');
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'SolverType', 'ode15s');
set(configsetObj, 'StopTime', 14*60*60);

tic
simData = txtl_runsim(Mobj,configsetObj);
toc

t_ode = simData.Time;
x_ode = simData.Data;

%% Plot results

figure(1)
plot(t_ode/60,x_ode(:,[34 42]))
xlabel('Time [min]')
ylabel('GFP/CFP [nM]')
legend('GFP','CFP')






