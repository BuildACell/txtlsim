% RFP_loading.m - leading effects due to RFP example
% VIpul Singhal 4 March 2013
% In collaboration with Shaobin Guo
%
% 
close all
clear all

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');


RFPlevels = [0 1 2 3 5 10 20 40];%3 4 5 6 7 8 9 10 
t_ode = cell(size(RFPlevels));
x_ode = cell(size(RFPlevels));

count = 1;
for RFPconc = RFPlevels
% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'tetR-lva', 0.1*(10/2.25),  'linear');
dna_deGFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'deGFP-lva', 1*(10/2.25), 'linear');
dna_RFP = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'RFP(647)-lva', RFPconc*(10/2.25), 'linear');
%dna_gamS = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 0, 'plasmid');
gamS = txtl_addspecies(tube3, 'protein gamS', 100*(10/2.25));
txtl_addspecies(tube3, 'aTc', 200*(10/2.25));
% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 6*60*60)
[t_ode{count}, x_ode{count}, mObj, simData] = txtl_runsim(Mobj, configsetObj);
count = count + 1;
end
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--', 'g--'};

figure(1)
for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,1); hold on;
  idx = findspecies(Mobj, 'protein deGFP-lva*');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('protein deGFP-lva*')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,2); hold on;
  idx = findspecies(Mobj, 'protein RFP-lva*');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('protein RFP-lva*')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,3); hold on;
  idx = findspecies(Mobj, 'aTc');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('aTc')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');


for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,4); hold on;
  idx = findspecies(Mobj, 'AA');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('AA')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
%print('-dtiff','-r200',[pwd '\examples\Shaobin\high aTc case 1000nM\' 'figure1'])

% 
figure(2)
for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,1); hold on;
  idx = findspecies(Mobj, 'AGTP');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('AGTP')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,2); hold on;
  idx = findspecies(Mobj, 'protein tetR-lvadimer');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('protein tetR-lvadimer')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,3); hold on;
  idx = findspecies(Mobj, 'DNA ptet--rbs--RFP-lva:protein tetR-lvadimer');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('DNA ptet--rbs--RFP-lva:protein tetR-lvadimer')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
 subplot(2,2,4); hold on;
  idx = findspecies(Mobj, 'DNA ptet--rbs--deGFP-lva:protein tetR-lvadimer');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('DNA ptet--rbs--deGFP-lva:protein tetR-lvadimer')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
%print('-dtiff','-r200',[pwd '\examples\Shaobin\high aTc case 1000nM\' 'figure2'])

figure(3)
for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,1); hold on;
  idx = findspecies(Mobj, 'Ribo');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('Ribo')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,2); hold on;
  idx = findspecies(Mobj, 'protein gamS');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('protein gamS')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,2,3); hold on;
  idx = findspecies(Mobj, 'RNA rbs--deGFP-lva');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('RNA rbs--deGFP-lva')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
 subplot(2,2,4); hold on;
  idx = findspecies(Mobj, 'RNA rbs--RFP-lva');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('RNA rbs--RFP-lva')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
%print('-dtiff','-r200',[pwd '\examples\Shaobin\high aTc case 1000nM\' 'figure3'])

figure(4)
for count = 1:size(RFPlevels,2)
  % Plot the time trace  
  subplot(2,1,1); hold on;
  idx = findspecies(Mobj, 'Ribo:RNA rbs--deGFP-lva');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('Ribo:RNA rbs--deGFP-lva')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:size(RFPlevels,2)
  % Plot the time trace  
 subplot(2,1,2); hold on;
  idx = findspecies(Mobj, 'Ribo:RNA rbs--RFP-lva');
  plot(t_ode{count}, x_ode{count}(:, idx), colors{count});
  labels{count} = [num2str(RFPlevels(count)) ' nM RFP plasmid'];
  title('Ribo:RNA rbs--RFP-lva')
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
%print('-dtiff','-r200',[pwd '\examples\Shaobin\high aTc case 1000nM\' 'figure4'])

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
