function [Mobj, t_ode, x_ode, names] = gen_switch(varargin)


% gen_switch(varargin) - Genetic Toggle Switch example
% Vipul Singhal, September 2012
% Modified from negautoreg.m by R. M. Murray, 8 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a genetic switch.  
%
%    It can be run just by hitting F5, and uses default inputs in that
%    case. 
%
numvarargs = length(varargin);

%{
if numvarargs > 6
    error('myfuns:gen_switch:TooManyInputs', ...
        'requires at most 6 optional inputs');
end
%}
default1 = 0; %tetR_initialConc
default2 = 0; %lacI_initialConc
default3 = 5*60*60; %stopTime in seconds
default4 = 'both'; % activeInducer
default5 = 5; %aTc initial conc
default6 = 5; %IPTG initial conc


optargs = {default1, default2, default3, default4, default5, default6};
optargs(1:numvarargs) = varargin;
[tetR_initialConc, lacI_initialConc, stopTime, activeInducer, aTc_initialConc, ...
    IPTG_initialConc] = optargs{:};



% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('c1');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');
% Define the DNA strands (defines TX-TL species + reactions)
% check the ptrc2 and lac lengths. In Gardener et al (2000), plasmids are
% used for tetR and lac. 
dna_lacI = txtl_adddna(tube3,'thio-junk(500)-ptet(50)', 'rbs(20)', 'lacI(647)-lav(40)-terminator(100)', 5, 'linear');
dna_tetR = txtl_adddna(tube3, 'thio-junk(500)-ptrc2(50)', 'rbs(20)', 'tetR(647)-lav(40)-terminator(100)', 5, 'linear');
dna_deGFP = txtl_adddna(tube3, 'p70(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_adddna(tube3,  'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

if strcmp(activeInducer,'both')
    txtl_addspecies(tube2, 'aTc', aTc_initialConc);
    txtl_addspecies(tube2, 'IPTG', IPTG_initialConc);
    else if strcmp( activeInducer,'aTc')
        txtl_addspecies(tube2, 'aTc', aTc_initialConc); % express lacI, repress tetR
        else if strcmp( activeInducer,'IPTG')
            txtl_addspecies(tube2, 'IPTG', IPTG_initialConc);% express tetR, repress lacI
            end
        end
end

lacIprotein = sbioselect(tube3, 'Name','protein lacI');
tetRprotein = sbioselect(tube3, 'Name','protein tetR');
set(lacIprotein, 'InitialAmount', lacI_initialConc);
set(tetRprotein, 'InitialAmount', tetR_initialConc);


%debug code:
%disp('flag1')
%get(lacIprotein, 'InitialAmount')
%get(tetRprotein, 'InitialAmount')
%activeInducer
%pause(10)

%
% Next we have to set up the reactions that describe how the circuit
% works.  Transcription and translation are already included above, so
% we just need to include protein-protein and protein-DNA interactions.
%
% Note that the commands in this section are standard Simbiology commands,
% so you can put anything you want here.
%

% No additional reactions required for this circuit
% tetR-DNA interactions are automatically included in tetR setup

%
% Describe the actual experiment that we want to run.  This includes 
% combining the various tubes and also adding any additional inducers
% or purified proteins that you want to include in the run.
%

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);

txtl_setup_parameters(Mobj);
%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', stopTime)
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end
[t_ode, x_ode, names] = sbiosimulate(Mobj, configsetObj);

% Top row: protein and RNA levels
figure(1); clf(); subplot(2,1,1);
iLacI = findspecies(Mobj, 'protein lacI-lav-terminator')
iTetR = findspecies(Mobj, 'protein tetR-lav-terminator')
iGamS = findspecies(Mobj, 'protein gamS');
iGFP = findspecies(Mobj, 'protein deGFP');
iGFPs = findspecies(Mobj, 'protein deGFP*');

plot(t_ode/60, x_ode(:, iTetR),'k-', t_ode/60, x_ode(:, iLacI), 'b-', t_ode/60, x_ode(:, iGamS), 'r-', ...
  t_ode/60, x_ode(:, iGFP) + x_ode(:, iGFPs), 'g--', ...
  t_ode/60, x_ode(:, iGFPs), 'g-');

title('Gene Expression');
lgh = legend({'tetR', 'lacI', 'GamS', 'GFPt', 'GFP*'}, 'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

% Second row, left: resource limits
subplot(2,2,3);
iNTP = findspecies(Mobj, 'NTP');
iAA  = findspecies(Mobj, 'AA');
iRNAP  = findspecies(Mobj, 'RNAP70');
iRibo  = findspecies(Mobj, 'Ribo');
mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
plot(...
  t_ode/60, x_ode(:, iAA)/x_ode(1, iAA), 'b-', ...
  t_ode/60, x_ode(:, iNTP)/x_ode(1, iNTP), 'r-', ...
  t_ode/60, x_ode(:, iRNAP)/x_ode(1, iRNAP), 'b--', ...
  t_ode/60, x_ode(:, iRibo)/x_ode(1, iRibo), 'r--');

title('Resource usage');
lgh = legend(...
  {'NTP [mM]', 'AA [mM]', 'RNAP70 [nM]', 'Ribo [nM]'}, ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [normalized]');
xlabel('Time [min]');

% Second row, right: DNA and mRNA
subplot(2,2,4);
iDNA_lacI = findspecies(Mobj, 'DNA thio-junk-ptrc2--rbs--tetR-lav-terminator');
iDNA_tetR = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--lacI-lav-terminator');
iDNA_gamS = findspecies(Mobj, 'DNA p70--rbs--gamS');
iRNA_tetR = findspecies(Mobj, 'RNA rbs--tetR-lav-terminator');
iRNA_gamS = findspecies(Mobj, 'RNA rbs--gamS');
plot(t_ode/60, x_ode(:, iDNA_tetR), 'b-', ...
  t_ode/60, x_ode(:, iDNA_gamS), 'r-', ...
  t_ode/60, x_ode(:, iRNA_tetR), 'b--', ...
  t_ode/60, x_ode(:, iRNA_gamS), 'r--');

title('DNA and mRNA');
lgh = legend(...
  names([iDNA_tetR, iDNA_lacI, iDNA_gamS, iRNA_tetR, iRNA_gamS]), ...
  'Location', 'Northeast');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
clc
disp(['Inducer = ' activeInducer '; lacI init = ' num2str(lacI_initialConc)...
    '; tetR init = ' num2str(tetR_initialConc)])
pause(3)
end
%Run a parameter sweep with different initial concentrations of the protein
%species (lacI and tetR), and plot their evolution to show the Bistability
%of this system. 



%
% Run a set of experiments to explore the effect of inducers
%
%! TODO: write up this section

% Mobj = txtl_combine(tube1, 6, tube2, 2, tube3, 2);
% Add inducer at a given concentration (must be nanomolar)
% txtl_addspecies(Mobj, 'aTc', 50);

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
