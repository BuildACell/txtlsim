function [t_ode, x_ode, Mobj] = RNA_repression(varargin)

% RNA_repression.m 
% Vipul Singhal, April 2013
%
%
if nargin == 0
    ATT_DEGFP_conc = 1;
    ANTI_conc = 0.1;
    CONTROL_DNA_conc = 0;
elseif nargin == 3
    ATT_DEGFP_conc = varargin{1};
    ANTI_conc = varargin{2};
    CONTROL_DNA_conc = varargin{3};
end




% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E9_RNAcascade');
tube2 = txtl_buffer('E9_RNAcascade');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('RNA_repression');

% DNA with attenuator and deGFP
dna_att_deGFP = txtl_add_dna(tube3, ...
  'pJ23119(35)', 'att(287)-rbs(20)', 'deGFP(1080)', ATT_DEGFP_conc*4.2, 'linear');					

dna_gamS = txtl_add_dna(tube3, ...
  'p70(35)', 'rbs(20)', 'gamS(1000)', 0*4.2, 'plasmid'); 

% DNA with antisense RNA and dummy protein
dna_anti_dummyprotein = txtl_add_dna(tube3, ...
  'pJ23119(35)', 'anti(91)', 'no_protein', ...	% promoter, utr, gene
ANTI_conc*4.2, ...					% concentration (nM)
  'linear');					% type

dna_control_dummyprotein = txtl_add_dna(tube3, ...
  'pJ23119(35)', 'control(91)', 'no_protein', ...	% promoter, utr, gene
CONTROL_DNA_conc*4.2, ...					% concentration (nM)
  'linear');					% type

txtl_addspecies(tube3, 'protein gamS', 10, 'Internal');
% DNA with control RNA and dummy protein (not translated)

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulation

simulationTime = 14*60*60;


[simData] = txtl_runsim(Mobj,simulationTime);
t_ode = simData.Time;
x_ode = simData.Data;
end


% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
