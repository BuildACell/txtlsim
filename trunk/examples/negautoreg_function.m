
function [t_ode, x_ode, Mobj] =  negautoreg_function(extract,dna_amount,tspan,atc_amount)

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract(extract);
tube2 = txtl_buffer(extract);

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');


% Define the DNA strands (defines TX-TL species + reactions)
%dna_tetR = txtl_add_dna(tube3, 'thio-junk(500)-ptet(50)', 'rbs(20)', 'tetR(1200)-lva(40)-terminator(100)', dna_amount*4.2, 'plasmid');
dna_tetR = txtl_add_dna(tube3, 'ptet(50)', 'rbs(20)', 'tetR(1200)', dna_amount*4.2, 'plasmid');



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
Mobj = txtl_combine([tube1, tube2, tube3]);
if ~isempty(atc_amount)
 txtl_addspecies(Mobj, 'aTc', atc_amount);
end
 

% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', tspan*60*60)
set(configsetObj, 'SolverType', 'ode15s');
[t_ode, x_ode] = txtl_runsim(Mobj, configsetObj);


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
