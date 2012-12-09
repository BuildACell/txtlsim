% txtl_protein_lacI.m - protein information for lacI
% RMM, 9 Sep 2012
%! TODO: header information incorrect (RMM didn't write this)
%! TODO: inconsistent naming - genes should start with lower case (lacI)
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

% Written by Richard Murray, 9 Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in the 
%      documentation and/or other materials provided with the distribution.
%
%   3. The name of the author may not be used to endorse or promote products 
%      derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function varargout = txtl_protein_lacI(mode, tube, protein, varargin)

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')

    geneData = [varargin{1};varargin{2}];
    defaultBasePairs = {'lacI','lva','terminator';647,40,100};
    geneData = txtl_setup_default_basepair_length(tube,geneData,...
        defaultBasePairs);
    
    varargout{1} = geneData(2,:);
    
    coreSpecies = {'IPTG',['IPTG:' protein.Name]};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)));
    
    
    % call additional functions to setup any other relevant species (like
    % multimers)
    txtl_protein_dimerization(mode, tube, protein); 
    txtl_protein_tetramerization(mode, tube, protein);
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode, 'Setup Reactions')
    
    % Parameters that describe this reaction
     kf_IPTG = 0.1; kr_IPTG = 0.01; 

    % Set up the binding reaction
    Robj1 = addreaction(tube, [protein.Name ' + IPTG <-> IPTG:' protein.Name]);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    Pobj1f = addparameter(Kobj1, 'TXTL_INDUCER_LACI_IPTG_F', kf_IPTG);
    Pobj1r = addparameter(Kobj1, 'TXTL_INDUCER_LACI_IPTG_R', kr_IPTG);
    set(Kobj1, 'ParameterVariableNames', {'TXTL_INDUCER_LACI_IPTG_F', 'TXTL_INDUCER_LACI_IPTG_R'});

    % degrade the IPTG inducer
    kf_IPTGdeg = 0.0001;
    Robj2 = addreaction(tube, 'IPTG -> null');
    Kobj2 = addkineticlaw(Robj2, 'MassAction');
    Pobj2 = addparameter(Kobj2, 'TXTL_INDUCER_DEGRADATION_IPTG', kf_IPTGdeg);
    set(Kobj2, 'ParameterVariableNames', {'TXTL_INDUCER_DEGRADATION_IPTG'});

    %Set up dimerization
    % Hsieh & Brenowitz 1997 JBC
    kf_dimer = 0.08637; % 1/(molecule*sec)
    kr_dimer =0.1; %0.00000001; % 1/sec

    txtl_protein_dimerization(mode, tube,protein, [kf_dimer, kr_dimer]); 

    %Set up tetramerization
    % Hsieh & Brenowitz 1997 JBC
    kf_tetramer = 0.00000602; % 1/(molecule*sec)
    kr_tetramer = 0.000001; %0.000001; % 1/sec
    txtl_protein_tetramerization(mode, tube,protein,[kf_tetramer,kr_tetramer]);

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_protein_lacI:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    





% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
