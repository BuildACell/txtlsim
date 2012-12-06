% txtl_protein_lambda.m - protein information for lambda
% VS Dec 2012
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

% Written by Vipul Singhal 2012
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

function varargout = txtl_protein_lambda(mode, tube, protein, varargin)
% in 'setup Species' mode, it returns an array of gene lengths, having
% added defaults in places where the lengths are missing. 

if strcmp(mode, 'Setup Species')

    geneFull = varargin{1};
    genelen = varargin{2};
    % set up gene default lengths of genes
    geneDefaultUsed = 0;
    for i = 1: length(geneFull)
        if isempty(genelen{i})        
            geneDefaultUsed = geneDefaultUsed+1;
            geneDefIdx(geneDefaultUsed) = i; %idx of segments to set defaults for
        end
    end
    if geneDefaultUsed ~= 0
        for i = 1:length(geneDefIdx)
            switch geneFull{geneDefIdx(i)}
                case 'lambda'
                    genelen{geneDefIdx(i)} = 647;
                case 'lva'
                    genelen{geneDefIdx(i)} = 40; 
                case 'terminator'
                    genelen{geneDefIdx(i)} = 100; 
            end
        end
    end

    % setup the species
%     foo = sbioselect(tube, 'Name', 'aTc');
%     if isempty(foo)
%         addspecies(tube, 'aTc');
%     end
%     foo = [];
%     foo = sbioselect(tube, 'Name', ['aTc:' protein.Name]);
%     if isempty(foo)
%         addspecies(tube, ['aTc:' protein.Name]);
%     end
%     foo = [];
    

    % call other functions in 'Setup Species' mode
    txtl_protein_dimerization('Setup Species', tube,protein);

    %return the gene lengths if not specified by the user
    varargout{1} = genelen;

elseif strcmp(mode, 'Setup Reactions')
    
    %!TODO: set up user defined parameters. 
%     kf_aTc = 1; kr_aTc = 0.1; 
%     % Set up the binding reaction
%     Robj1 = addreaction(tube, ['[' protein.Name '] + aTc <-> [aTc:' protein.Name ']']);
%     Kobj1 = addkineticlaw(Robj1, 'MassAction');
%     Pobj1f = addparameter(Kobj1, 'TXTL_INDUCER_TETR_ATC_F', kf_aTc);
%     Pobj1r = addparameter(Kobj1, 'TXTL_INDUCER_TETR_ATC_R', kr_aTc);
%     set(Kobj1, 'ParameterVariableNames', {'TXTL_INDUCER_TETR_ATC_R', 'TXTL_INDUCER_TETR_ATC_R'});

    % degrade the aTc inducer
%     kf_aTcdeg = 0.0001;
%     Robj2 = addreaction(tube, ['aTc -> null']);
%     Kobj2 = addkineticlaw(Robj2, 'MassAction');
%     Pobj2 = addparameter(Kobj2, 'TXTL_INDUCER_DEGRADATION_ATC', kf_aTcdeg);
%     set(Kobj2, 'ParameterVariableNames', {'TXTL_INDUCER_DEGRADATION_ATC'});


    txtl_protein_dimerization('Setup Reactions', tube,protein, [0.001, 0.0001]);

    % % ! TODO: Check if tetR undergoes tertramerization
    % % Hsieh & Brenowitz 1997 JBC
    % %! TODO: these may be strain/excract dependent
    % kf_tetramer = 0.000602; % 1/(molecule*sec)
    % kr_tetramer = 0.000001; % 1/sec
    % Rlist(end+1) = txtl_protein_tetramerization(tube,protein,[kf_tetramer,kr_tetramer]);
else
    error('txtltoolbox:txtl_protein_lambda:undefinedmode', 'The possible modes are ''Setup Species'' and ''Setup Reactions''.')
end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
