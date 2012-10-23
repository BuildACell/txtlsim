% txtl_protein_sigma28.m - protein information for sigma28 factor
% Zoltan A Tuza,  Sep 2012
%
% This file contains a description of the protein produced by sigma28.
% Calling the function txtl_protein_sigma28() will set up the reactions for
% sequestration by the inducer aTc.

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

function [Rlist, genelen] = txtl_protein_sigma28(tube, protein, geneFull, genelen)

% set up gene default lengths
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
            case 'sigma28'
                genelen{geneDefIdx(i)} = 1000;
            case 'lva'
                genelen{geneDefIdx(i)} = 40; 
            case 'terminator'
                genelen{geneDefIdx(i)} = 100; 
        end
    end
end
%sequestration of RNAP by sigma28 factor
Kf = 100; % nM^-1s^-1
Kr = 0.1; % s^-1

%RNAP + Sigma70 <-> RNAP70
% Set up the reaction
Robj1 = addreaction(tube, ['RNAP + ' protein.Name ' <-> RNAP28']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', Kf);
Pobj1r = addparameter(Kobj1, 'kr', Kr);
set(Kobj1, 'ParameterVariableNames', {'kf','kr'});


% Return the list of reactions that we set up
Rlist = [Robj1];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
