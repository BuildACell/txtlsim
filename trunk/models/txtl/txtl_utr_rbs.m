% txtl_utr_rbs.m - promoter information for standard RBS
% RMM, 8 Sep 2012
%
% This file contains a description of a standard ribosome binding site.
% Calling the function txtl_utr_rbs() will set up the reactions for
% translation with the measured binding rates and translation rates.

% Written by Richard Murray, Sep 2012
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

function [Sobj, rbslen] = txtl_utr_rbs(tube, rna, protein, rbsFull, rbslen)


% set up promoter default lengths
rbsDefaultUsed = 0;
for i = 1: length(rbsFull)
    if isempty(rbslen{i})
        rbsDefaultUsed = rbsDefaultUsed+1;
        rbsDefIdx(rbsDefaultUsed) = i; %idx of segments to set defaults for
    end
end

if rbsDefaultUsed ~= 0
    for i = 1:length(rbsDefIdx)
        switch rbsFull{rbsDefIdx(i)}
            case 'rbs'
                rbslen{rbsDefIdx(i)} = 20;
            case 'spacer'
                rbslen{rbsDefIdx(i)} = 200; 
        end
    end
end

% Parameters that describe this RBS
%! TODO: replace these values with correct values
kf_rbs = log(2)/0.1;			% 100 ms bind rate
kr_rbs = 0.05 * kf_rbs;			% Km of ~0.05 (from VN model)

% Create strings for the reactants and products
RNA = ['[' rna.Name ']'];
Protname = protein.Name;

% Set up species for bound reaction
Sobj = addspecies(tube, ['Ribo:' rna.Name]);

% Set up the binding reaction
Robj = addreaction(tube, [RNA ' + Ribo <-> [' Sobj.Name ']']);
Kobj = addkineticlaw(Robj, 'MassAction');
Pobjf = addparameter(Kobj, 'kf', kf_rbs);
Pobjr = addparameter(Kobj, 'kr', kr_rbs);
set(Kobj, 'ParameterVariableNames', {'kf', 'kr'});

% Return the list of reactions that we set up
%! TODO: optional return list of all reactions
% Rlist = Robj;

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
