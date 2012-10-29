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

function [Rlist, genelen] = txtl_protein_lacI(tube, protein, geneFull, genelen)
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
            case 'lacI'
                genelen{geneDefIdx(i)} = 647;
            case 'lva'
                genelen{geneDefIdx(i)} = 40; 
            case 'terminator'
                genelen{geneDefIdx(i)} = 100; 
        end
    end
end
% Parameters that describe this RBS
 kf_IPTG = 0.1; kr_IPTG = 0.01; 

% Set up the binding reaction
Robj1 = addreaction(tube, [protein.Name ' + IPTG <-> IPTG:' protein.Name]);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_IPTG);
Pobj1r = addparameter(Kobj1, 'kr', kr_IPTG);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

% degrade the IPTG inducer
kf_IPTGdeg = 0.0001;
Robj2 = addreaction(tube, ['IPTG -> null']);
Kobj2 = addkineticlaw(Robj2, 'MassAction');
Pobj2 = addparameter(Kobj2, 'kf', kf_IPTGdeg);
set(Kobj2, 'ParameterVariableNames', {'kf'});

Rlist = [Robj1, Robj2];

%Set up dimerization
% Hsieh & Brenowitz 1997 JBC
kf_dimer = 0.08637; % 1/(molecule*sec)
kr_dimer =0.1; %0.00000001; % 1/sec

Rlist(end+1) = txtl_protein_dimerization(tube,protein,[kf_dimer,kr_dimer]); 



%Set up tetramerization
% Hsieh & Brenowitz 1997 JBC
kf_tetramer = 0.00000602; % 1/(molecule*sec)
kr_tetramer = 0.000001; %0.000001; % 1/sec
Rlist(end+1) = txtl_protein_tetramerization(tube,protein,[kf_tetramer,kr_tetramer]);








% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
