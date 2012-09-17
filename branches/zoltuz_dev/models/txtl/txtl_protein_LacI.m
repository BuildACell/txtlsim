% txtl_protein_LacI.m - protein information for LacI
% RMM, 9 Sep 2012
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

function Rlist = txtl_protein_LacI(tube, protein)

% Parameters that describe this RBS
kf_aTc = 1; kr_aTc = 0.1; 

% Set up the binding reaction
%Robj1 = addreaction(tube, [protein.Name ' + aTc <-> aTc:' protein.Name]);
%Kobj1 = addkineticlaw(Robj1, 'MassAction');
%Pobj1f = addparameter(Kobj1, 'kf', kf_aTc);
%Pobj1r = addparameter(Kobj1, 'kr', kr_aTc);
%set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

% Set up degradation
Robj2 = addreaction(tube, [protein.Name ' -> null']);
Kobj2 = addkineticlaw(Robj2,'MassAction');
Pobj2 = addparameter(Kobj2,  'kf', 0.001);
set(Kobj2, 'ParameterVariableNames','kf');

Rlist = [Robj2];
%Rlist = [];


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
