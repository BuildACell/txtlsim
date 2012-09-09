% txtl_prom_ptet.m - promoter information for ptet promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.

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

function Rlist = txtl_prom_ptet(tube, dna, rna)

% Parameters that describe this promoter
%! TODO: replace these values with correct values
kf_ptet = log(2)/0.1;			% 100 ms bind rate
kr_ptet = 10/kf_ptet;			% Km of 10 (same as p70, from VN)
ktx_ptet = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
DNAbound = ['[RNAP70:' dna.Name ']'];

% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> ' DNAbound]);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_ptet);
Pobj1r = addparameter(Kobj1, 'kr', kr_ptet);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

% Set up the transcription reaction
Robj2 = addreaction(tube, [DNAbound ' + NTP -> ' DNA ' + ' RNA ' + ' RNAP]);
Kobj2 = addkineticlaw(Robj2, 'MassAction');
Pobj2 = addparameter(Kobj2, 'ktx', ktx_ptet);
set(Kobj2, 'ParameterVariableNames', {'ktx'});

% Add reactions for sequestration of promoter by TetR 
kf_tetR = 0.2; kr_tetR = 1;		% reaction rates (from sbio)
Robj3 = addreaction(tube, ...
  [DNA ' + [protein tetR] <-> [DNA tetR:protein tetR]']);
Kobj3 = addkineticlaw(Robj3,'MassAction');
Pobj3 = addparameter(Kobj3, 'k3', kf_tetR);
Pobj3r = addparameter(Kobj3, 'k3r', kr_tetR);
set(Kobj3, 'ParameterVariableNames', {'k3', 'k3r'});

Rlist = [Robj1, Robj2, Robj3];

% Local variables:
% mode: matlab
% End:
