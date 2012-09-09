% txtl_prom_p70.m - promoter information for p70 promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the standard p70 promoter.
% Calling the function txtl_prom_p70() will set up the reactions for
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

function Rlist = txtl_prom_p70(tube, dna, rna)

% Parameters that describe this promoter
%! TODO: replace these values with correct values
kf70 = log(2)/0.1;			% 100 ms bind rate
kr70 = 10/kf70;				% Km of 10 nM (from VN model)
ktx70 = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
DNAbound = ['[RNAP70:' dna.Name ']'];

% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> ' DNAbound]);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf70);
Pobj1r = addparameter(Kobj1, 'kr', kr70);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

% Set up the transcription reaction
Robj2 = addreaction(tube, [DNAbound ' + NTP -> ' DNA ' + ' RNA ' + ' RNAP]);
Kobj2 = addkineticlaw(Robj2, 'MassAction');
Pobj2 = addparameter(Kobj2, 'ktx', ktx70);
set(Kobj2, 'ParameterVariableNames', {'ktx'});

Rlist = [Robj1, Robj2];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
