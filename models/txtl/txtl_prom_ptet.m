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

function [Rlist, promlen] = txtl_prom_ptet(tube, dna, rna, promFull, promlen)

% set up promoter default lengths
promDefaultUsed = 0;
for i = 1: length(promFull)
    if isempty(promlen{i})
        promDefaultUsed = promDefaultUsed+1;
        promDefIdx(promDefaultUsed) = i; %idx of segments to set defaults for
    end
end

if promDefaultUsed ~= 0
    for i = 1:length(promDefIdx)
        switch promFull{promDefIdx(i)}
            case 'ptet'
                promlen{promDefIdx(i)} = 50;
            case 'junk'
                promlen{promDefIdx(i)} = 500; 
            case 'thio'
                promlen{promDefIdx(i)} = 0; 
        end
    end
end
% Parameters that describe this promoter
%! TODO: replace these values with correct values
kf_ptet = log(2)/0.1;			% 100 ms bind rate
kr_ptet = 10 * kf_ptet;			% Km of 10 (same as p70, from VN)
ktx_ptet = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
RNAPbound = ['RNAP70:' dna.Name];

% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_ptet);
Pobj1r = addparameter(Kobj1, 'kr', kr_ptet);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

%
% Now put in the reactions for the utilization of NTPs
% Use an enzymatic reaction to proper rate limiting
%

Rlist1 = txtl_rnap_rnap70(tube, dna, rna, RNAPbound);

%
% Add reactions for sequestration of promoter by tetRdimer 
%

%! TODO: RMM, 29 Sep 2012
%! TODO: txtl_protein_tetR defines tetramers, which aren't used
%! TODO: proper implementation for tetR is via two operator sites (I think)
% VS: yes, there are 2 operators, see for example, 
% C Berens and W. Hillen, Gene regulation by tetracyclines, Eur. J. Biochem. 270, 3109�3121 (2003)

kf1_tetR = 0.2; kr1_tetR = 1;		% reaction rates (from sbio)
%! TODO: the 'DNA tetR needs to be generalized
Robj4 = addreaction(tube, ...
  [DNA ' + [protein tetRdimer] <-> [' dna.name ':protein tetRdimer]']);
Kobj4 = addkineticlaw(Robj4,'MassAction');
Pobj4 = addparameter(Kobj4, 'k4', kf1_tetR);
Pobj4r = addparameter(Kobj4, 'k4r', kr1_tetR);
set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

kf2_tetR = 0.2; kr2_tetR = 1;		% reaction rates (from sbio)
%! TODO: the 'DNA tetR needs to be generalized
Robj5 = addreaction(tube, ...
  ['[' dna.name ':protein tetRdimer] <-> [' dna.name ':protein tetRdimer:protein tetRdimer]']);
Kobj5 = addkineticlaw(Robj5,'MassAction');
Pobj5 = addparameter(Kobj5, 'k5', kf2_tetR);
Pobj5r = addparameter(Kobj5, 'k5r', kr2_tetR);
set(Kobj5, 'ParameterVariableNames', {'k5', 'k5r'});

Rlist = [Robj1, Rlist1, Robj4, Robj5];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End: