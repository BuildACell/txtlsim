% txtl_prom_p28_ptet.m - promoter information for p28 and ptet combinatorial promoter
% Zoltan Tuza, Oct 2012
%
% This file contains a description of the p28 and ptet combinatorial promoter.
% Calling the function txtl_prom_p28_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% 
% 

% Written by Zoltan Tuza, Oct 2012
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

function [Rlist, promlen] = txtl_prom_p28_ptet(tube, dna, rna, promFull, promlen)


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
            case 'p28_ptet'
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
RNAP = 'RNAP28';			% RNA polymerase name for reactions
RNAPbound = ['RNAP28:' dna.Name];
P1 = 'protein sigma28';
P2 = 'protein tetRdimer';

% %Set up binding reaction for sigma28
% Robj0 = addreaction(tube, [DNA ' + ' P1 ' <-> ' dna.Name ':' P1 ]);
% Kobj0 = addkineticlaw(Robj0, 'MassAction');
% Pobj0f = addparameter(Kobj0, 'kf', kf_ptet);
% Pobj0r = addparameter(Kobj0, 'kr', kr_ptet);
% set(Kobj0, 'ParameterVariableNames', {'kf', 'kr'});
% 
% % RNAP complex
% Robj1 = addreaction(tube, [dna.Name ':' P1 ' + ' RNAP ' <-> [' RNAPbound ':' P1 ']']);
% Kobj1 = addkineticlaw(Robj1, 'MassAction');
% Pobj1f = addparameter(Kobj1, 'kf', kf_ptet);
% Pobj1r = addparameter(Kobj1, 'kr', kr_ptet);
% set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

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

% TX
Rlist1 = txtl_rnap_rnap28(tube, dna, rna, RNAPbound);


% 
% Set up bindig reaction for both sigma28 and tetR
Robj3 = addreaction(tube, [dna.Name ' + ' P2 ' <-> ' dna.Name ':' P2 ]);
Kobj3 = addkineticlaw(Robj3, 'MassAction');
Pobj3f = addparameter(Kobj3, 'kf', kf_ptet/2);
Pobj3r = addparameter(Kobj3, 'kr', kr_ptet/2);
set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
% 
% Robj4 = addreaction(tube, [dna.Name ':' P1 ' + ' P2 ' <-> ' dna.Name ':' P1 ':' P2 ]);
% Kobj4 = addkineticlaw(Robj4, 'MassAction');
% Pobj4f = addparameter(Kobj4, 'kf', kf_ptet*2);
% Pobj4r = addparameter(Kobj4, 'kr', kr_ptet);
% set(Kobj4, 'ParameterVariableNames', {'kf', 'kr'});
% 
% Robj5 = addreaction(tube, [dna.Name ':' P2 ' + ' P1 ' <-> ' dna.Name ':' P1 ':' P2 ]);
% Kobj5 = addkineticlaw(Robj5, 'MassAction');
% Pobj5f = addparameter(Kobj5, 'kf', kf_ptet*2);
% Pobj5r = addparameter(Kobj5, 'kr', kr_ptet);
% set(Kobj5, 'ParameterVariableNames', {'kf', 'kr'});
% 
% 
% Set up binding reaction for tetR
Robj2 = addreaction(tube, [dna.Name ':' P2 ' + ' RNAP ' <-> [' RNAPbound ':' P2 ']' ]);
Kobj2 = addkineticlaw(Robj2, 'MassAction');
Pobj2f = addparameter(Kobj2, 'kf', kf_ptet*2);
Pobj2r = addparameter(Kobj2, 'kr', kr_ptet);
set(Kobj2, 'ParameterVariableNames', {'kf', 'kr'});

% 
% decrease the kcat rate
kcat = log(2)/(rna.UserData/30);
kcat = kcat*0.001;
Rlist2 = txtl_rnap_rnap28(tube, dna, rna, [RNAPbound ':' P2 ],{P2},kcat);

 
Rlist = [Robj1, Rlist1];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
