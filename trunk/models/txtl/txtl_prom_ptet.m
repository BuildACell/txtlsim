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

txtl_transcription(tube, dna, rna, RNAP, RNAPbound);

%
% Add reactions for sequestration of promoter by tetRdimer 
%

%! TODO: RMM, 29 Sep 2012
%! TODO: txtl_protein_tetR defines tetramers, which aren't used
%! TODO: proper implementation for tetR is via two operator sites (I think)
% VS: yes, there are 2 operators, see for example, 
% C Berens and W. Hillen, Gene regulation by tetracyclines, Eur. J. Biochem. 270, 3109�3121 (2003)
%{
kf1_tetR = 0.2; kr1_tetR = 1;		% reaction rates (from sbio)
Robj4 = addreaction(tube, ...
  [DNA ' + [protein tetRdimer] <-> [' dna.name ':protein tetRdimer]']);
Kobj4 = addkineticlaw(Robj4,'MassAction');
Pobj4 = addparameter(Kobj4, 'k4', kf1_tetR);
Pobj4r = addparameter(Kobj4, 'k4r', kr1_tetR);
set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

kf2_tetR = 0.2; kr2_tetR = 1;		
Robj5 = addreaction(tube, ...
  ['[' dna.name ':protein tetRdimer] + [protein tetRdimer] <-> [' dna.name ':protein tetRdimer:protein tetRdimer]']);
Kobj5 = addkineticlaw(Robj5,'MassAction');
Pobj5 = addparameter(Kobj5, 'k5', kf2_tetR);
Pobj5r = addparameter(Kobj5, 'k5r', kr2_tetR);
set(Kobj5, 'ParameterVariableNames', {'k5', 'k5r'});

kf3_tetR = 0.2; kr3_tetR = 1;		% reaction rates (from sbio)
Robj7 = addreaction(tube, ...
  [DNA ' + [protein tetR-lvadimer] <-> [' dna.name ':protein tetR-lvadimer]']);
Kobj7 = addkineticlaw(Robj7,'MassAction');
Pobj7 = addparameter(Kobj7, 'k7', kf3_tetR);
Pobj7r = addparameter(Kobj7, 'k7r', kr3_tetR);
set(Kobj7, 'ParameterVariableNames', {'k7', 'k7r'});

kf4_tetR = 0.2; kr4_tetR = 1;		
Robj6 = addreaction(tube, ...
  ['[' dna.name ':protein tetR-lvadimer] + [protein tetR-lvadimer] <-> [' dna.name ':protein tetR-lvadimer:protein tetR-lvadimer]']);
Kobj6 = addkineticlaw(Robj6,'MassAction');
Pobj6 = addparameter(Kobj6, 'k6', kf4_tetR);
Pobj6r = addparameter(Kobj6, 'k6r', kr4_tetR);
set(Kobj6, 'ParameterVariableNames', {'k6', 'k6r'});
%}
kf5_tetR = 3; kr5_tetR = 0.5;		% reaction rates
Robj8 = addreaction(tube, ...
  [DNA ' + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer1]']);
Kobj8 = addkineticlaw(Robj8,'MassAction');
Pobj8 = addparameter(Kobj8, 'k8', kf5_tetR);
Pobj8r = addparameter(Kobj8, 'k8r', kr5_tetR);
set(Kobj8, 'ParameterVariableNames', {'k8', 'k8r'});

kf6_tetR = 0.000005; kr6_tetR = 0.000005;
Robj10 = addreaction(tube, ...
  [DNA ' + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer2]']);
Kobj10 = addkineticlaw(Robj10,'MassAction');
Pobj10 = addparameter(Kobj10, 'k10', kf6_tetR);
Pobj10r = addparameter(Kobj10, 'k10r', kr6_tetR);
set(Kobj10, 'ParameterVariableNames', {'k10', 'k10r'});

kf7_tetR = 0.0000006; kr7_tetR = 0.0000005;		% effectively nonexistent
Robj9 = addreaction(tube, ...
  ['[' dna.name ':protein tetR-lva-terminatordimer1] + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer]']);
Kobj9 = addkineticlaw(Robj9,'MassAction');
Pobj9 = addparameter(Kobj9, 'k9', kf7_tetR);
Pobj9r = addparameter(Kobj9, 'k9r', kr7_tetR);
set(Kobj9, 'ParameterVariableNames', {'k9', 'k9r'});
	
kf8_tetR = 0.0000006; kr8_tetR = 0.0000005;
Robj11 = addreaction(tube, ...
  ['[' dna.name ':protein tetR-lva-terminatordimer2] + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer]']);
Kobj11 = addkineticlaw(Robj11,'MassAction');
Pobj11 = addparameter(Kobj11, 'k11', kf8_tetR);
Pobj11r = addparameter(Kobj11, 'k11r', kr8_tetR);
set(Kobj11, 'ParameterVariableNames', {'k11', 'k11r'});

Rlist = [Robj1, Robj8, Robj9, Robj10, Robj11];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
