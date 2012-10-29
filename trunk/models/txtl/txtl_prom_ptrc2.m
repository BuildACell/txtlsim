% txtl_prom_ptrc2.m - promoter information for ptrc2 promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the ptrc2 promoter.
% Calling the function txtl_prom_ptrc2() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% The binding of the promoter to the tetR repressor is used in the
% gen_switch example. 

% VS Sep 2012
% Adapted from Richard Murray's original code. 
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

function [Rlist, promlen] = txtl_prom_ptrc2(tube, dna, rna, promFull, promlen)

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
            case 'ptrc2'
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
kf_ptrc2 = log(2)/0.1;			% 100 ms bind rate
kr_ptrc2 = 10 * kf_ptrc2;			% Km of 10 (same as p70, from VN)
ktx_ptrc2 = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
RNAPbound = ['RNAP70:' dna.Name];


% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_ptrc2);
Pobj1r = addparameter(Kobj1, 'kr', kr_ptrc2);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

%
% Now put in the reactions for the utilization of NTPs
% Use an enzymatic reaction to proper rate limiting
%

Rlist1 = txtl_rnap_rnap70(tube, dna, rna, RNAPbound);

%
% Add reactions for sequestration of promoter by lacIdimer and lacItetramer
% Ptrc-2 only has 1 operator site: Olac. I think both lacIdimer and lacItetramer
% should be able to bind to this. -VS, 10/22
% See supplementary info in the genetic toggle switch paper, gardener and
% collins, 2000. 

% ** IMPORTANT issue in implementation and modularity **
% We avoid selecting for a protein species containing lacI (and tetR in 
% the ptet file), using sbioselect as commented out below because of
% the order in which the lacI and tetR proteins are declared: in
% gen_switch.m, ptrc2 is declared AFTER lacI has been created, so lacI
% exists in the tube, and sbioselect can find it and its dimers. But if
% this was tried for tetR and ptet, it would fail because ptet is declared
% BEFORE tetR is created. So for now we avoid the sbioselect based search using 
% regular expressions, and just list all the 3 possibilities. In the
% future, when all the species and reactions can be declared in an external
% file, these problems will be taken care of, since all the species will be
% declared simultaneously. 

%{
Obj = sbioselect (tube, 'Where',...
 'Name', 'regexp', 'protein ' && 'lacI')
protein = sbioselect(tube, 'Type', 'species', 'Name', protstr);
PROTEIN = ['[' protein.Name ']']; %PROTEIN species name for reactions
%}
%{
kf1_lacI = 1; kr1_lacI = 0.1;		% 
Robj4 = addreaction(tube, ...
  [DNA ' + [protein lacIdimer] <-> [' dna.name ':protein lacIdimer]']);
Kobj4 = addkineticlaw(Robj4,'MassAction');
Pobj4 = addparameter(Kobj4, 'k4', kf1_lacI);
Pobj4r = addparameter(Kobj4, 'k4r', kr1_lacI);
set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

kf2_lacI = 1; kr2_lacI = 0.1;		% 
Robj5 = addreaction(tube, ...
  [DNA ' + [protein lacItetramer] <-> [' dna.name ':protein lacItetramer]']); 
Kobj5 = addkineticlaw(Robj5,'MassAction');
Pobj5 = addparameter(Kobj5, 'k5', kf2_lacI);
Pobj5r = addparameter(Kobj5, 'k5r', kr2_lacI);
set(Kobj5, 'ParameterVariableNames', {'k5', 'k5r'});


kf3_lacI = 1; kr3_lacI = 0.1;		% 
Robj6 = addreaction(tube, ...
  [DNA ' + [protein lacI-lvadimer] <-> [' dna.name ':protein lacI-lvadimer]']);
Kobj6 = addkineticlaw(Robj6,'MassAction');
Pobj6 = addparameter(Kobj6, 'k6', kf3_lacI);
Pobj6r = addparameter(Kobj6, 'k6r', kr3_lacI);
set(Kobj6, 'ParameterVariableNames', {'k6', 'k6r'});

kf4_lacI = 1; kr4_lacI = 0.1;		% 
Robj7 = addreaction(tube, ...
  [DNA ' + [protein lacI-lvatetramer] <-> [' dna.name ':protein lacI-lvatetramer]']); 
Kobj7 = addkineticlaw(Robj7,'MassAction');
Pobj7 = addparameter(Kobj7, 'k7', kf4_lacI);
Pobj7r = addparameter(Kobj7, 'k7r', kr4_lacI);
set(Kobj7, 'ParameterVariableNames', {'k7', 'k7r'});
%}
kf5_lacI = 0.50; kr5_lacI = 0.03;		% 
Robj8 = addreaction(tube, ...
  [DNA ' + [protein lacI-lva-terminatordimer] <-> [' dna.name ':protein lacI-lva-terminatordimer]']);
Kobj8 = addkineticlaw(Robj8,'MassAction');
Pobj8 = addparameter(Kobj8, 'k8', kf5_lacI);
Pobj8r = addparameter(Kobj8, 'k8r', kr5_lacI);
set(Kobj8, 'ParameterVariableNames', {'k8', 'k8r'});

kf6_lacI = 0.00000012; kr6_lacI = 0.000005;		% effectively ignored. 
Robj9 = addreaction(tube, ...
  [DNA ' + [protein lacI-lva-terminatortetramer] <-> [' dna.name ':protein lacI-lva-terminatortetramer]']); 
Kobj9 = addkineticlaw(Robj9,'MassAction');
Pobj9 = addparameter(Kobj9, 'k9', kf6_lacI);
Pobj9r = addparameter(Kobj9, 'k9r', kr6_lacI);
set(Kobj9, 'ParameterVariableNames', {'k9', 'k9r'});

Rlist = [Robj1, Rlist1, Robj8, Robj9];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
