% txtl_prom_placI.m - promoter information for lacI promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_placI() will set up the reactions for
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

function [Rlist, promlen] = txtl_prom_placI(tube, dna, rna, promFull, promlen)

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
            case 'placI'
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
kf_placI = log(2)/0.1;			% 100 ms bind rate
kr_placI = 10 * kf_placI;			% Km of 10 (same as p70, from VN)
ktx_placI = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
RNAPbound = ['RNAP70:' dna.Name];

% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_placI);
Pobj1r = addparameter(Kobj1, 'kr', kr_placI);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

%
% Now put in the reactions for the utilization of NTPs
% Use an enzymatic reaction to proper rate limiting
%

Rlist1 = txtl_rnap_rnap70(tube, dna, rna, RNAPbound);

%
%% Add reactions for sequestration of promoter by lacI 
%
% bind lacI tetramer to one site of the dna:
% DNA + tetramer <-> DNA:tetramer
% then bind the tetramer to the other site on the dna
% DNA:tetramer <-> DNA:tetramerLoop

%there are two sites where the tetramer binds: Om and Oa (see Alberts p437,
%5th ed). We model these as equivalent for now. 

%bind tetramer to a single side
kf_lacI1 = 5; kr_lacI1 = 1;		
Robj4 = addreaction(tube, ...
  ['[' dna.Name '] + [protein lacItetramer] <-> ' ...
   '[' dna.Name ':protein lacItetramer]']);
Kobj4 = addkineticlaw(Robj4,'MassAction');
Pobj4 = addparameter(Kobj4, 'k4', kf_lacI1);
Pobj4r = addparameter(Kobj4, 'k4r', kr_lacI1);
set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

%DNA looping
kf_lacI2 = 20; kr_lacI2 = 1;		
Robj5 = addreaction(tube, ...
  ['[' dna.Name ':protein lacItetramer] <-> ' ...
  '[' dna.Name ':protein lacItetramerLoop]']);
Kobj5 = addkineticlaw(Robj5,'MassAction');
Pobj5 = addparameter(Kobj5, 'k5', kf_lacI2);
Pobj5r = addparameter(Kobj5, 'k5r', kr_lacI2);
set(Kobj5, 'ParameterVariableNames', {'k5', 'k5r'});

%{
% bind dimers individually, ad then tetramerize (we set reaction rates to 0
% for these first. Can be activate by changing rr's to nonzero values. 

% single dimer binding to the DNA
kf_lacI3 = 0; kr_lacI3 = 0;		
Robj6 = addreaction(tube, ...
  ['[' dna.Name '] + [protein lacIdimer] <-> ' ...
   '[' dna.Name ':protein lacIdimer]']);
Kobj6 = addkineticlaw(Robj6,'MassAction');
Pobj6 = addparameter(Kobj6, 'k6', kf_lacI3);
Pobj6r = addparameter(Kobj6, 'k6r', kr_lacI3);
set(Kobj6, 'ParameterVariableNames', {'k6', 'k6r'});

% 2 dimers bound individually to the DNA
kf_lacI4 = 0; kr_lacI4 = 0;		
Robj7 = addreaction(tube, ...
  ['[' dna.Name ':protein lacIdimer] + [protein lacIdimer] <-> ' ...
   '[' dna.Name ':protein lacIdimer:protein lacIdimer]']);
Kobj7 = addkineticlaw(Robj7,'MassAction');
Pobj7 = addparameter(Kobj7, 'k7', kf_lacI4);
Pobj7r = addparameter(Kobj7, 'k7r', kr_lacI4);
set(Kobj7, 'ParameterVariableNames', {'k7', 'k7r'});

% DNA bound to two dimer -> dimers joining to form tetramer-DNA loop
kf_lacI5 = 0; kr_lacI5 = 0;		
Robj8 = addreaction(tube, ...
  [ '[' dna.Name ':protein lacIdimer:protein lacIdimer] <-> ' ...
   '[' dna.Name ':protein lacItetramerLoop]']);
Kobj8 = addkineticlaw(Robj8,'MassAction');
Pobj8 = addparameter(Kobj8, 'k8', kf_lacI5);
Pobj8r = addparameter(Kobj8, 'k8r', kr_lacI5);
set(Kobj8, 'ParameterVariableNames', {'k8', 'k8r'});
%}

Rlist = [Robj1, Rlist1, Robj4, Robj5];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
