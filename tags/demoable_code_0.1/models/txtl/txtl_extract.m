% txtl_extract.m - function to create a tube of TX-TL extract
%! TODO: add documentation

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

function tube = txtl_extract(name)
tube = txtl_newtube(name);

% Add in ribosomes and RNAP70
%! TODO: update these numbers based on measurements
df = 1000;				% dilution factor of TX-TL mix
addspecies(tube, 'RNAP', 25/df);	% 25 nM based on simulac
sigma70 = addspecies(tube, 'protein sigma70', 25/df);	% 25 nM based on simulac
addspecies(tube, 'Ribo', 300/df);	% 300 nM based on simulac

% Add RNAP+Sigma70 <-> RNAP70 reaction
%! TODO: figure out the correct reaction rates; may need to implement
%    as quasi-steady state??  Original versions generated very stiff equations
% Kf = 1.7e6;% M^-1s^-1
% Kr = 4.3e-4; % s^-1
Kf = 100; Kr = 0.01;
% Set up the reaction
Robj1 = addreaction(tube, ['RNAP + ' sigma70.Name ' <-> RNAP70']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', Kf);
Pobj1r = addparameter(Kobj1, 'kr', Kr);
set(Kobj1, 'ParameterVariableNames', {'kf','kr'});

% Add in exonuclease + protection reactions (if [protein gamS] > 0)
%! TODO: update these numbers based on measurements
kgamS = 1;				% gamS binding rate
addspecies(tube, 'RecBCD', 25/df);	% 25 nM to match RNAP
Robj = addreaction(tube, 'RecBCD + [protein gamS] -> RecBCD:gamS');
Kobj = addkineticlaw(Robj,'MassAction');
Pobj = addparameter(Kobj, 'kf', kgamS);
set(Kobj, 'ParameterVariableNames', {'kf'});

% Add in RNA degradation
addspecies(tube, 'RNase', 25/df);	% 25 nM to match RNAP

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
