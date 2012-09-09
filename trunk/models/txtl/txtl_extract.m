% txtl_extract.m - function to create a tube of TX-TL extract
%! TODO: add documentation

% Written by Richard Murray, Sep 2012
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
addspecies(tube, 'RNAP70', 100);
addspecies(tube, 'Ribo', 100);

% Add in exonuclease + protection reactions (if [protein gamS] > 0)
kgamS = 1;				% gamS binding rate
addspecies(tube, 'RecBCD', 100);
Robj = addreaction(tube, 'RecBCD + [protein gamS] -> RecBCD:gamS');
Kobj = addkineticlaw(Robj,'MassAction');
Pobj = addparameter(Kobj, 'kf', kgamS);
set(Kobj, 'ParameterVariableNames', {'kf'});

% Add in RNA degradation
addspecies(tube, 'RNase', 100);
