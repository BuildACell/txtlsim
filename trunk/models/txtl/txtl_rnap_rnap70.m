% txtl_rnap_rnap70.m - reactions for RNAP with sigma70 bound
% RMM, 9 Sep 2012
%
% This file sets up the transcription reactions for RNAP bound to
% sigma70 ("RNAP70").  It can be called by promoter files that need to
% set up the approriate transcription reactions.

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

function Rlist = txtl_rnap_rnap70(tube, dna, rna, RNAPbound)
RNAP = 'RNAP70';			% RNA polymerase name for reactions

% Compute the number of NTPs required, in 100 NTP blocks
ntpcnt = floor(rna.UserData/100);	% get number of NTP blocks
if (ntpcnt == 0) 
  ntpstr = '';
else
  ntpstr = int2str(ntpcnt);
end

% Set up the transcription reaction
Robj1 = addreaction(tube, ...
  ['[' RNAPbound '] + ' ntpstr ' NTP <-> [NTP:' RNAPbound ']']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
set(Kobj1, 'ParameterVariableNames', {'TXTL_NTP_RNAP_F', 'TXTL_NTP_RNAP_R'});

Robj2 = addreaction(tube, ...
  ['[NTP:' RNAPbound '] -> ' dna.Name ' + ' rna.Name ' + ' RNAP]);
Kobj2 = addkineticlaw(Robj2, 'MassAction');
%generating unique parameter name for the current RNA
rN = regexprep(rna.Name, {'( )'}, {''});
uniqueName = sprintf('TXTL_TX_rate_%s',rN);
set(Kobj2, 'ParameterVariableNames', uniqueName);


Rlist = [Robj1, Robj2];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
% Parameters describing the enzymatic process
