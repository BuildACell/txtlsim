% txtl_protein_tetramerization-m - general protein tetramerization
% Zoltan A. Tuza Sep 2012
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

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

function Robj = txtl_protein_tetramerization(tube,protein,reactionRates)
% function for protein tetramerization.
% tube: sbiomodel object, where the reaction occurs
% protein: SimBiology Species Array
% reacctionRates: 2x1 or 1x2 vector contains the forward and reverse
% reaction rates. (For forward reaction only, set the reverse reaction rete to zero!)
%
% Return: SimBiology Reaction Array

if reactionRates(2) == 0
   Robj = addreaction(tube, ['2 ' protein.Name 'dimer ->' protein.Name 'tetramer']);
   Kobj = addkineticlaw(Robj,'MassAction');
   Pobj = addparameter(Kobj,  'kf', reactionRates(1));
   set(Kobj, 'ParameterVariableNames','kf');
else
   Robj = addreaction(tube, ['2 ' protein.Name 'dimer <->' protein.Name 'tetramer']); 
   Kobj = addkineticlaw(Robj,'MassAction');
   Pobjf = addparameter(Kobj, 'kf', reactionRates(1));
   Pobjr = addparameter(Kobj, 'kr', reactionRates(2));
   set(Kobj, 'ParameterVariableNames', {'kf', 'kr'});
end
   
end



