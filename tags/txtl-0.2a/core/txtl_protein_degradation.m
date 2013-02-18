% txtl_protein_degradation.m - general protein degradation model
% Zoltan A. Tuza Sep 2012
% Vipul Singhal, Oct 2012

% Protein Degradation using LVA tag and ClpXP. 
% reactionRates is an input argument that has 3 elements: 
% forward and backward binding rates for ClpXP and Protein, and the forward
% degradation rate of the ClpXP-Protein complex. 

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

function txtl_protein_degradation(mode, tube,protein,varargin)
% function for protein degradation.
% tube: sbiomodel object, where the reaction occurs
% protein: simBiology object
% reacctionRate: degration rate

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')
    %{
    a so-called LVA tag seemed to be the most efficient tag to make GFP unstable. 
    This tag consists of a short peptide sequence (AANDENYALVA) and is attached 
    to the C-terminal end of GFP.
    Source: http://2008.igem.org/Team:KULeuven/Data/GFP which cites:
    J.B. Andersen et al., �New Unstable Variants of Green Fluorescent Protein for 
    Studies of Transient Gene Expression in Bacteria,� Applied and Environmental 
    Microbiology, vol. 64, Jun. 1998, pp. 2240�2246.
    %}

    
    coreSpecies = {'ClpXP',[protein.Name ':ClpXP']};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)));
    
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode, 'Setup Reactions')

    % Protein monomer binds with ClpXP protease
    reactionRate = varargin{1};
    
    Robj = addreaction(tube, ['[' protein.Name '] + ClpXP <-> [' protein.Name ':ClpXP]']);
    Kobj = addkineticlaw(Robj,'MassAction');
    rN = regexprep(protein.Name, {'( )'}, {''});
    uniqueNameF = sprintf('TXTL_PROT_DEGRAD_COMPLEX_%s_F',rN);
    uniqueNameR = sprintf('TXTL_PROT_DEGRAD_COMPLEX_%s_R',rN);
    Pobjf = addparameter(Kobj, uniqueNameF, reactionRate(1));
    Pobjr = addparameter(Kobj, uniqueNameR, reactionRate(2));
    set(Kobj, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});
    

    % Degradation
       Robj2 = addreaction(tube, ['[' protein.Name ':ClpXP] -> ClpXP']);
       Kobj2 = addkineticlaw(Robj2,'MassAction');
       rN = regexprep(protein.Name, {'( )'}, {''});
       uniqueName = sprintf('TXTL_PROT_DEGRAD_%s',rN);
       Pobj2 = addparameter(Kobj2, uniqueName, reactionRate(3));
       set(Kobj2, 'ParameterVariableNames', {uniqueName});

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_protein_lacI:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 

end



