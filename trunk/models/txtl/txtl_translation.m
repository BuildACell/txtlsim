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


function [] = txtl_translation(mode, tube, dna, rna, protein, Ribobound)



%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')

     % Set up the species for translation 
    coreSpecies = {'AA',['AA:' Ribobound.Name],'Ribo','RNase'};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)));
    
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode, 'Setup Reactions')
    
    AA_model = 1;
    if AA_model == 1
        %tube.UserData.AAmodel == 1

        aacnt = floor(protein.UserData/100);	% get number of K amino acids
        if (aacnt == 0) 
          aastr = '';
        else
          aastr = int2str(aacnt);
        end
        Robj = addreaction(tube, ...
          ['[' Ribobound.Name '] + ' aastr ' AA <-> [AA:' Ribobound.Name ']']);
        Kobj = addkineticlaw(Robj, 'MassAction');
        set(Kobj, 'ParameterVariableNames', {'TXTL_AA_F', 'TXTL_AA_R'});
    else
        Robj = addreaction(tube, ...
          ['[' Ribobound.Name '] + AA <-> [AA:' Ribobound.Name ']']);
        Kobj = addkineticlaw(Robj, 'MassAction');
        set(Kobj, 'ParameterVariableNames', {'TXTL_AA_F', 'TXTL_AA_R'});

        Robj3 = addreaction(tube, ...
         ['[AA:' Ribobound.Name '] -> ' rna.Name ' +  Ribo']);
        Kobj3 = addkineticlaw(Robj3, 'MassAction');
        %generate unique parameter name for the current protein
        rN = regexprep(protein.Name, {'( )'}, {''});
        uniqueName = sprintf('TXTL_TL_rate_%s_AA_consumption',rN);
        set(Kobj3, 'ParameterVariableNames', uniqueName);


    end
    Robj1 = addreaction(tube, ...
      ['[AA:' Ribobound.Name '] -> ' rna.Name ' + ' protein.Name ' +  Ribo']);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    %generate unique parameter name for the current protein
    rN = regexprep(protein.Name, {'( )'}, {''});
    uniqueName = sprintf('TXTL_TL_rate_%s',rN);
    set(Kobj1, 'ParameterVariableNames', uniqueName);
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_protein_lacI:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    
    

end