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


function txtl_translation(mode, tube, dna, rna, protein, Ribobound)
% basically, this function does nothing if the RNA does not have a ribosome
% binding site. so for instance if the RNA is: [RNA att] or [RNA anti].


%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
     % Set up the species for translation 
    coreSpecies = {'AA',['AA:AGTP:' Ribobound.Name],'Ribo'};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
    AAparameters = {'TXTL_AA_F',tube.UserData.ReactionConfig.AA_Forward;
                  'TXTL_AA_R',tube.UserData.ReactionConfig.AA_Reverse};
    
    % translation rate             
    ktlExpression =  strrep(tube.UserData.ReactionConfig.Translation_Rate,...
            'Protein_Length','protein.UserData');             
    ktl_rbs = eval(ktlExpression);              
              
    % AA consumption models              
    if tube.UserData.ReactionConfig.AAmodel == 1

        aacnt = floor(protein.UserData/100);  % get number of K amino acids
        if (aacnt == 0) 
          aastr = '';
        else
          aastr = int2str(aacnt);
        end
        
        txtl_addreaction(tube,...
            ['[' Ribobound.Name '] + ' aastr ' AA <-> [AA:' Ribobound.Name ']'],...
            'MassAction',AAparameters);
    else
        
        txtl_addreaction(tube, ...
            ['[' Ribobound.Name '] + AA + AGTP <-> [AA:AGTP:' Ribobound.Name ']'],...
            'MassAction',AAparameters);
        
        aacnt = floor(protein.UserData/100);
        aa_consump_rate = (aacnt-1)*ktl_rbs;
        
        txtl_addreaction(tube, ...
            ['[AA:AGTP:' Ribobound.Name '] -> ' rna.Name ' +  Ribo + ATP'],...
            'MassAction',{'TXTL_TL_AA_consumption',aa_consump_rate});
    end
    
    % Translation
    txtl_addreaction(tube, ...
     ['[AA:AGTP:' Ribobound.Name '] -> ' rna.Name ' + ' protein.Name ' +  Ribo'],...
     'MassAction',{'TXTL_TL_rate',ktl_rbs});
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_translation:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    
    

end