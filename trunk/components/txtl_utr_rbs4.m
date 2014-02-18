% txtl_utr_rbs4.m - promoter information for variant 4 of the rbs
% VS Dec 2013
%
% This file contains a description of a ribosome binding site.
% It has a different ribosome binding affinity.

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

function varargout = txtl_utr_rbs4(mode, tube, rna, protein, varargin)
RBSnum = '4';
ribo_f = 0.000005;
ribo_r = 1;
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')

    utrRbsData = varargin{1};
    defaultBasePairs = {['rbs' RBSnum],'spacer'; 20,200};
    utrRbsData = txtl_setup_default_basepair_length(tube,utrRbsData,defaultBasePairs);
    

    RiboBound = ['Ribo:' rna.Name];
    coreSpecies = {'Ribo',RiboBound};
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    varargout{1} = sbioselect(tube, 'Name', RiboBound);
    varargout{2} = utrRbsData;

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')

    listOfSpecies = get(tube.species, 'name');
    
    % Set up the binding reaction
        txtl_addreaction(tube,['[' rna.Name '] + Ribo <-> [Ribo:' rna.Name ']'],...
            'MassAction',{['TXTL_UTR_RBS' RBSnum '_F'],ribo_f*tube.UserData.ReactionConfig.Ribosome_Binding_F;
            ['TXTL_UTR_RBS' RBSnum '_R'],ribo_r*tube.UserData.ReactionConfig.Ribosome_Binding_R});

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error(['txtltoolbox:txtl_utr_rbs' RBSnum ':undefinedmode'], ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
