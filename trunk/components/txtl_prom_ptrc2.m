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

function varargout = txtl_prom_ptrc2(mode, tube, dna, rna, varargin)

    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];
    % importing the corresponding parameters
    paramObj = txtl_component_config('lacI');    

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')

    
    promoterData = varargin{1};
    defaultBasePairs = {'placI','junk','thio';...
        paramObj.Promoter_Length,paramObj.Junk_Length,paramObj.Thio_Length};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');

    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);



%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    listOfSpecies = varargin{1};
    
    % Parameters that describe this promoter
    parameters = {'TXTL_PTRC2_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_PTRC2_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);

    % nominal transcription
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);

    
    matchStr = regexp(listOfSpecies,'(^protein lacI.*dimer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
    listOflacIdimers = vertcat(matchStr{:});
    % repression of ptrc2 by lacIdimer
    if ~isempty(listOflacIdimers)
        for i = 1:size(listOflacIdimers,1)
             txtl_addreaction(tube,...
                [DNA ' + ' listOflacIdimers{i} ' <-> [' dna.name ':' listOflacIdimers{i} ']'],...
            'MassAction',{'ptrc2_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                          'ptrc2_sequestration_R',getDNASequestrationRates(paramObj,'R')});
        end
    end

       % Un-repression of lacIdimer repressed ptrc2 in the presence of IPTG inducer.
    if ~isempty(listOflacIdimers)
        for k = 1:size(listOflacIdimers,1)
            txtl_addreaction(tube,...
                [DNA ' + IPTG:' listOflacIdimers{k} ' <-> [' dna.name ':IPTG:' listOflacIdimers{k} ']'],...
            'MassAction',{'ptet_IPTG_dna_F',0.2;...
                          'ptet_IPTG_dna_R',1}); 
            dnaIPTGlacI = sbioselect(tube, 'Type', 'species', 'Name', [dna.name ':IPTG:' listOflacIdimers{k}]);   
            txtl_addreaction(tube,[dnaIPTGlacI.name ' + ' RNAP ' <-> [' RNAP ':' dnaIPTGlacI.name ']'],...
        'MassAction',parameters);     
            txtl_transcription(mode, tube, dnaIPTGlacI, rna, RNAP, [RNAP ':' dnaIPTGlacI.name]);          
        end
    end

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_ptrc2:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
