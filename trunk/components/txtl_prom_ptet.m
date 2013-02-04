% txtl_prom_ptet.m - promoter information for ptet promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_ptet() will set up the reactions for
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

function varargout = txtl_prom_ptet(mode, tube, dna, rna,varargin)

    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];
    % importing the corresponding parameters
    paramObj = txtl_component_config('tetR');


%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')
    
    promoterData = varargin{1};
    defaultBasePairs = {'ptet','junk','thio';...
        paramObj.Promoter_Length,paramObj.Junk_Length,paramObj.Thio_Length};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)));
  
    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
   
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'Setup Reactions')
    listOfSpecies = varargin{1};
    
    % Parameters that describe this promoter
    parameters = {'TXTL_PTET_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_PTET_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    %
    % nominal transcription
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);
    
    matchStr = regexp(listOfSpecies,'(^protein tetR.*dimer$)','tokens','once'); 
    listOftetRdimer = vertcat(matchStr{:});
    
   
    if ~isempty(listOftetRdimer)
        for k = 1:size(listOftetRdimer,1)
            txtl_addreaction(tube,...
                [DNA ' + aTc:' listOftetRdimer{k} ' <-> [' dna.name ':aTc:' listOftetRdimer{k} ']'],...
            'MassAction',{'ptet_atc_dna_F',0.2;...
                          'ptet_atc_dna_R',1}); 
            dnaaTcTetR = sbioselect(tube, 'Type', 'species', 'Name', [dna.name ':aTc:' listOftetRdimer{k}]);          
            txtl_transcription(mode, tube, dnaaTcTetR, rna, RNAP, RNAPbound);          
        end
    end

    %
    % Add reactions for sequestration of promoter by tetRdimer 
    %
    
 
    %! TODO make all these reactions conditional on specie availability
    
    if ~isempty(listOftetRdimer)
        for k = 1:size(listOftetRdimer,1)
            txtl_addreaction(tube,...
                [DNA ' + ' listOftetRdimer{k} ' <-> [' dna.name ':' listOftetRdimer{k} ']'],...
            'MassAction',{'ptet_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                          'ptet_sequestration_R',getDNASequestrationRates(paramObj,'R')});
            %{
            Robj10 = addreaction(tube, ...
              [DNA ' + ' listOftetR{i} ' <-> [' dna.name ':' listOftetR{i} '2]']);
            Kobj10 = addkineticlaw(Robj10,'MassAction');
            uniqueNameF = sprintf('TXTL_PTET_REPRESSION2_%s_F',rN);
            uniqueNameR = sprintf('TXTL_PTET_REPRESSION2_%s_R',rN);
            set(Kobj10, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});
            %}
            
        end
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
