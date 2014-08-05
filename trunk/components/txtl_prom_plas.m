% txtl_prom_plas.m - promoter information for plas
% VS Dec 2013
%
% This file contains a description of the plasR promoter, which is activated by the lasR protein. .
% Calling the function txtl_prom_plas() will set up the reactions for
% transcription with the measured binding rates and transription rates.
%
%
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

function varargout= txtl_prom_plas(mode, tube, dna, rna, varargin)


% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
RNAPbound = ['RNAP70:' dna.Name];
P1 = 'protein sigma70';

%P3 = 'OC12HSL:protein lasR'; %wont yet know if the lva version of this
%was added or not.


% importing the corresponding parameters
paramObj = txtl_component_config('plas');

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    promoterData = varargin{1};
    if nargin==8
        prom_spec = varargin{2};
        rbs_spec = varargin{3};
        gene_spec = varargin{4};
    elseif nargin~=5
        error('the number of argument should be 5 or 8, not %d',nargin);
    end
    defaultBasePairs = {'plas','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,P1};
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    %{
    % !TODO: VS Aug 2014: dont yet know what form of P3 exists: is lasR lva/ssrA bound or not?
    % (and in the future more variations could exist) So we cannot actually
    % do txtl_transcription just yet! Functionally this is not a problem
    % because: Lets say we do not introduce the activated RNAPbound at this
    % (Setup Species) stage. because P3, which is the protein bound to the
    % inducer, *will* get set up in the 'setup species' of the protein file
    % (if the file is present and the protein SHOULD be set up), and so
    % will be available for the setup reactions stage to search and call in
    % the transcription step. If there is no protein file, then the search
    % returns no results, and so the transcription is not set up, or only
    % leaky transcription is set up. One thing to note is that, to set up
    % the species (the AGTP:CUTP:RNAPbound species), might need to call 'setup species' mode of
    % txtl_transcription in the Setup Reactions mode of this promoter file!
    % But maybe if we dont set up that species no error will be trown.
    % We'll see.
    
    % VS Aug 2014: Comment txtl_transcription for now.
    %     if mode.utr_attenuator_flag
    % %         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec );
    %         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P3 ],prom_spec, rbs_spec, gene_spec,{P3} );
    %     else
    % %         txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
    %         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3});
    %     end
    %}
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    listOfSpecies = varargin{1};
    if nargin==8
        prom_spec = varargin{2};
        rbs_spec = varargin{3};
        gene_spec = varargin{4};
    elseif nargin~=5
        error('the number of argument should be 5 or 8, not %d',nargin);
    end
    
    %   Parameters that describe this promoter. No leaky expression for now.
    %     parameters = {'TXTL_PLAS_RNAPbound_F',paramObj.RNAPbound_Forward;...
    %                   'TXTL_PLAS_RNAPbound_R',paramObj.RNAPbound_Reverse};
    %     txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
    %             'MassAction',parameters);
    
    p = regexp(listOfSpecies,'^OC12HSL:protein lasR(-lva)?$', 'match');
    TF = vertcat(p{:});
    for k = 1:size(TF,1)
        setup_promoter_reactions(mode, tube, dna, rna, RNAP,RNAPbound,...
            prom_spec, rbs_spec, gene_spec,TF{k}, paramObj)
    end
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_plas:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end
end


function setup_promoter_reactions(mode, tube, dna,rna, RNAP, RNAPbound,...
    prom_spec, rbs_spec, gene_spec, TF, paramObj)

txtl_addreaction(tube, ...
    [dna.Name ' + ' TF ' <-> [' dna.Name ':' TF ']' ],...
    'MassAction',{'TXTL_PLAS_TFBIND_F',paramObj.activation_F;...
    'TXTL_PLAS_TFBIND_R',paramObj.activation_R});

txtl_addreaction(tube, ...
    [dna.Name ':' TF ' + ' RNAP ' <-> [' RNAPbound ':' TF ']' ],...
    'MassAction',{'TXTL_PLAS_TFRNAPbound_F',paramObj.RNAPbound_Forward_actv;...
    'TXTL_PLAS_TFRNAPbound_R',paramObj.RNAPbound_Reverse_actv});


if mode.utr_attenuator_flag
    mode.add_dna_driver = 'Setup Species';
    txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
    mode.add_dna_driver = 'Setup Reactions';
    txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
else
    mode.add_dna_driver = 'Setup Species';
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
    mode.add_dna_driver = 'Setup Reactions';
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
end

end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
