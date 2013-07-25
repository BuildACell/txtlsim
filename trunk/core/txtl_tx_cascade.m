function varargout = txtl_tx_cascade(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec)
% txtl_tx_cascade.m - RNA circuit set-up file
% VS 7-25-2013
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

% cant do this because it fails in reaction mode. 
% prom_spec = tube.userdata.DNAinfo{end}{1};
% rbs_spec = tube.userdata.DNAinfo{end}{2};
% gene_spec = tube.userdata.DNAinfo{end}{3};

% get the various strings specified in the promoter, utr, gene
[promData, promStr] = txtl_parsespec(prom_spec);
[utrData, utrStr] = txtl_parsespec(rbs_spec);
[geneData, geneStr] = txtl_parsespec(gene_spec);

att = utrData{1,1};
if ~strcmp(att, {'att1','att2'}) %can be extended to any number of att's
    error('Something went wrong. We think there is Attenuator RNA present');
end

%% Setup Species
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    if nargin == 9
        coreSpecies = {'NTP',RNAPbound,['NTP:' RNAPbound], ['RNA ' att], RNAP};
    else
        error('the number of argument should be at 9, not %d',nargin);
    end
    varargout{1} = coreSpecies;
%% Setup Reactions
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
  
    
    % Set up tx reaction string. reactions themselves later
    if nargin ==9
                transcriptionEq1 = ...
            ['[NTP:' RNAPbound '] -> [' RNAPbound ':RNA ' att ']']; 
        transcriptionEq2 = ...
            ['[NTP:' RNAPbound ':RNA ' att '] -> '  dna.Name ' + ' rna.Name ' + ' RNAP];
    else
        error('the number of argument should be at 9, not %d',nargin);
    end
    ktxExpression1 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
        'RNA_Length','rna.UserData.att'); % rna.UserData.att defined in txtl_add_dna
    ktxExpression2 =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
        'RNA_Length','rna.UserData.remaining'); % rna.UserData.remaining defined in txtl_add_dna
    ktx1 = eval(ktxExpression1);
    ktx2 = eval(ktxExpression2);
    
    % binding NTP
    NTPparameters = {'TXTL_NTP_RNAP_F', tube.UserData.ReactionConfig.NTP_Forward;
        'TXTL_NTP_RNAP_R', tube.UserData.ReactionConfig.NTP_Reverse};
    txtl_addreaction(tube,['[' RNAPbound '] + NTP <-> [NTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);
    txtl_addreaction(tube,['[' RNAPbound ':RNA ' att '] + NTP <-> [NTP:' RNAPbound ':RNA ' att ']'],...
        'MassAction',NTPparameters);

    %dummy reaction for NTP consumption
    ntpcnt1 = ceil(rna.UserData.att/100);	% get number of NTP blocks
    ntpcnt2 = ceil(rna.UserData.remaining/100);	% get number of NTP blocks
    NTPConsumptionRate1 = {'TXTL_NTP_consumption1',(ntpcnt1-1)*ktx1};
    NTPConsumptionRate2 = {'TXTL_NTP_consumption2',(ntpcnt2-1)*ktx2};
    txtl_addreaction(tube,['[NTP:' RNAPbound '] -> ' dna.Name ' +  ' RNAP],...
        'MassAction',NTPConsumptionRate1);
    txtl_addreaction(tube,['[NTP:' RNAPbound ':RNA ' att '] -> ' dna.Name ' +  ' RNAP ' + RNA ' att ],...
        'MassAction',NTPConsumptionRate2);
    
    %re bind RNAP to DNA (and RNA att) --- check this with Richard
    promoterName = promData{1,end};
    paramObj = txtl_component_config(promoterName);
    parameters = {'TXTL_P70_RNAPbound_F',paramObj.RNAPbound_Forward;...
        'TXTL_P70_RNAPbound_R',paramObj.RNAPbound_Reverse};
    txtl_addreaction(tube,[dna.Name ' + ' RNAP ' + [RNA ' att  '] <-> ' RNAPbound ':RNA ' att ],...
        'MassAction',parameters);
    
    %att anti reactions
    attnumber = att(end);
    [~,listOfSpecies] = getstoichmatrix(tube);
    matchStr = regexp(listOfSpecies,['(^RNA .*anti' attnumber '$)'],'tokens','once'); 
    listOfantisense = vertcat(matchStr{:});
    
    complex1_rate = {'TXTL_RNA_C1_F', tube.UserData.ReactionConfig.complex1_F_rate; ...
                'TXTL_RNA_C1_R', tube.UserData.ReactionConfig.complex1_R_rate};
    att_anti_termination_rate = {'TXTL_RNA_ATTANTI_TERM', tube.UserData.ReactionConfig.att_anti_termination_rate};
    
    if ~isempty(listOfantisense)
        for k = 1:size(listOfantisense,1)
            txtl_addreaction(tube,...
                [RNAPbound ':RNA ' att ' + ' listOfantisense{k} ' <-> [' RNAPbound ':RNA ' att ':' listOfantisense{k} ']'],...
            'MassAction',complex1_rate);
            txtl_addreaction(tube,['[' RNAPbound ':RNA ' att ':' listOfantisense{k} '] -> ' dna.Name ' + [' listOfantisense{k} ':RNA ' att '] + ' RNAP],...
                'MassAction',att_anti_termination_rate);
        end
    end    
    
    % transcription
    txtl_addreaction(tube,transcriptionEq1,'MassAction',{'TXTL_transcription_rate1',ktx1});
    txtl_addreaction(tube,transcriptionEq2,'MassAction',{'TXTL_transcription_rate2',ktx2});

    % Auto-Termination
    auto_termination_rate = {'TXTL_RNA_ATT_AUTOTERM', tube.UserData.ReactionConfig.auto_termination_rate};
    txtl_addreaction(tube,['[NTP:' RNAPbound ':RNA ' att '] -> ' dna.Name ' + [RNA ' att '] + ' RNAP],...
        'MassAction',auto_termination_rate);
    
    
else
    error('txtltoolbox:txtl_tx_cascade:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end


