% txtl_transcription.m - sigma factor independent implementation for gene
% transcription in the TXTL system
% RMM, 9 Sep 2012
%
% It can be called by promoter files that need to
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

function txtl_transcription(mode, varargin)
tube = varargin{1};
dna = varargin{2};
rna = varargin{3};
RNAP = varargin{4}; % RNA polymerase name for reactions
RNAPbound = varargin{5};

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    if ~mode.utr_attenuator_flag
        if nargin < 6
            error('the number of argument should be at least 6, not %d',nargin);
        elseif nargin > 6
            extraSpecies = varargin{6};
            coreSpecies = {'NTP',RNAPbound,['NTP:' RNAPbound],RNAP,extraSpecies{:}};
        else
            coreSpecies = {'NTP',RNAPbound,['NTP:' RNAPbound],RNAP};
        end
    else
        if nargin < 6
            error('the number of argument should be at least 6, not %d',nargin);
        elseif nargin > 6
            extraSpecies = varargin{6};
            coreSpecies = {'NTP',RNAPbound,['NTP:' RNAPbound],RNAP, 'RNA att', extraSpecies{:}};
        else
            coreSpecies = {'NTP',RNAPbound,['NTP:' RNAPbound], 'RNA att', RNAP};
        end
    end
    
    
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    if ~mode.utr_attenuator_flag
        if nargin < 6
            error('the number of argument should be at least 6, not %d',nargin);
        elseif nargin > 6
            extraSpecies = varargin{6};
            % processing the extraSpecies
            extraStr = extraSpecies{1};
            for k=2:size(extraSpecies,1)
                extraStr = [extraStr '+' extraSpecies{k}];
            end
            transcriptionEq = ...
                ['[NTP:' RNAPbound '] -> ' dna.Name ' + ' rna.Name ' + ' RNAP ' + ' extraStr];
        else
            transcriptionEq = ...
                ['[NTP:' RNAPbound '] -> ' dna.Name ' + ' rna.Name ' + ' RNAP];
        end
        ktxExpression =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
            'RNA_Length','rna.UserData');
        ktx = eval(ktxExpression);
        % parameter values
        NTPparameters = {'TXTL_NTP_RNAP_F', tube.UserData.ReactionConfig.NTP_Forward;
            'TXTL_NTP_RNAP_R', tube.UserData.ReactionConfig.NTP_Reverse};
    else
                
        if nargin < 6
            error('the number of argument should be at least 7, not %d',nargin);
        elseif nargin > 6
            extraSpecies = varargin{6};
            % processing the extraSpecies
            extraStr = extraSpecies{1};
            for k=2:size(extraSpecies,1)
                extraStr = [extraStr '+' extraSpecies{k}];
            end
            transcriptionEq1 = ...
                ['[NTP:' RNAPbound '] -> [' RNAPbound ':RNA att]']; % att, when present, is always the first element in the RNA.
            transcriptionEq2 = ...
                ['[NTP:' RNAPbound ':RNA att] -> '  dna.Name ' + ' rna.Name ' + ' RNAP ' + ' extraStr];
        else
            transcriptionEq1 = ...
                ['[NTP:' RNAPbound '] -> [' RNAPbound ':RNA att]']; % att, when present, is always the first element in the RNA.
            transcriptionEq2 = ...
                ['[NTP:' RNAPbound ':RNA att] -> '  dna.Name ' + ' rna.Name ' + ' RNAP];
        end
        %replace the string for RNA_length variable in the expression with the
        %string for the actual rna length. Tx rate the same for both
        %reactions. can be made different if needed. 
        ktxExpression =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
            'RNA_Length','rna.UserData');
        ktx1 = eval(ktxExpression);
        ktxExpression =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
            'RNA_Length','rna.UserData');
        ktx2 = eval(ktxExpression);
        % parameter values for the binding of NTP to RNAPbound or 
        % [RNAPbound ':RNA att']. We assume same rate. can be made
        % different
        NTPparameters = {'TXTL_NTP_RNAP_F', tube.UserData.ReactionConfig.NTP_Forward;
            'TXTL_NTP_RNAP_R', tube.UserData.ReactionConfig.NTP_Reverse}; 
        
    end
    
    % NTP consumption models
    if tube.UserData.ReactionConfig.NTPmodel == 1
        % Compute the number of NTPs required, in 100 NTP blocks
        ntpcnt = floor(rna.UserData/100);	% get number of NTP blocks
        if (ntpcnt == 0)
            ntpstr = '';
        else
            ntpstr = int2str(ntpcnt);
        end
        
        txtl_addreaction(tube,['[' RNAPbound '] + ' ntpstr ' NTP <-> [NTP:' RNAPbound ']'],...
            'MassAction',NTPparameters);
    else
        if ~mode.utr_attenuator_flag
            % to deal with stiffness due to high reaction-order
            txtl_addreaction(tube,['[' RNAPbound '] + NTP <-> [NTP:' RNAPbound ']'],...
                'MassAction',NTPparameters);
            
            %dummy raction
            ntpcnt = floor(rna.UserData/100);	% get number of NTP blocks
            NTPConsumptionRate = {'TXTL_NTP_consumption',(ntpcnt-1)*ktx};
            
            txtl_addreaction(tube,['[NTP:' RNAPbound '] -> ' dna.Name ' +  ' RNAP],...
                'MassAction',NTPConsumptionRate);
        else
            % to deal with stiffness due to high reaction-order
            txtl_addreaction(tube,['[' RNAPbound '] + NTP <-> [NTP:' RNAPbound ']'],...
                'MassAction',NTPparameters);
            txtl_addreaction(tube,['[' RNAPbound ':RNA att] + NTP <-> [NTP:' RNAPbound ':RNA att]'],...
                'MassAction',NTPparameters);            
            %dummy raction
            ntpcnt = floor(rna.UserData/100);	% get number of NTP blocks
            NTPConsumptionRate1 = {'TXTL_NTP_consumption1',(ntpcnt-1)*ktx1};
            NTPConsumptionRate2 = {'TXTL_NTP_consumption2',(ntpcnt-1)*ktx2};
            txtl_addreaction(tube,['[NTP:' RNAPbound '] -> ' dna.Name ' +  ' RNAP],...
                'MassAction',NTPConsumptionRate1);
            txtl_addreaction(tube,['[NTP:' RNAPbound ':RNA att] -> ' dna.Name ' +  ' RNAP ' + RNA att'],...
                'MassAction',NTPConsumptionRate2);            
        end
        
    end
    
    % transcription and att-anti reactions
    if mode.utr_attenuator_flag
        txtl_addreaction(tube,transcriptionEq1,'MassAction',{'TXTL_transcription_rate1',ktx1});
        txtl_addreaction(tube,transcriptionEq2,'MassAction',{'TXTL_transcription_rate2',ktx2});
        % Auto-Termination
        auto_termination_rate = {'TXTL_RNA_ATT_AUTOTERM', tube.UserData.ReactionConfig.auto_termination_rate};
        txtl_addreaction(tube,['[NTP:' RNAPbound ':RNA att] -> ' dna.Name ' + [RNA att] + ' RNAP],...
            'MassAction',auto_termination_rate);
        % Complex 1 Formation (prevents extension of att, thereby
        % inhibiting transcription and subsequent trasnlation
        complex1_rate = {'TXTL_RNA_C1_F', tube.UserData.ReactionConfig.complex1_F_rate; ...
            'TXTL_RNA_C1_R', tube.UserData.ReactionConfig.complex1_R_rate};
        txtl_addreaction(tube,['[' RNAPbound ':RNA att] + [RNA anti] <-> [' RNAPbound ':RNA att:RNA anti]'],...
            'MassAction',complex1_rate);         
        % termination
        att_anti_termination_rate = {'TXTL_RNA_ATTANTI_TERM', tube.UserData.ReactionConfig.att_anti_termination_rate};
        txtl_addreaction(tube,['[' RNAPbound ':RNA att:RNA anti] -> ' dna.Name ' + [RNA att:RNA anti] + ' RNAP],...
            'MassAction',att_anti_termination_rate);
    else
        txtl_addreaction(tube,transcriptionEq,'MassAction',{'TXTL_transcription_rate1',ktx});        
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_transcription:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end




% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
% Parameters describing the enzymatic process
