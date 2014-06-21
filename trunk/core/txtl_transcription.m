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
% if called with extra species, then RNAPbound is more than just RNAP70:DNA. there are other sopecies also attached.

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    RNAPbound_term = ['term_' RNAPbound];
    if nargin < 6
        error('the number of argument should be at least 6, not %d',nargin);
    elseif nargin > 6
        extraSpecies = varargin{6};
        coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term,RNAP,extraSpecies{:}};
    else
        coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term,RNAP};
    end
    
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    ktxExpression =  strrep(tube.Userdata.ReactionConfig.Transcription_Rate,...
        'RNA_Length','rna.UserData');
    ktx = eval(ktxExpression); %kt/rna_length = 1.5(ntps^-1) / rnalength(ntp)
    
    ntpcnt = round(rna.UserData/2); 
    %is this divided by 4 or is it to be divided by 2? Since the AGTP conc = ATP conc + GTP
    %conc, so that at every step, when AGTP:CUTP gets used, actually only 2
    %nucleotides are bing used. Maybe the solution is to set AGTP conc =
    %ATP conc, and CUTP conc = CTP conc? think about this. Or actually I
    %will just have ntpcnt = rna.length/2 above. That way we can keep AGTP
    %= atp + gtp. 
    NTPConsumptionRate = {'TXTL_NTP_consumption',(ntpcnt-1)*ktx};
    
    RNAPbound_term = ['term_' RNAPbound];
    transcriptionEq = ...
        ['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound_term ' + ' rna.Name];
    if nargin < 6
        error('the number of argument should be at least 6, not %d',nargin);
    elseif nargin > 6
        extraSpecies = varargin{6};
        % processing the extraSpecies
        extraStr = extraSpecies{1};
        for k=2:size(extraSpecies,2)
            extraStr = [extraStr ' + ' extraSpecies{k}];
        end
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound '] -> ' RNAP  ' + ' dna.Name ' + ' extraStr],...
        'MassAction',NTPConsumptionRate); 

        txtl_addreaction(tube,['[' RNAPbound_term '] -> ' RNAP  ' + ' dna.Name ' + ' extraStr],...
            'MassAction',{'TXTL_RNAPBOUND_TERMINATION_RATE', tube.UserData.ReactionConfig.RNAPbound_termination_rate});
        
    else
        txtl_addreaction(tube,['[CUTP:AGTP:' RNAPbound '] -> ' RNAP  ' + ' dna.Name],...
        'MassAction',NTPConsumptionRate);
    %notice that this is still the separation method. in the next release we will do 
    %['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound], and recharacterize
    %params. 
        txtl_addreaction(tube,['[' RNAPbound_term '] -> ' RNAP  ' + ' dna.Name],...
            'MassAction',{'TXTL_RNAPBOUND_TERMINATION_RATE', tube.UserData.ReactionConfig.RNAPbound_termination_rate});
    end
    
    
    
    
    
    % parameter values
    NTPparameters = {'TXTL_NTP_RNAP_F', tube.UserData.ReactionConfig.NTP_Forward;
        'TXTL_NTP_RNAP_R', tube.UserData.ReactionConfig.NTP_Reverse};
    NTPparameters_fast = {'TXTL_NTP_RNAP_F', 1000*tube.UserData.ReactionConfig.NTP_Forward;
        'TXTL_NTP_RNAP_R', 1000*tube.UserData.ReactionConfig.NTP_Reverse};
    
    % bind nucleotides
    txtl_addreaction(tube,['[' RNAPbound '] + AGTP <-> [AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_fast);
    txtl_addreaction(tube,['[' RNAPbound '] + CUTP <-> [CUTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_fast);
    txtl_addreaction(tube,['[AGTP:' RNAPbound '] + CUTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);
    txtl_addreaction(tube,['[CUTP:' RNAPbound '] + AGTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters);
    
    
    
    
    txtl_addreaction(tube,transcriptionEq,'MassAction',{'TXTL_transcription_rate1',ktx});
    
    
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
