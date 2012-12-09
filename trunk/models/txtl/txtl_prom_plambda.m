% txtl_prom_plambda.m - promoter information for plambda promoter
% VS Dec 2012
%
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.

% Written by Vipul Singhal, Dec 2012
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

function varargout = txtl_prom_plambda(mode, tube, dna, rna,varargin)

    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')
    
    
    promoterData = [varargin{1};varargin{2}];
    defaultBasePairs = {'plambda','junk','thio';50,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData(2,:);

    coreSpecies = {RNAP,RNAPbound};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)));
    
    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
    %
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'Setup Reactions')
    listOfSpecies = varargin{1};
    
    % Parameters that describe this promoter
    %! TODO: replace these values with correct values
    
    %generating unique parameter name for the current RNA

    kf_plambda = log(2)/0.1;			% 100 ms bind rate
    kr_plambda = 10 * kf_plambda;			% Km of 10 (same as p70, from VN)
    %ktx_ptet = log(2)/(rna.UserData/30);	% 30 base/second transcription

    

    % Set up binding reaction
    Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    Pobj1f = addparameter(Kobj1, 'TXTL_PLAMBDA_RNAPbound_F', kf_plambda);
    Pobj1r = addparameter(Kobj1, 'TXTL_PLAMBDA_RNAPbound_R', kr_plambda);
    set(Kobj1, 'ParameterVariableNames', {'TXTL_PLAMBDA_RNAPbound_F', 'TXTL_PLAMBDA_RNAPbound_R'});

    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
    %

    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound); 
    
    plambdaRepression = false;
    %! TODO make all these reactions conditional on specie availability
    matchStr = regexp(listOfSpecies,'(^protein lambda.*dimer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
    listOflambdadimer = vertcat(matchStr{:});
    if ~isempty(listOflambdadimer)
        plambdaRepression = true;
    end
    
    if plambdaRepression
        for i = 1:size(listOflambdadimer,1)
            Robj8 = addreaction(tube, ...
              [DNA ' + ' listOflambdadimer{i} ' <-> [' dna.name ':' listOflambdadimer{i} ']']);
            Kobj8 = addkineticlaw(Robj8,'MassAction');
            rN = regexprep(listOflambdadimer{i}, {'( )'}, {''});
            uniqueNameF = sprintf('TXTL_PLAMBDA_REPRESSION_%s_F',rN);
            uniqueNameR = sprintf('TXTL_PLAMBDA_REPRESSION_%s_R',rN);
            set(Kobj8, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});            
        end
    end 
        
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_plambda:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
