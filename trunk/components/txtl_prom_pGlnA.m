% txtl_prom_pGlnA.m - promoter information for glnA promoter file
% Zoltan Tuza, Sep 2013
%
% 
% 

% Written by Zoltan Tuza, Sep 2013
%
% Copyright (c) 2013 by California Institute of Technology
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

function varargout= txtl_prom_pGlnA(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP54';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP54:' dna.Name];
    P1 = 'protein sigma54';
    
    P3 = 'protein NRI-p';
    
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('pGlnA');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')

    promoterData = varargin{1};
    defaultBasePairs = {'pGlnA','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,P1, P3};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
 
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3}); 
     
    %(check agains shaobin results. the parameters here should be tuned to
    %get the shaobin curves. translation/degradation etc should be standard. also, the params modifies are the 
    %RNAP binding affinities, given below. so, nothing in tx: if RNAP is bound, tx proceeds as normal 
    %(what effects will this have for future work?))

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    % Parameters that describe this promoter (this is where the variation
    % in the promoter strength comes in. 
    parameters = {'TXTL_pGlnA_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_pGlnA_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    %
    % basal transcription
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);


    %% NRI-p 
    %
    % set up binding reactions for AraC:arabinose. 
    txtl_addreaction(tube,[dna.Name ' + ' P3 ' <-> ' dna.Name ':' P3 ],...
     'MassAction',{'TXTL_DNA_NRI-p_F',2.86e-3;'TXTL_DNA_NRI-p_R',0.11e-4});
    
    
    % the binding of P2 to the DNA-RNAP complex. note that due to the reaction
    % below, this binding means that RNAP will soon leave the DNA.  Hence,
    % P2 binding to RNAPbound expels the RNAP, redusing transcription.
     txtl_addreaction(tube,[RNAPbound  ' + ' P3 ' <-> [' RNAPbound ':' P3 ']'],...
     'MassAction',{'TXTL_RNAPbound_NRI-p_F',2.86e-3;'TXTL_RNAPbound_NRI-p_R',0.11e-4});
    
    % 
    % Set up binding reaction for tetR. notice that the DNA-RNAP-P2 complex
    % is v unstable, and expels the RNAP readily. 
     txtl_addreaction(tube,[dna.Name ':' P3 ' + ' RNAP ' <-> [' RNAPbound ':' P3 ']' ],...
     'MassAction',{'TXTL_DNA_NRI-p_RNAPbound_F',paramObj.RNAPbound_Forward*50;'TXTL_DNA_NRI-p_RNAPbound_R',paramObj.RNAPbound_Reverse});
    

    
        
    
    %%
   
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
    
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3}); %highest rate
   

    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_pGlnA:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
