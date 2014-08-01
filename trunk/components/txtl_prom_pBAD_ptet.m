% txtl_prom_pBAD_ptet.m - promoter information for pBAD and ptet combinatorial promoter
% Zoltan Tuza, Oct 2012
% Vipul Singhal Jun 2014
%
% This file contains a description of the pBAD and ptet combinatorial promoter.
% Calling the function txtl_prom_pBAD_ptet() will set up the reactions for
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

function varargout= txtl_prom_pBAD_ptet(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];
    P1 = 'protein sigma70';
    P2 = 'protein tetRdimer';
    P3 = 'protein AraC';
    AraCbound = ['arabinose:' P3];
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('pBAD_tetR');
    
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
    defaultBasePairs = {'pBAD_ptet','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,P1,P2, AraCbound, P3};
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)),'Internal');
    % empty cellarray for amount => zero amount
    if mode.utr_attenuator_flag
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ], prom_spec, rbs_spec, {P2} );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{AraCbound} );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{P2, AraCbound} );
    else
        txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ],{P2}); %lowest rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' AraCbound ],{AraCbound}); %highest rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ':' AraCbound ],{P2, AraCbound}); %slightly higher than 1.
    end

    %(check agains shaobin results. the parameters here should be tuned to
    %get the shaobin curves. translation/degradation etc should be standard.)

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    if nargin==8
    prom_spec = varargin{2};
    rbs_spec = varargin{3};
    gene_spec = varargin{4};
    elseif nargin~=5
        error('the number of argument should be 5 or 8, not %d',nargin);
    end
    
    % parameters for leaky transcription
    parameters = {'TXTL_PBADPTET_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_PBADPTET_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    %


    %% Set up bindig reaction for tetR
    %
    % promoter binds to tetR
    Robj1 = addreaction(tube, [dna.Name ' + ' P2 ' <-> ' dna.Name ':' P2 ]);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    Pobj1f = addparameter(Kobj1, 'kf', 8.86e-1);
    Pobj1r = addparameter(Kobj1, 'kr', 0.11e-4);
    set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});
    
    % the binding of P2 to the DNA-RNAP complex. note that due to the reaction
    % below, this binding means that RNAP will soon leave the DNA.  Hence,
    % P2 binding to RNAPbound expels the PNAP, redusing transcription.
    Robj2 = addreaction(tube, [RNAPbound  ' + ' P2 ' <-> ' RNAPbound ':' P2 ]);
    Kobj2 = addkineticlaw(Robj2, 'MassAction');
    Pobj2f = addparameter(Kobj2, 'kf', 8.86e-1);
    Pobj2r = addparameter(Kobj2, 'kr', 0.11e-4);
    set(Kobj2, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    % Set up binding reaction for tetR. notice that the DNA-RNAP-P2 complex
    % is v unstable, and expels the RNAP readily. 
    Robj3 = addreaction(tube, [dna.Name ':' P2 ' + ' RNAP ' <-> [' RNAPbound ':' P2 ']' ]);
    Kobj3 = addkineticlaw(Robj3, 'MassAction');
    Pobj3f = addparameter(Kobj3, 'kf', paramObj.RNAPbound_Forward);
    Pobj3r = addparameter(Kobj3, 'kr', paramObj.RNAPbound_Reverse*1000);
    set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
    
    %% AraC
    %
    % set up binding reactions for AraC:arabinose. 
    Robj4 = addreaction(tube, [dna.Name ' + ' AraCbound ' <-> ' dna.Name ':' AraCbound ]);
    Kobj4 = addkineticlaw(Robj4, 'MassAction');
    Pobj4f = addparameter(Kobj4, 'kf', 2.86e-1);
    Pobj4r = addparameter(Kobj4, 'kr', 0.11e-4);
    set(Kobj4, 'ParameterVariableNames', {'kf', 'kr'});
    
    % the binding of P2 to the DNA-RNAP complex. note that due to the reaction
    % below, this binding means that RNAP will soon leave the DNA.  Hence,
    % P2 binding to RNAPbound expels the PNAP, redusing transcription.
    Robj5 = addreaction(tube, [RNAPbound  ' + ' AraCbound ' <-> [' RNAPbound ':' AraCbound ']']);
    Kobj5 = addkineticlaw(Robj5, 'MassAction');
    Pobj5f = addparameter(Kobj5, 'kf', 2.86e-2);
    Pobj5r = addparameter(Kobj5, 'kr', 0.11e-4);
    set(Kobj5, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    % Set up binding reaction for tetR. notice that the DNA-RNAP-P2 complex
    % is v unstable, and expels the RNAP readily. 
    Robj6 = addreaction(tube, [dna.Name ':' AraCbound ' + ' RNAP ' <-> [' RNAPbound ':' AraCbound ']' ]);
    Kobj6 = addkineticlaw(Robj6, 'MassAction');
    Pobj6f = addparameter(Kobj6, 'kf', paramObj.RNAPbound_Forward*300);
    Pobj6r = addparameter(Kobj6, 'kr', paramObj.RNAPbound_Reverse);
    set(Kobj6, 'ParameterVariableNames', {'kf', 'kr'});
    
    %% AraC and tetR
   
    % set up binding reactions for AraC:arabinose. 
    Robj7 = addreaction(tube, [dna.Name ':' AraCbound ' + ' P2 ' <-> ' dna.Name ':' P2 ':' AraCbound ]);
    Kobj7 = addkineticlaw(Robj7, 'MassAction');
    Pobj7f = addparameter(Kobj7, 'kf', 2.86e-3);
    Pobj7r = addparameter(Kobj7, 'kr', 0.11e-4);
    set(Kobj7, 'ParameterVariableNames', {'kf', 'kr'});
    
    Robj8 = addreaction(tube, [dna.Name ':' P2 ' + ' AraCbound ' <-> ' dna.Name ':' P2 ':' AraCbound ]);
    Kobj8 = addkineticlaw(Robj8, 'MassAction');
    Pobj8f = addparameter(Kobj8, 'kf', 2.86e-1);
    Pobj8r = addparameter(Kobj8, 'kr', 0.11e-4);
    set(Kobj8, 'ParameterVariableNames', {'kf', 'kr'});
    
    % the binding of P2 to the DNA-RNAP complex. note that due to the reaction
    % below, this binding means that RNAP will soon leave the DNA.  Hence,
    % P2 binding to RNAPbound expels the PNAP, redusing transcription.
Robj9 = addreaction(tube, [RNAPbound ':' AraCbound  ' + ' P2 ' <-> [' RNAPbound ':' P2 ':' AraCbound ']']);
    Kobj9 = addkineticlaw(Robj9, 'MassAction');
    Pobj9f = addparameter(Kobj9, 'kf', 2.86e-3);
    Pobj9r = addparameter(Kobj9, 'kr', 0.11e-4);
    set(Kobj9, 'ParameterVariableNames', {'kf', 'kr'});
    
        Robj10 = addreaction(tube, [RNAPbound ':' P2  ' + ' AraCbound ' <-> [' RNAPbound ':' P2 ':' AraCbound ']']);
    Kobj10 = addkineticlaw(Robj10, 'MassAction');
    Pobj10f = addparameter(Kobj10, 'kf', 2.86e-1);
    Pobj10r = addparameter(Kobj10, 'kr', 0.11e-4);
    set(Kobj10, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    % Set up binding reaction for tetR. notice that the DNA-RNAP-P2 complex
    % is v unstable, and expels the RNAP readily. 
    Robj11 = addreaction(tube, [ dna.Name ':' P2 ':' AraCbound ' + ' RNAP ' <-> [' RNAPbound ':' P2 ':' AraCbound ']' ]);
    Kobj11 = addkineticlaw(Robj11, 'MassAction');
    Pobj11f = addparameter(Kobj11, 'kf', paramObj.RNAPbound_Forward);
    Pobj11r = addparameter(Kobj11, 'kr', paramObj.RNAPbound_Reverse);
    set(Kobj11, 'ParameterVariableNames', {'kf', 'kr'});    
    
    %%

    if mode.utr_attenuator_flag
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ], prom_spec, rbs_spec, {P2} );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{AraCbound} );
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{P2, AraCbound} );
    else
        txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ],{P2}); %lowest rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' AraCbound ],{AraCbound}); %highest rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ':' AraCbound ],{P2, AraCbound}); %slightly higher than 1.
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_pBAD_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
