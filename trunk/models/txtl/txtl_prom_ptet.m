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

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'Setup Species')
    
    promoterData = varargin{1};
    defaultBasePairs = {'ptet','junk','thio';50,500,0};
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
    %! TODO: replace these values with correct values
    
    %generating unique parameter name for the current RNA

    kf_ptet = log(2)/0.1;			% 100 ms bind rate
    kr_ptet = 10 * kf_ptet;			% Km of 10 (same as p70, from VN)
    %ktx_ptet = log(2)/(rna.UserData/30);	% 30 base/second transcription

    % Set up binding reaction
    Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    Pobj1f = addparameter(Kobj1, 'TXTL_PTET_RNAPbound_F', kf_ptet);
    Pobj1r = addparameter(Kobj1, 'TXTL_PTET_RNAPbound_R', kr_ptet);
    set(Kobj1, 'ParameterVariableNames', {'TXTL_PTET_RNAPbound_F', 'TXTL_PTET_RNAPbound_R'});

    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
    %

    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);

    %
    % Add reactions for sequestration of promoter by tetRdimer 
    %

    %! TODO: RMM, 29 Sep 2012
    %! TODO: txtl_protein_tetR defines tetramers, which aren't used
    %! TODO: proper implementation for tetR is via two operator sites (I think)
    % VS: yes, there are 2 operators, see for example, 
    % C Berens and W. Hillen, Gene regulation by tetracyclines, Eur. J. Biochem. 270, 3109�3121 (2003)
    %{
    kf1_tetR = 0.2; kr1_tetR = 1;		% reaction rates (from sbio)
    Robj4 = addreaction(tube, ...
      [DNA ' + [protein tetRdimer] <-> [' dna.name ':protein tetRdimer]']);
    Kobj4 = addkineticlaw(Robj4,'MassAction');
    Pobj4 = addparameter(Kobj4, 'k4', kf1_tetR);
    Pobj4r = addparameter(Kobj4, 'k4r', kr1_tetR);
    set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

    kf2_tetR = 0.2; kr2_tetR = 1;		
    Robj5 = addreaction(tube, ...
      ['[' dna.name ':protein tetRdimer] + [protein tetRdimer] <-> [' dna.name ':protein tetRdimer:protein tetRdimer]']);
    Kobj5 = addkineticlaw(Robj5,'MassAction');
    Pobj5 = addparameter(Kobj5, 'k5', kf2_tetR);
    Pobj5r = addparameter(Kobj5, 'k5r', kr2_tetR);
    set(Kobj5, 'ParameterVariableNames', {'k5', 'k5r'});

    kf3_tetR = 0.2; kr3_tetR = 1;		% reaction rates (from sbio)
    Robj7 = addreaction(tube, ...
      [DNA ' + [protein tetR-lvadimer] <-> [' dna.name ':protein tetR-lvadimer]']);
    Kobj7 = addkineticlaw(Robj7,'MassAction');
    Pobj7 = addparameter(Kobj7, 'k7', kf3_tetR);
    Pobj7r = addparameter(Kobj7, 'k7r', kr3_tetR);
    set(Kobj7, 'ParameterVariableNames', {'k7', 'k7r'});

    kf4_tetR = 0.2; kr4_tetR = 1;		
    Robj6 = addreaction(tube, ...
      ['[' dna.name ':protein tetR-lvadimer] + [protein tetR-lvadimer] <-> [' dna.name ':protein tetR-lvadimer:protein tetR-lvadimer]']);
    Kobj6 = addkineticlaw(Robj6,'MassAction');
    Pobj6 = addparameter(Kobj6, 'k6', kf4_tetR);
    Pobj6r = addparameter(Kobj6, 'k6r', kr4_tetR);
    set(Kobj6, 'ParameterVariableNames', {'k6', 'k6r'});
    %}
    
    ptetRepression = false;
    %! TODO make all these reactions conditional on specie availability
    matchStr = regexp(listOfSpecies,'(^protein tetR.*dimer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
    listOftetRdimer = vertcat(matchStr{:});
    if ~isempty(listOftetRdimer)
        ptetRepression = true;
    end
    
    if ptetRepression
        for i = 1:size(listOftetRdimer,1)
            Robj8 = addreaction(tube, ...
              [DNA ' + ' listOftetRdimer{i} ' <-> [' dna.name ':' listOftetRdimer{i} '1]']);
            Kobj8 = addkineticlaw(Robj8,'MassAction');
            rN = regexprep(listOftetRdimer{i}, {'( )'}, {''});
            uniqueNameF = sprintf('TXTL_PTET_REPRESSION1_%s_F',rN);
            uniqueNameR = sprintf('TXTL_PTET_REPRESSION1_%s_R',rN);
            set(Kobj8, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});
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
    
        

    % We comment out the second dimer binding for now. We can add these
    % affects at a later date. 
    %{
    kf7_tetR = 0.0000006; kr7_tetR = 0.0000005;		% effectively nonexistent
    Robj9 = addreaction(tube, ...
      ['[' dna.name ':protein tetR-lva-terminatordimer1] + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer]']);
    Kobj9 = addkineticlaw(Robj9,'MassAction');
    Pobj9 = addparameter(Kobj9, 'k9', kf7_tetR);
    Pobj9r = addparameter(Kobj9, 'k9r', kr7_tetR);
    set(Kobj9, 'ParameterVariableNames', {'k9', 'k9r'});

    kf8_tetR = 0.0000006; kr8_tetR = 0.0000005;
    Robj11 = addreaction(tube, ...
      ['[' dna.name ':protein tetR-lva-terminatordimer2] + [protein tetR-lva-terminatordimer] <-> [' dna.name ':protein tetR-lva-terminatordimer:protein tetR-lva-terminatordimer]']);
    Kobj11 = addkineticlaw(Robj11,'MassAction');
    Pobj11 = addparameter(Kobj11, 'k11', kf8_tetR);
    Pobj11r = addparameter(Kobj11, 'k11r', kr8_tetR);
    set(Kobj11, 'ParameterVariableNames', {'k11', 'k11r'});
    %}
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
