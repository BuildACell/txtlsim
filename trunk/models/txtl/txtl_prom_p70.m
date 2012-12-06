% txtl_prom_p70.m - promoter information for p70 promoter
% RMM, 8 Sep 2012
%
% This file contains a description of the standard p70 promoter.
% Calling the function txtl_prom_p70() will set up the reactions for
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

function varargout = txtl_prom_p70(mode, tube, dna, rna, varargin)

if strcmp(mode, 'Setup Species')
    
    promFull = varargin{1};
    promlen = varargin{2};


    % set up promoter default lengths
    promDefaultUsed = 0;
    for i = 1: length(promFull)
        if isempty(promlen{i})
            promDefaultUsed = promDefaultUsed+1;
            promDefIdx(promDefaultUsed) = i; %idx of segments to set defaults for
        end
    end

    if promDefaultUsed ~= 0
        for i = 1:length(promDefIdx)
            switch promFull{promDefIdx(i)}
                case 'p70'
                    promlen{promDefIdx(i)} = 50;
                case 'junk'
                    promlen{promDefIdx(i)} = 500; 
                case 'thio'
                    promlen{promDefIdx(i)} = 0; 
            end
        end
    end

    % Create strings for reactants and products
    RNAPbound = ['RNAP70:' dna.Name];	% Name of bound complex
    RNAP = 'RNAP70';
    foo = sbioselect(tube, 'Name', 'RNAP70');
    if isempty(foo)
        addspecies(tube, 'RNAP70');
    end
    foo = [];
    foo = sbioselect(tube, 'Name', RNAPbound);
    if isempty(foo)
        addspecies(tube, RNAPbound);
    end
    foo = [];
    varargout{1} = promlen;
    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
    % 
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);

    
elseif strcmp(mode, 'Setup Reactions')
    
    listOfSpecies = varargin{1};
    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];	% Name of bound complex

    %
    % Set up binding reaction
    %
    Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> ' RNAPbound]);
    Kobj1 = addkineticlaw(Robj1, 'MassAction');
    set(Kobj1, 'ParameterVariableNames', {'TXTL_P70_RNAPbound_F', 'TXTL_P70_RNAPbound_R'});
    %
    % Now put in the reactions for the utilization of NTPs
    % Use an enzymatic reaction to proper rate limiting
    % 
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);

    
else 
        error('txtltoolbox:txtl_prom_p70:undefinedmode', 'The possible modes are ''Setup Species'' and ''Setup Reactions''.')
end 

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
