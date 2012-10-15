% txtl_default_lengths.m - template file for MATLAB functions
% VS, Oct 2012
%
% This file provides the default lengths in # of NTPs of the various DNA
% fragments if none are provided by the user. It returns a cell array of 

% Written by Vipul Singhal 7 OCT 2012
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

function [promFull, promlen, rbsFull, rbslen, geneFull, genelen] = txtl_default_lengths(promFull, promlen, rbsFull, rbslen, geneFull, genelen)

% I think this error never occurs, can remove when this is confirmed
if ~isequal(length(promFull), length(promlen)) || ~isequal(length(rbsFull), length(rbslen)) || ~isequal(length(geneFull), length(genelen))
    error('number of promoter fragments must match the number of promoter lengths');
end

promDefaultUsed = 0;
rbsDefaultUsed = 0;
geneDefaultUsed = 0;

% set up promoter default lengths
for i = 1: length(promFull)
    if promlen{i} == 0
        promDefaultUsed = promDefaultUsed+1;
        promDefIdx(promDefaultUsed) = i;
    end
end
if promDefaultUsed ~= 0
    for i = 1:length(promDefIdx)
        switch promFull{promDefIdx(i)}
            case 'ptet'
                promlen{promDefIdx(i)} = 50;
            case 'p70'
                promlen{promDefIdx(i)} = 50;
            case 'ptrc2'
                promlen{promDefIdx(i)} = 50;
            case 'placI'
                promlen{promDefIdx(i)} = 50;
            case 'junk'
                promlen{promDefIdx(i)} = 500;           
        end
    end
end

% defaults for the UTR
rbsDefIdx = [];
for i = 1: length(rbsFull)
    if rbslen{i} == 0
        rbsDefaultUsed = rbsDefaultUsed+1
        rbsDefIdx(rbsDefaultUsed) = i
    end
end

if rbsDefaultUsed ~= 0
    for i = 1:length(rbsDefIdx)
        switch rbsFull{rbsDefIdx(i)}
            case 'rbs'
                rbslen{rbsDefIdx(i)} = 20;
            case 'spacer'
                rbslen{rbsDefIdx(i)} = 200;
        end
    end
end

% gene default lengths
geneDefIdx = 0;
for i = 1: length(geneFull)
    if genelen{i} == 0
        geneDefaultUsed = geneDefaultUsed+1;
        geneDefIdx(geneDefaultUsed) = i;
    end
end
if geneDefaultUsed ~= 0
    for i = 1:length(geneDefIdx)
        switch geneFull{geneDefIdx(i)}
            case 'tetR'
                genelen{geneDefIdx(i)} = 647;
            case 'lacI'
                genelen{geneDefIdx(i)} = 647;
            case 'deGFP'
                genelen{geneDefIdx(i)} = 1000;
            case 'gamS'
                genelen{geneDefIdx(i)} = 1000;
            % any others
        end
    end        
end



            
            



end


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
