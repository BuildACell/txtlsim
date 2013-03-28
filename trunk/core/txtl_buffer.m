% txtl_extract.m - function to create a tube of TX-TL buffer
%! TODO: add documentation

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

function tube = txtl_buffer(name)
tube = txtl_newtube(name);
TXTLconfig = txtl_reaction_config(name);

% Add in NTPs and amino acids
% Buffer contents
%    NTP: ATP 1.2mM, GTP 1.2mM, CTP 1.2mM, UTP 1.2mM
%    AA: 1.5nM (for each amino acid)
% due to limiting factors only 20% percent can be utilized, c.f. V Noireaux
% 2003.
%
% Buffer is 41% of the 10ul reaction volume
stockMulti = 10/4.1666; 
 addspecies(tube, 'NTP', stockMulti*TXTLconfig.NTP_Concentration);		% 100 ntp's/unit (NTP = ATP+GTP+CTP+UTP)
 addspecies(tube, 'AA',  stockMulti*TXTLconfig.AA_Concentration);		% 100 aa's/unit (20 amino acid)



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
