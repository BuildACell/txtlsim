%TXTL_DNA   Set up species and reactions for a DNA segment
%
% dna = txtl_dna(tube, promspec, rbsspec, genespec, amount, type)
%
% tube = Simbiology model object
% promspec = spec of the form 'prom(nn)' where 'prom' is the promoter name
%   and 'len' is the length of the promoter.
% rbsspec = spec of the form 'rbs(nn)' where 'rbs' is the RBS name
%   and 'len' is the length of the RBS.
% genespec = spec of the form 'gene(nn)' where 'gene' is the gene name
%   and 'len' is the length of the gene.
% amount = amount of DNA to put in the tube (in nM)
% type = 'linear' if you want to include degradation reactions
% 

% Written by Richard Murray, Sep 2012
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

function dna = txtl_dna(tube, promspec, rbsspec, genespec, amount, type)

% Parameters used in this file
kDNA_deg = 0.1;				% Binding rate of RecBCD

% Extract out the names and lengths of the promoter, RBS and gene
[prom, promlen] = txtl_parsespec(promspec);
[rbs, rbslen] = txtl_parsespec(rbsspec);
[gene, genelen] = txtl_parsespec(genespec);

%
% Create the species for the DNA, RNA and protein
%
% Store the length of the DNA, transcript or protein in userdata.  Used
% to calculate out rate constants.

dnastr = ['DNA ' prom '=' rbs '=' gene];
dna = addspecies(tube, dnastr, amount);
dna.UserData = promlen + rbslen + genelen;

rnastr = ['RNA ' rbs '=' gene];
rna = addspecies(tube, rnastr);
rna.UserData = rbslen + genelen;	% length in NTPs

protstr = ['protein ' gene];
protein = addspecies(tube, protstr);
protein.UserData = genelen / 3;		% length in amino acids

% Transcription
if exist(['txtl_prom_' prom]) == 2
  % Run the promoter specific setup
  Rlist = eval(['txtl_prom_' prom '(tube, dna, rna)']);
else
  % Issue a warning and run the default promoter
  warning(['TXTL: can''t find txtl_prom_' prom ...
      '; using default promoter params']);
  Rlist = txtl_prom_p70(tube, dna, rna);
end

% DNA degradation
%! TODO: eventually, we should allow file-based reactions
if strcmp(type, 'linear')
  Robj = addreaction(tube, ...
    [dna.Name ' + RecBCD -> RecBCD']);
  Kobj = addkineticlaw(Robj,'MassAction');
  %! TODO: update to include DNA protection
  Pobj = addparameter(Kobj, 'kf', kDNA_deg);
  set(Kobj, 'ParameterVariableNames', {'kf'});
end

% RNA + translation + degradation
if exist(['txtl_rbs_' rbs]) == 2
  % Run the RBS specific setup
  Rlist = eval(['txtl_rbs_' rbs '(tube, rna, protein)']);
else
  % Issue a warning and run the default promoter
  warning(['TXTL: can''t find txtl_rbs_' prom ...
      '; using default promoter params']);
  txtl_rbs_rbs(tube, rna, protein)
end

% Protein + degradation (if tagged)
%! TODO: add protein degradation capability

% All done!
return

% Utility function for parsing out a specification string
function [name, len] = txtl_parsespec(spec)
  %! TODO: use 'split' instead?
  tokens = regexp(spec, '(\w*)\((\d*)\)', 'tokens');
  %! TODO: add error checking
  name = tokens{1}{1};
  len = str2num(tokens{1}{2});
  return
