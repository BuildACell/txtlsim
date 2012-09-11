function dna = txtl_dna(tube, promspec, rbsspec, genespec, amount, type)
%TXTL_DNA   Set up species and reactions for a DNA segment
%
%   dna = TXTL_DNA(tube, promspec, rbsspec, genespec, amount, type)
%   constructs the species and reactions required for transcription,
%   translation and degradation of DNA, mRNA and proteins in the 
%   TX-TL system.
%
%   * tube = Simbiology model object
%   * promspec = spec of the form 'prom(nn)' where 'prom' is the 
%     promoter name and 'len' is the length of the promoter.
%   * rbsspec = spec of the form 'rbs(nn)' where 'rbs' is the RBS 
%     name and 'len' is the length of the RBS.
%   * genespec = spec of the form 'gene(nn)' where 'gene' is the 
%     gene name and 'len' is the length of the gene.
%   * amount = amount of DNA to put in the tube (in nM)
%   * type = 'linear' if you want to include degradation reactions

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

% Parameters used in this file
%! TODO: update these parameters to something reasonable
kDNA_deg = 0.1;				% Binding rate of RecBCD
kRNA_deg = log(2)/10;			% mRNA degradation: 10 sec half life

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
%! TODO: allow protection strands by comparing promoter length to ~35
if strcmp(type, 'linear')
  Robj1 = addreaction(tube, ...
    [dna.Name ' + RecBCD -> RecBCD']);
  Kobj1 = addkineticlaw(Robj1,'MassAction');
  %! TODO: update to include DNA protection
  Pobj1 = addparameter(Kobj1, 'kf', kDNA_deg);
  set(Kobj1, 'ParameterVariableNames', {'kf'});
end

% Translation: setup file should return pointer to RBS bound species
if exist(['txtl_utr_' rbs]) == 2
  % Run the RBS specific setup
  Ribobound = eval(['txtl_utr_' rbs '(tube, rna, protein)']);
else
  % Issue a warning and run the default RBS
  warning(['TXTL: can''t find txtl_utr_' rbs ...
      '; using default promoter params']);
  Ribobound = txtl_utr_rbs(tube, rna, protein);
end

% Now put in the reactions for the utilization of amino acids
% Use an enzymatic reaction to proper rate limiting
kf_aa = log(2) / 0.001;			% binding rate of 1 ms
kr_aa = 1 * kf_aa;			% Km of 100 for amino acid usage
ktl_rbs = log(2)/(protein.UserData/10);	% 10 AA/second translation

% Compute the number of amino acids required, in 100 AA blocks
aacnt = floor(protein.UserData/100);	% get number of K amino acids
if (aacnt == 0) 
  aastr = '';
else
  aastr = int2str(aacnt);
end

% Set up the translation reaction
Robj = addreaction(tube, ...
  ['[' Ribobound.Name '] + ' aastr ' AA <-> [AA:' Ribobound.Name ']']);
Kobj = addkineticlaw(Robj, 'MassAction');
Pobjf = addparameter(Kobj, 'kf', kf_aa);
Pobjr = addparameter(Kobj, 'kr', kr_aa);
set(Kobj, 'ParameterVariableNames', {'kf', 'kr'});

Robj1 = addreaction(tube, ...
  ['[AA:' Ribobound.Name '] -> ' rna.Name ' + ' protein.Name ' +  Ribo']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1 = addparameter(Kobj1, 'ktl', ktl_rbs);
set(Kobj1, 'ParameterVariableNames', {'ktl'});

% Add in mRNA degradation reactions
Robj2 = addreaction(tube, [rna.Name ' + RNase -> RNase']);
Kobj2 = addkineticlaw(Robj2,'MassAction');
Pobj2 = addparameter(Kobj2, 'kf', kRNA_deg);
set(Kobj2, 'ParameterVariableNames', {'kf'});

% Protein reactions + degradation (if tagged)
if exist(['txtl_protein_' gene]) == 2
  % Run the protein specific setup
  Rlist = eval(['txtl_protein_' gene '(tube, protein)']);
end
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

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
