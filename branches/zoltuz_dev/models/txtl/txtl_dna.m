function dna = txtl_dna(varargin)
%TXTL_DNA   Set up species and reactions for a DNA segment
%
%   dna = TXTL_DNA(tube, promspec, rbsspec, genespec, amount, type)
%   constructs the species and reactions required for transcription,
%   translation and degradation of DNA, mRNA and proteins in the 
%   TX-TL system.
%
%   * tube = Simbiology model object
%   * prepromspec = Cell array of nucleatide sequences and corresponding
%   sizes. One example of their use is as a protection from exonucleases. 
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
numvarargs = length(varargin);
if numvarargs < 6
    error('myfuns:txtl_dna:TooFewInputs', ...
        'requires at least (tube, promspec, rbsspec, genespec, amount, type)');
end
[tube, promspec, rbsspec, genespec, amount, type] = varargin{1:6};
if numvarargs == 7
    prepromspec = varargin{7};
else prepromspec = {}
end



% Parameters used in this file
%! TODO: update these parameters to something reasonable
kDNA_recbcd_f = 0.4;	% forward rr for DNA + RecBCD <-> DNA:RecBCD
kDNA_recbcd_r = 0.1;	% backward rr for DNA + RecBCD <-> DNA:RecBCD
kRNA_deg = log(2)/10;			% mRNA degradation: 10 sec half life

% there can be multiple nucleotide sequences before a promoter, and if any
% of them is a protection secuence, the last one (if there are multiple,
% there shouldnt be!) will enable a reduction in DNA degradation rate. 
junklength = 0;
junkDNAflag = false;
thioDNAflag = false;
junk = '';
thio = '';
if ~isempty(prepromspec)
    preprom = cell(1,length(prepromspec));
    prepromlen =  cell(1,length(prepromspec));
    for i = 1:length(prepromspec)
        [preprom{i}, prepromlen{i}] = txtl_parsespec(prepromspec{i});
        if preprom{i} == 'junk'
            junkDNAflag = true;
            iJunk = i;
            junklength = prepromlen{i};
            junk = 'junk-';
        end
        if preprom{i} == 'thio'
            iThio = i;
            thioDNAflag = true;
            thio = 'thio-';
            %thiolength = prepromlen{i};
        end
    end
end
% forward rr for DNA:RecBCD -> RecBCD
% !TODO: Read BFS to learn how the number of protection NTPs relate to rr. 
% For now, use this simple expression. 
if junkDNAflag
    kDNA_complex_deg = log(2)/(1+junklength/100);	
else
    kDNA_complex_deg = 0.5;
end
% currently, we have a simple protection due to thiosulfate. 
% !TODO: come up with something less arbitrary
if thioDNAflag
    kDNA_complex_deg = 0.5*kDNA_complex_deg;
end


  

% Extract out the names and lengths of the promoter, RBS and gene
[prom, promlen] = txtl_parsespec(promspec);
[rbs, rbslen] = txtl_parsespec(rbsspec);
[gene, genelen] = txtl_parsespec(genespec);

%
% Create the species for the DNA, RNA and protein
%
% Store the length of the DNA, transcript or protein in userdata.  Used
% to calculate out rate constants.

dnastr = ['DNA ' thio junk prom '=' rbs '=' gene];
dna = addspecies(tube, dnastr, amount);
dna.UserData = junklength + promlen + rbslen + genelen;

rnastr = ['RNA ' rbs '=' gene];
rna = addspecies(tube, rnastr);
rna.UserData = rbslen + genelen;	% length in NTPs

protstr = ['protein ' gene];
protein = sbioselect(tube, 'Type', 'species', 'Name', protstr);
if isempty(protein)
protein = addspecies(tube, protstr);
end
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
   Rlist = txtl_dna_degradation(tube, dna, [kDNA_recbcd_f, kDNA_recbcd_r, kDNA_complex_deg]); 
   % get reaction rates accordingsly, from user data or extraccted from the properties. 
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
  if isempty(tokens{1}{1}) | isempty(tokens{1}{2})  
     error('wrong string format: %s, it should be: name(length), e.g. p70(50)',spec) 
  end
  name = tokens{1}{1};
  len = str2num(tokens{1}{2});
  return

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
