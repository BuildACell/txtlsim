function dna = txtl_adddna(tube, promspec, rbsspec, genespec, amount, type)
%TXTL_ADDDNA   Set up species and reactions for a DNA segment
%
%   dna = TXTL_ADDDNA(tube, promspec, rbsspec, genespec, amount, type)
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
% Modifications: VS, Oct 2012
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
%{
% leave this here for now, in case variable arguments are needed in the
future. 
numvarargs = length(varargin);
if numvarargs < 6
    error('myfuns:txtl_dna:TooFewInputs', ...
        'requires at least (tube, promspec, rbsspec, genespec, amount, type)');
end
[tube, promspec, rbsspec, genespec, amount, type] = varargin{1:6};
if numvarargs == 7
    prepromspec = varargin{7};
else prepromspec = {};
end
%}

% Parameters used in this file
%! TODO: update these parameters to something reasonable
kDNA_recbcd_f = 0.4;	% forward rr for DNA + RecBCD <-> DNA:RecBCD
kDNA_recbcd_r = 0.1;	% backward rr for DNA + RecBCD <-> DNA:RecBCD
kRNA_deg = log(2)/10;			% mRNA degradation: 10 sec half life
% Extract out the names and lengths of the promoter, RBS and gene as cell
% arrays
[promFull, promlen] = txtl_parsespec(promspec);
[rbsFull, rbslen] = txtl_parsespec(rbsspec);
[geneFull, genelen] = txtl_parsespec(genespec);

% forward rr for DNA:RecBCD -> RecBCD
% !TODO: Find out how the number of protection NTPs relate to rr. 
junklength = 0;
thiolength = 0;
junkDNAflag = false;
thioDNAflag = false;
% The format of the DNA is: DNA thio-junk-ptet--rbs-spacer--gene-lav-terminator 
for i = 1:length(promFull)
    if strcmp(promFull{i}, 'junk')
        junkDNAflag = true;
        junklength = promlen{i};
    end
    if strcmp(promFull{i},'thio')
        thioDNAflag = true;
        thiolength = promlen{i};
    end
end

if junkDNAflag
    kDNA_complex_deg = log(2)/(1+junklength/100);	
else
    kDNA_complex_deg = 0.5;
end
% currently, we have a simple protection due to thiosulfate. 
% !TODO: come up with something less arbitrary
if thioDNAflag
    kDNA_complex_deg = 0.5*kDNA_complex_deg;
    % thiolength does not do anything for now. 
end
%
% Create the species for the DNA, RNA and protein
%
% Store the length of the DNA, transcript or protein in userdata.  Used
% to calculate out rate constants.

if length(promFull)>=1
    promstr = promFull{1};
    promlenTot = promlen{1};
    justProm = promFull{end}; % assuming thio-junk-...-prom
end
if length(promFull)>=2
    for i = 2:length(promFull)
        promstr = [promstr '-' promFull{i}]; % expect: thio-junk-prom
        promlenTot = promlenTot+ promlen{i};
    end
    
end
if length(rbsFull)>=1
    rbsstr = rbsFull{1};
    rbslenTot = rbslen{1};
    justRbs = rbsFull{1};
end
if length(rbsFull)>=2
    for i = 2:length(rbsFull)
        rbsstr = [rbsstr '-' rbsFull{i}]; %expect: rbs-spacer
        rbslenTot = rbslenTot+ rbslen{i};
    end
end

% !TODO check out what terminators and lavs are. 
protDEGflag = false;
protTERMflag = false;
for i = 1:length(geneFull)
    if strcmp(geneFull{i}, 'lav')
        protDEGflag = true;
    end
    if strcmp(geneFull{i},'terminator')
        protTERMflag = true;
    end
end
if length(geneFull)>=1
    genestr = geneFull{1}; 
    genelenTot = genelen{1};
    justGene = geneFull{1}; %assuming the format is gene-lav-...-terminator
end
if length(geneFull)>=2
    for i = 2:length(geneFull)
        genestr = [genestr '-' geneFull{i}]; %expect: gene-lav-terminator
        genelenTot = genelenTot+ genelen{i};
    end
end

dnastr = ['DNA ' promstr '--' rbsstr '--' genestr];
dna = addspecies(tube, dnastr, amount);
dna.UserData = promlenTot + rbslenTot + genelenTot;

rnastr = ['RNA ' rbsstr '--' genestr];
rna = addspecies(tube, rnastr);
rna.UserData = rbslenTot + genelenTot;	% length in NTPs

% protstr looks something like 'protein tetR-lav-terminator'
% !Question: Should we change it to just 'protein tetR'?
protstr = ['protein ' genestr]; 
protein = sbioselect(tube, 'Type', 'species', 'Name', protstr);
if isempty(protein)
protein = addspecies(tube, protstr);
end
protein.UserData = genelenTot / 3;		% length in amino acids
% should this be total gene length / 3 or some part of it? Depends on what
% the gene will be made of. apart from degradation tags, what is there?

% Transcription
if exist(['txtl_prom_' justProm]) == 2 
    % instead of using the full thio-junk-prom to name the promoter file, we use just the prom part
    % and pass the dna and rna strings as a arguments. the thio and junk
    % are available in the name of the dna object. 
    
  Rlist = eval(['txtl_prom_' justProm '(tube, dna, rna)']);
else
  % Issue a warning and run the default promoter
  warning(['TXTL: can''t find txtl_prom_' justProm ...
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
if exist(['txtl_utr_' justRbs]) == 2
  % Run the RBS specific setup
  Ribobound = eval(['txtl_utr_' justRbs '(tube, rna, protein)']);
else
  % Issue a warning and run the default RBS
  warning(['TXTL: can''t find txtl_utr_' justRbs ...
      '; using default rbs params']);
  Ribobound = txtl_utr_rbs(tube, rna, protein);
end

% Now put in the reactions for the utilization of amino acids
% Use an enzymatic reaction to proper rate limiting
kf_aa = log(2) / 0.001;			% binding rate of 1 ms
kr_aa = 1 * kf_aa;			% Km of 100 for amino acid usage
% !TODO: check how the adding of degradation tag and termination site affects this.
ktl_rbs = log(2)/(protein.UserData/10);	% 10 AA/second translation 


% Compute the number of amino acids required, in 100 AA blocks
% !TODO: think about some sequential method for doing this. May have to run
% comparative simulations to see if the increased computation in the
% sequential case is offset by a more accurate simulation result. 
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
if exist(['txtl_protein_' justGene]) == 2
  % Run the protein specific setup
  Rlist = eval(['txtl_protein_' justGene '(tube, protein)']);
end
%! TODO: add protein degradation capability
if protDEGflag
  % Run the protein degradation
  degradationRate = 0.0001; % !TODO: Find a reasonable value
  Rlist = txtl_protein_degradation(tube, protein,degradationRate);
end
% All done!
return

% Utility function for parsing out a specification string
function [names, lengths] = txtl_parsespec(spec)
  
  indivRegions = regexp(spec, '-','split'); %cell array of individual xyz(123) strings
  namesAndLengths = regexp(indivRegions, '\w+','match'); %cell array of cells containing the names and lengths of dna regions
  names = cell(1,length(namesAndLengths));
  lengths = cell(1,length(namesAndLengths));
    
  %error checking followed by returning parsed strings
  for i = 1:length(namesAndLengths)
      if isempty(namesAndLengths{i})
          error('txtl_adddna:wrongStringFormat',...
              ['the string %s should be: name(length)-name2(length2)-' ...
              '...-nameN(lengthN), where the lengths are optional. eg: thio-junk(500)-ptet(50)'...
              'the name must start with an alphabet'], spec)
      else
          A = isstrprop(namesAndLengths{i}{1},'alpha');
          if ~A(1) 
              error('txtl_adddna:wrongSpeciesName',...
                  ['species named %s should start with an alphabet. Format is' ...
                  ' name(length). Where the lengths are optional. eg: thio or junk(500)'],indivRegions{i})
          end
      end
      % return the parsed name and optional length
      names{i} = namesAndLengths{i}{1};
      if length(namesAndLengths{i}) == 1
          lengths{i} = 0;
      else if length(namesAndLengths{i})>2
              error('txtl_adddna:tooManyElements',...
                  ['the string %s is not of the format name(length). '...
                  'It has unwanted elements after '')'''],...
                  indivRegions{i});
          else if length(namesAndLengths{i})==2
                  lengths{i} = str2double(namesAndLengths{i}{2});
              end
          end
      end
  end
  
  return

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
