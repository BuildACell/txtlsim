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

% Parameters used in this file
kDNA_recbcd_f = 0.4;% forward rr for DNA + RecBCD <-> DNA:RecBCD
kDNA_recbcd_r = 0.1;% backward rr for DNA + RecBCD <-> DNA:RecBCD
% Extract out the names and lengths
[promFull, promlen] = txtl_parsespec(promspec);
[rbsFull, rbslen] = txtl_parsespec(rbsspec);
[geneFull, genelen] = txtl_parsespec(genespec);

% set up protein reactions and data, followed by utr followed by promoter
% (promoter reactions require the lengths of the rna, therefore need to be
% set up after the protein and utr files are called.

%% Protein properties, parameters and reactions 
protDEGflag = false; % check for degradation tag and terminator
protTERMflag = false;
for i = 1:length(geneFull)
    if strcmp(geneFull{i}, 'lva')
        protDEGflag = true;
    end
    if strcmp(geneFull{i},'terminator') % does not do anything yet
        protTERMflag = true;
    end
end
if length(geneFull)>=1 % construct gene string
    genestr = geneFull{1}; 
    justGene = geneFull{1}; %assuming the format is gene-lva-...-terminator
end
if length(geneFull)>=2
    for i = 2:length(geneFull)
        genestr = [genestr '-' geneFull{i}]; %expect: gene-lva-terminator
    end
end
protstr = ['protein ' genestr]; % protstr looks something like 'protein tetR-lva-terminator'
protein = sbioselect(tube, 'Type', 'species', 'Name', protstr);
if isempty(protein)
protein = addspecies(tube, protstr);
end
if exist(['txtl_protein_' justGene]) == 2
  % Run the protein specific setup
  [Rlist, genelen] = eval(['txtl_protein_' justGene '(tube, protein, geneFull, genelen)']);
end
% protein lengths
if length(geneFull)>=1
    genelenTot = genelen{1};
end
if length(geneFull)>=2 
    for i = 2:length(geneFull)
        genelenTot = genelenTot+ genelen{i};
    end
end
protein.UserData = genelenTot / 3;
%% Untranslated Region properties, parameters and reactions 
if length(rbsFull)>=1 % construct rbs string
    rbsstr = rbsFull{1};
    justRbs = rbsFull{1};
end
if length(rbsFull)>=2
    for i = 2:length(rbsFull)
        rbsstr = [rbsstr '-' rbsFull{i}]; %expect: rbs-spacer
    end
end
rnastr = ['RNA ' rbsstr '--' genestr];
rna = addspecies(tube, rnastr);
% Translation: setup file should return pointer to RBS bound species
if exist(['txtl_utr_' justRbs]) == 2
  % Run the RBS specific setup
  [Ribobound, rbslen] = eval(['txtl_utr_' justRbs '(tube, rna, protein, rbsFull, rbslen)']);
else
  % Issue a warning and run the default RBS
  warning(['TXTL: can''t find txtl_utr_' justRbs ...
      '; using default rbs params']);
  [Ribobound, rbslen] = txtl_utr_rbs(tube, rna, protein, rbsFull, rbslen);
end
% utr lengths
if length(rbsFull)>=1
    rbslenTot = rbslen{1};
end
if length(rbsFull)>=2 
    for i = 2:length(rbsFull)
        rbslenTot = rbslenTot+ rbslen{i};
    end
end

rna.UserData = rbslenTot + genelenTot;

%% Promoter properties, parameters and reactions 
if length(promFull)>=1 
    promstr = promFull{1};
    justProm = promFull{end}; % assuming {'thio','junk','prom'}
end
if length(promFull)>=2
    for i = 2:length(promFull)
        promstr = [promstr '-' promFull{i}]; % expect: 'thio-junk-prom'
    end
end
dnastr = ['DNA ' promstr '--' rbsstr '--' genestr];
dna = addspecies(tube, dnastr, amount);

% Transcription
if exist(['txtl_prom_' justProm], 'file') == 2    
  [~, promlen] = eval(['txtl_prom_' justProm '(tube, dna, rna, promFull, promlen)']);
else
  warning(['TXTL: can''t find txtl_prom_' justProm ...
      '; using default promoter params']);
  [~, promlen] = txtl_prom_p70(tube, dna, rna, promFull, promlen);
end

% promoter lengths
if length(promFull)>=1
    promlenTot = promlen{1};
end
if length(promFull)>=2 
    for i = 2:length(promFull)
        promlenTot = promlenTot+ promlen{i};
    end
end

%junk and thio dna
junklength = 0;
thiolength = 0;
junkDNAflag = false;
thioDNAflag = false;
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
if thioDNAflag
    kDNA_complex_deg = 0.5*kDNA_complex_deg;
end

% total dna length
dna.UserData = promlenTot + rbslenTot + genelenTot;

%% General purpose reactions
% DNA degradation
%! TODO: eventually, we should allow file-based reactions
%! TODO: allow protection strands by comparing promoter length to ~35
if strcmp(type, 'linear')
   Rlist = txtl_dna_degradation(tube, dna, [kDNA_recbcd_f, kDNA_recbcd_r, kDNA_complex_deg]); 
end


% Now put in the reactions for the utilization of amino acids 

% Set up the translation reaction
AA_model = 1;

if AA_model == 1
    
    aacnt = floor(protein.UserData/100);	% get number of K amino acids
    if (aacnt == 0) 
      aastr = '';
    else
      aastr = int2str(aacnt);
    end
    Robj = addreaction(tube, ...
      ['[' Ribobound.Name '] + ' aastr ' AA <-> [AA:' Ribobound.Name ']']);
    Kobj = addkineticlaw(Robj, 'MassAction');
    set(Kobj, 'ParameterVariableNames', {'TXTL_AA_F', 'TXTL_AA_R'});
else
    Robj = addreaction(tube, ...
      ['[' Ribobound.Name '] + AA <-> [AA:' Ribobound.Name ']']);
    Kobj = addkineticlaw(Robj, 'MassAction');
    set(Kobj, 'ParameterVariableNames', {'TXTL_AA_F', 'TXTL_AA_R'});
    
    Robj3 = addreaction(tube, ...
     ['[AA:' Ribobound.Name '] -> ' rna.Name ' +  Ribo']);
    Kobj3 = addkineticlaw(Robj3, 'MassAction');
    %generate unique parameter name for the current protein
    rN = regexprep(protein.Name, {'( )'}, {''});
    uniqueName = sprintf('TXTL_TL_rate_%s_AA_consumption',rN);
    set(Kobj3, 'ParameterVariableNames', uniqueName);
    
    
end
Robj1 = addreaction(tube, ...
  ['[AA:' Ribobound.Name '] -> ' rna.Name ' + ' protein.Name ' +  Ribo']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
%generate unique parameter name for the current protein
rN = regexprep(protein.Name, {'( )'}, {''});
uniqueName = sprintf('TXTL_TL_rate_%s',rN);
set(Kobj1, 'ParameterVariableNames', uniqueName);

% Add in mRNA degradation reactions
Robj2 = addreaction(tube, [rna.Name ' + RNase -> RNase']);
Kobj2 = addkineticlaw(Robj2,'MassAction');
set(Kobj2, 'ParameterVariableNames', {'TXTL_RNAdeg_F'});

% Protein degradation (if tagged)
if protDEGflag
  degradationRate = [0.005 0.001 0.001]; 
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
          if ~A(1) % this happens when the name of the dna fragment does not start with an alphabet
              error('txtl_adddna:wrongSpeciesName',...
                  ['species named %s should start with an alphabet. Format is' ...
                  ' name(length). Where the lengths are optional. eg: thio or junk(500)'],indivRegions{i})
          end
      end
      % return the parsed name and optional length
      names{i} = namesAndLengths{i}{1};
      if length(namesAndLengths{i}) == 1
          lengths{i} = [];
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
  % !TODO add error checking for numerical values for the lengths. 
  return

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
