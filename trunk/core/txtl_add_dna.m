%TXTL_ADD_DNA   Set up species and reactions for a DNA segment
%
%
%   dna = TXTL_ADD_DNA(tube, prom_spec, rbs_spec, gene_spec, amount, type)
%   constructs the species and reactions required for transcription,
%   translation and degradation of DNA, mRNA and proteins in the
%   TX-TL system.
%
%   * tube = Simbiology model object
%   * preprom_spec = Cell array of nucleatide sequences and corresponding
%   sizes. One example of their use is as a protection from exonucleases.
%   * prom_spec = spec of the form 'pre_prom(nn)'-'prom(nn)' where 'prom' is the
%     promoter name and 'len' is the length of the promoter. pre_prom cound
%     consist of nucleatide sequences and corresponding
%   sizes. One example of their use is as a protection from exonucleases.
%   * rbs_spec = spec of the form 'rbs(nn)' where 'rbs' is the RBS
%     name and 'len' is the length of the RBS.
%   * gene_spec = spec of the form 'gene(nn)-lva(nn)-terminator(nn)' where 'gene' is the
%     gene name and 'len' is the length of the gene.
%   * amount = amount of DNA to put in the tube (in nM)
%   * type = 'linear' if you want to include degradation reactions

% Written by Richard Murray, Sep 2012
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
%%
function dna = txtl_add_dna(tube, prom_spec, rbs_spec, gene_spec, dna_amount, type, varargin)
mode = struct('add_dna_driver', {[]},...
              'prot_deg_flag',{false},...
              'no_protein_flag',{false}, ...
              'prot_term_flag',{false},...
              'utr_attenuator_flag',{false},...
              'utr_antisense_flag',{false},...
              'utr_rbs_flag',{false}, ...
              'prom_junk_flag',{false},'prom_thio_flag',{false}, ...
              'sim_module_exception', {false}, 'double_antisense', {false}); 
          %generalise code to remove sim module. 
          
          
% Extract out the names and lengths
[promData, promStr] = txtl_parsespec(prom_spec);
[utrData, utrStr] = txtl_parsespec(rbs_spec);
[geneData, geneStr] = txtl_parsespec(gene_spec);
%utrData is a cell array, 2 x n, n = num of utr domains. 1st row:
%names ('att', 'rbs'). second: lengths.

% check for variations in DNA, used to select specific code
mode.prot_deg_flag = checkForStringInACellList(geneData(1,:),'lva');
mode.prot_term_flag = checkForStringInACellList(geneData(1,:),'terminator');
mode.no_protein_flag = checkForStringInACellList(geneData(1,:),'no_protein');
mode.utr_attenuator_flag = checkForStringInACellList(utrData(1,:),{'att1', 'att2'});
mode.utr_antisense_flag = checkForStringInACellList(utrData(1,:),{'anti1', 'anti2'});
mode.utr_rbs_flag = checkForStringInACellList(utrData(1,:),'rbs');
[mode.prom_junk_flag, junkIndex] = checkForStringInACellList(promData(1,:),'junk');
mode.prom_thio_flag = checkForStringInACellList(promData(1,:),'thio');

[BooleanAtt1Present,indexesOfAtt1Present] = checkForStringInACellList(utrData(1,:),'att1');
[BooleanAnti1Present,indexesOfAnti1Present] = checkForStringInACellList(utrData(1,:),'anti1');
if length(indexesOfAtt1Present)==2 
    % %! TODO: VS 11/13/13 in the more general setup, the sim module exception should not need to be handled so exceptionally
    % how to implement this generally? may be worth discussing.
    if indexesOfAtt1Present == [1 2]
        mode.sim_module_exception = true;
    end
end
if length(indexesOfAnti1Present)==2
    if indexesOfAnti1Present == [2 3]
        mode.double_antisense = true;
    end
end

%! TODO: (VS 4/17/13) need to add to this if other new types of domains are added to which ribisome can bind.

% species name string building
geneName = geneData{1,1}; %assuming the format is gene-lva-...-terminator
%! TODO: VS 4/17/13 Is rbs always the last entry of utr? what about
%spacer?

rbsName = utrData{1,end}; % format is att-...-rbs.
if mode.no_protein_flag
    protstr = ['protein ' geneStr];
    rnastr = ['RNA ' utrStr];
    dnastr = ['DNA ' promStr '--' utrStr];

else
    protstr = ['protein ' geneStr]; % protstr looks something like 'protein tetR-lva-terminator'
    rnastr = ['RNA ' utrStr '--' geneStr];
    dnastr = ['DNA ' promStr '--' utrStr '--' geneStr];
end

    %! TODO: VS 11/13/13 Enable having anything after antisense. 
% the steps here are to modify these flags here so the antisense can have
% stuff after it, and in the tx_cascade file, the antisense att interaction
% need not have the antisense the last thing on the RNA. as long as
% antisense is present. this can be recognized using regexp. 


promoterName = promData{1,end}; % assuming {'thio','junk','prom'}




%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    mode.add_dna_driver = 'Setup Species';
    tubeUser = get(tube, 'UserData');
    dnaList = tubeUser.DNAinfo;
    dnaList{end+1} = {prom_spec, rbs_spec, gene_spec, dna_amount, type, 'rxns_not_set_up', mode};
    tubeUser.DNAinfo = dnaList;
    set(tube, 'UserData', tubeUser)
    clear dnaList tubeUser
    
    % set up protein reactions and data, followed by utr followed by promoter
    % (promoter reactions require the lengths of the rna, therefore need to be
    % set up after the protein and utr files are called.
    
    
    %% Protein properties, parameters and reactions %
    %! TODO: VS 11/13/13 implement multiple proteins in genestring, which
    %means this needs to change. 
    protein = txtl_addspecies(tube, protstr, 0, 'Internal');
    if exist(['txtl_protein_' geneName], 'file') == 2
        geneData = eval(['txtl_protein_' geneName '(mode, tube, protein, geneData)']);
    end
    genelenTot = sum(cell2mat(geneData(2,:)));
    protein.UserData = genelenTot / 3;
    
    %% Untranslated Region properties, parameters and reactions %
    %! TODO: VS 11/13/13 all this will need to be changed in the more
    % general setting. need a file based system. 
    rna = txtl_addspecies(tube, rnastr, 0, 'Internal');
    if exist(['txtl_utr_' rbsName], 'file') == 2
        if mode.utr_rbs_flag 
            %ribosome must bind
            [Ribobound, utrlen] = eval(['txtl_utr_' rbsName '(mode, tube, rna, protein, utrData)']);
        else
            % no rbs present, no ribosome binds. (so far there are no files apart from txtl_utr_rbs
            %! TODO: VS 11/13/13 can do this better, considering all the
            %possible cases, and that rbs may not always be the last entry.
            [utrlen] = eval(['txtl_utr_' rbsName '(mode, tube, rna, protein, utrData)']);
        end
    else 
        if mode.utr_rbs_flag 
            %file doesnt exist and rbs flag is set. basically never happens, since rbs file is present. 
            warning('txtltoolbox:txtl_add_dna:fileNotFound', ['TXTL: can''t find txtl_utr_' rbsName ...
                '; using default rbs params']);
            [Ribobound, utrlen] = txtl_utr_rbs(mode, tube, rna, protein, utrData);
        else
            if ~(mode.utr_antisense_flag || strcmp(rbsName, 'control'))
            warning('txtltoolbox:txtl_add_dna:fileNotFound', ['TXTL: can''t find txtl_utr_' rbsName ...
                '; using default rbs params']);
            end
            [utrlen] = txtl_utr_rbs(mode, tube, rna, protein, utrData);
        end
    end
    
    utrlenTot = sum(cell2mat(utrlen(2,:)));
    if mode.utr_attenuator_flag
        if mode.sim_module_exception % SIM module exception
                attFirstLength = cell2mat(utrlen(2,1)); 
                attSecondLength = cell2mat(utrlen(2,2)); 
                restOfRNALength = utrlenTot + genelenTot - attFirstLength - attSecondLength;
                RNAlengthData.attFirst = attFirstLength;
                RNAlengthData.attSecond = attSecondLength;
                RNAlengthData.remaining = restOfRNALength;
            
        else
        attenuatorLength = cell2mat(utrlen(2,1)); % assume att is always the first thing in the UTR 
        %! TODO: VS 11/13/13 the abose assumption also not necessarily true
        %anymore. all this needs changing!
        restOfRNALength = utrlenTot + genelenTot - attenuatorLength;
        RNAlengthData.att = attenuatorLength;
        RNAlengthData.remaining = restOfRNALength;
        end
    else
        RNAlengthData = utrlenTot + genelenTot;
    end
    
    rna.UserData = RNAlengthData;
    clear RNAlengthData attenuatorLength restOfRNALength
    
    %% Promoter properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%
    % DNA solution is 22.5% of the 10ul reaction volume
    stockMulti = 10/2.25;
    dna_amount = dna_amount*stockMulti;
    dna = txtl_addspecies(tube, dnastr, dna_amount, 'Internal');
    % Transcription %%
    if exist(['txtl_prom_' promoterName], 'file') == 2 
        if mode.utr_attenuator_flag
        promData = eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, promData, prom_spec, rbs_spec, gene_spec)']);
        else
            promData = eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, promData)']);
        end
    else
        warning(['TXTL: can''t find txtl_prom_' promoterName ...
            '; using default promoter params']);
        promData = txtl_prom_p70(mode, tube, dna, rna, promData);
    end
    promlenTot = sum(cell2mat(promData(2,:)));
    %     if mode.utr_antisense_flag
    %         dna.UserData = promlenTot + utrlenTot;
    %     else
    %         dna.UserData = promlenTot + utrlenTot + genelenTot;
    %     end
    dna.UserData = promlenTot + utrlenTot + genelenTot;
    
    %% Translation %%
    if mode.utr_rbs_flag % no translation if there is no ribosome binding site present
        txtl_translation(mode, tube, dna, rna, protein, Ribobound);
    end
    
    %% DNA, protein degradation
    
    % DNA degradation
    if strcmp(type, 'linear') 
        % set up species for linear DNA degradation
        % Extract is 1/3 of the 10ul reaction volume
        extractStockMulti = 10/(10/3);
        % Add in exonuclease 
        txtl_addspecies(tube, 'RecBCD', extractStockMulti*5);	% % 0.3 from Clare Chen (Jongmin will provide ref)
        
        txtl_dna_degradation(mode, tube, dna);
    end
    
    % Protein degradation (if tagged)
    if mode.prot_deg_flag
        txtl_protein_degradation(mode, tube, protein);
    end
    
    % All done!
    return
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(varargin{1}, 'Setup Reactions')
    mode.add_dna_driver = varargin{1};
    % get a list of the species to search through before setting up
    % certain reactions
    listOfSpecies = get(tube.species, 'name');
    
    
    % set up protein reactions and data, followed by utr followed by promoter
    % (promoter reactions require the lengths of the rna, therefore need to be
    % set up after the protein and utr files are called.
    
    %% Protein properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%%
    
    protein = sbioselect(tube, 'Name', protstr);
    
    if exist(['txtl_protein_' geneName], 'file') == 2
        % Run the protein specific setup
        eval(['txtl_protein_' geneName '(mode, tube, protein, listOfSpecies)']);
    end
    
    %% Untranslated Region properties, parameters and reactions %%%%%%%%%%%
    
    rna = sbioselect(tube, 'Name', rnastr);
    if exist(['txtl_utr_' rbsName], 'file') == 2
        eval(['txtl_utr_' rbsName '(mode, tube, rna, protein)']);
    else
        if ~(mode.utr_antisense_flag || strcmp(rbsName, 'control'))
            warning('txtltoolbox:txtl_add_dna:fileNotFound', ['TXTL: can''t find txtl_utr_' rbsName ...
                '; using default rbs params']);
        end
        
        txtl_utr_rbs(mode, tube, rna, protein, utrData);
    end
    
    %% Promoter properties, parameters and reactions %%%%%%%%%%%%%%%%%%%%%%
    
    dna = sbioselect(tube, 'Name', dnastr);
    % Transcription %%
    if exist(['txtl_prom_' promoterName], 'file') == 2 
        if mode.utr_attenuator_flag
        eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, listOfSpecies,prom_spec, rbs_spec, gene_spec)']);
        else
            eval(['txtl_prom_' promoterName '(mode, tube, dna, rna, listOfSpecies)']);
        end
    else
        warning(['TXTL: can''t find txtl_prom_' promoterName ...
            '; using default promoter params']);
        txtl_prom_p70(mode, tube, dna, rna, listOfSpecies);
    end
    
    % Translation %%
    if mode.utr_rbs_flag % no translation unless rbs present. (no rbs => no Ribobound)
        Ribobound = sbioselect(tube, 'Name', ['Ribo:' rna.Name]);
        txtl_translation(mode, tube, dna, rna, protein, Ribobound);
    end
    %% DNA, mRNA, protein degradation
    
    % DNA degradation
    if strcmp(type, 'linear')
        if mode.prom_junk_flag
            junkLength = promData{2,junkIndex};
            kDNA_complex_deg = log(2)/(1+junkLength/100);
        else
            kDNA_complex_deg = tube.UserData.ReactionConfig.DNA_RecBCD_complex_deg;
        end
        if mode.prom_thio_flag
            kDNA_complex_deg = 0.5*kDNA_complex_deg;
        end
        
        % forward rr for DNA + RecBCD <-> DNA:RecBCD
        kDNA_recbcd_f = tube.UserData.ReactionConfig.DNA_RecBCD_Forward;
        % backward rr for DNA + RecBCD <-> DNA:RecBCD
        kDNA_recbcd_r = tube.UserData.ReactionConfig.DNA_RecBCD_Reverse;
        
        % linear DNA protection by gamS
        txtl_addreaction(tube,'RecBCD + [protein gamS] <-> RecBCD:gamS',...
            'MassAction',{'GamS_RecBCD_f',tube.UserData.ReactionConfig.GamS_RecBCD_F;...
            'GamS_RecBCD_r',tube.UserData.ReactionConfig.GamS_RecBCD_F/5});
        
        txtl_dna_degradation(mode, tube, dna, [kDNA_recbcd_f, kDNA_recbcd_r, kDNA_complex_deg]);
    end
    
    % Add in mRNA degradation reactions
    txtl_mrna_degradation(mode, tube, dna, rna, rbs_spec);
    
    % Protein degradation (if tagged)
    if mode.prot_deg_flag
        degradationRate = ...
            [tube.UserData.ReactionConfig.Protein_ClpXP_Forward tube.UserData.ReactionConfig.Protein_ClpXP_Reverse...
            tube.UserData.ReactionConfig.Protein_ClpXP_complex_deg];
        txtl_protein_degradation(mode, tube, protein,degradationRate);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_add_dna:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end


end % end of function


function [binVariable,indexes] = checkForStringInACellList(cellList,matchStr)
FlagVector = cellfun(@(x) strcmp(x,matchStr),cellList,'UniformOutput',false);
indexes = find(cell2mat(FlagVector) > 0);
if sum(cell2mat(FlagVector)) >= 1
    binVariable = true;
else
    binVariable = false;
end
end

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
