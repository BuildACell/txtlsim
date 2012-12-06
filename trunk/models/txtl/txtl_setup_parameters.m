% Standard parameters for the TXTL system
function txtl_setup_parameters(modelObj)

% Set up parameters
kf_aa = log(2) / 0.001;			% binding rate of 1 ms
kr_aa = 1 * kf_aa;			% Km of 100 for amino acid usage
addparameter(modelObj, 'TXTL_AA_F', kf_aa);
addparameter(modelObj, 'TXTL_AA_R', kr_aa);


% Parameters used in this file
%! TODO: update these parameters to something reasonable

kRNA_deg = log(2)/10;			% mRNA degradation: 10 sec half life
addparameter(modelObj, 'TXTL_RNAdeg_F', kRNA_deg);


% Parameters that describe this promoter
%! TODO: replace these values with correct values

kf_ntp = log(2) / 0.001;		% binding rate of 1 ms
kr_ntp = 1 * kf_ntp;			% Km of 100 for NTP usage
addparameter(modelObj, 'TXTL_NTP_RNAP_F', kf_ntp);
addparameter(modelObj, 'TXTL_NTP_RNAP_R', kr_ntp);


% Add RNAP+Sigma70 <-> RNAP70 reaction
Kf = 100; Kr = 0.01;
addparameter(modelObj, 'TXTL_RNAP_S70_F', Kf);
addparameter(modelObj, 'TXTL_RNAP_S70_R', Kr);

% TX
% sweep through the RNAs and set up reaction rates
%! TODO find a better way to access the list of species
[~,listOfSpecies] = getstoichmatrix(modelObj);
matchStr = regexp(listOfSpecies,'(^RNA .*)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
listOfRNAs = vertcat(matchStr{:});
speciesIndex = findspecies(modelObj,listOfRNAs');
for k = 1:size(speciesIndex,2)
    % TX parameters
    ktx = log(2)/(modelObj.Species(speciesIndex(k)).UserData/30);		% 30 NTP/second transcription
    rN = regexprep(modelObj.Species(speciesIndex(k)).Name, {'( )'}, {''});
    uniqueName = sprintf('TXTL_TX_rate_%s',rN);
    addparameter(modelObj, uniqueName, ktx);
    % dummy reaction for NTP consumption 
    ntpcnt = floor(modelObj.Species(speciesIndex(k)).UserData/100);	% get number of NTP blocks
    ntp_consump_rate = (ntpcnt-1)*ktx;
    uniqueName = sprintf('TXTL_TX_rate_%s_NTP_consumption',rN);
    addparameter(modelObj, uniqueName, ntp_consump_rate);
end

%TL
% same procedure for proteins
matchStr = regexp(listOfSpecies,'(^protein .*)','tokens','once');
listOfProteins = vertcat(matchStr{:});
speciesIndex = findspecies(modelObj,listOfProteins');
for k = 1:size(speciesIndex,2)
    if ~isempty(modelObj.Species(speciesIndex(k)).UserData)
        % !TODO: check how the adding of degradation tag and termination site affects this.
        ktl_rbs = log(2)/(modelObj.Species(speciesIndex(k)).UserData/10);	% 10 AA/second translation 
        rN = regexprep(modelObj.Species(speciesIndex(k)).Name, {'( )'}, {''});
        uniqueName = sprintf('TXTL_TL_rate_%s',rN);
        addparameter(modelObj, uniqueName, ktl_rbs);
        
        % dummy raction for AA consumption
        aacnt = floor(modelObj.Species(speciesIndex(k)).UserData/100);
        aa_consump_rate = (aacnt-1)*ktl_rbs;
        uniqueName = sprintf('TXTL_TL_rate_%s_AA_consumption',rN);
        addparameter(modelObj, uniqueName, aa_consump_rate);
        
    end
end

%% promoter specific parameter setup
% STANDARD PROMOTERS: ptet, p70, placI, ptrc2

    % DNA + RNAP <-> RNAPbound
        % ptet
        kf_ptet = log(2)/0.1;			% 100 ms bind rate
        kr_ptet = 10 * kf_ptet;			% Km of 10 (same as p70, from VN)
        addparameter(modelObj, 'TXTL_PTET_RNAPbound_F', kf_ptet);
        addparameter(modelObj, 'TXTL_PTET_RNAPbound_R', kr_ptet);
        
        % p70       
        kf_p70 = log(2)/0.1;			% 100 ms bind rate
        kr_p70 = 10 * kf_p70;			% Km of 10 nM (from VN model)
        addparameter(modelObj, 'TXTL_P70_RNAPbound_F', kf_p70);
        addparameter(modelObj, 'TXTL_P70_RNAPbound_R', kr_p70);
        
        % placI       
        kf_placI = log(2)/0.1;			% 100 ms bind rate
        kr_placI = 10 * kf_placI;		% Km of 10 (same as p70, from VN)
        addparameter(modelObj, 'TXTL_PLACI_RNAPbound_F', kf_placI);
        addparameter(modelObj, 'TXTL_PLACI_RNAPbound_R', kr_placI);
        
        % ptrc2        
        kf_ptrc2 = log(2)/0.1;			% 100 ms bind rate
        kr_ptrc2 = 10 * kf_ptrc2;		% Km of 10 (same as p70, from VN)
        addparameter(modelObj, 'TXTL_PTRC2_RNAPbound_F', kf_ptrc2);
        addparameter(modelObj, 'TXTL_PTRC2_RNAPbound_R', kr_ptrc2);
        
    % Promoter protein interactions
        % ptet
        ptetRepression = false;
        kf_ptet_repression_singledimer1 = 0.0001; % do sweeps over these models to see
        kr_ptet_repression_singledimer1 = 0.0001;
        matchStr = regexp(listOfSpecies,'(^protein tetR.*dimer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
        listOftetRdimer = vertcat(matchStr{:});
        if ~isempty(listOftetRdimer)
            ptetRepression = true;
        end
        
        if ptetRepression
            for i = 1:size(listOftetRdimer,1)
                rN = regexprep(listOftetRdimer{i}, {'( )'}, {''});
                uniqueNameF = sprintf('TXTL_PTET_REPRESSION1_%s_F',rN);
                uniqueNameR = sprintf('TXTL_PTET_REPRESSION1_%s_R',rN);
                addparameter(modelObj, uniqueNameF, kf_ptet_repression_singledimer1);
                addparameter(modelObj, uniqueNameR, kr_ptet_repression_singledimer1);
                %the remaining repression reactions are commented out in
                %the ptet file to keep things simple for now. Adding them
                %involves simply copying the existing template. 
            end
        end
        
        % p70
        
        % placI
        placIRepression = false;
        kf_placI_repression_tetramer = 0.00001;
        kr_placI_repression_tetramer = 0.0001;
        matchStr = regexp(listOfSpecies,'(^protein lacI.*tetramer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
        listOflacItetramers = vertcat(matchStr{:});
        if ~isempty(listOflacItetramers)
            placIRepression = true;
        end
        
        %bind tetramer to a single side
        if placIRepression
            for i = 1:size(listOflacItetramers,1)                
                rN = regexprep(listOflacItetramers{i}, {'( )'}, {''});
                uniqueNameF = sprintf('TXTL_PLACI_REPRESSION_%s_F',rN);
                uniqueNameR = sprintf('TXTL_PLACI_REPRESSION_%s_R',rN);
                addparameter(modelObj, uniqueNameF, kf_placI_repression_tetramer);
                addparameter(modelObj, uniqueNameR, kr_placI_repression_tetramer);
            end
        end
        
        % ptrc2
        ptrc2Repression = false;
        kf_ptrc2_repression_tetramer = 0.50;
        kr_ptrc2_repression_tetramer = 0.03;
        %! TODO make all these reactions conditional on specie availability
        matchStr = regexp(listOfSpecies,'(^protein lacI.*dimer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
        listOflacIdimer = vertcat(matchStr{:});
        if ~isempty(listOflacIdimer)
            ptrc2Repression = true;
        end
        
        if ptrc2Repression
            for i = 1:size(listOflacIdimer,1)
                rN = regexprep(listOflacIdimer{i}, {'( )'}, {''});
                uniqueNameF = sprintf('TXTL_PTRC2_REPRESSION_%s_F',rN);
                uniqueNameR = sprintf('TXTL_PTRC2_REPRESSION_%s_R',rN);
                addparameter(modelObj, uniqueNameF, kf_ptrc2_repression_tetramer);
                addparameter(modelObj, uniqueNameR, kr_ptrc2_repression_tetramer);
            end
        end
      


% ADDITIONAL PROMOTERS (may be user defined)


%% utr specific parameter setup
% STANDARD UTR: rbs
kf_rbs = log(2)/0.1;			% 100 ms bind rate
kr_rbs = 0.05 * kf_rbs;			% Km of ~0.05 (from VN model)
addparameter(modelObj, 'TXTL_UTR_RBS_F', kf_rbs);
addparameter(modelObj, 'TXTL_UTR_RBS_R', kr_rbs);

% ADDITIONAL UTR (may be user defined)


%% gene specific parameter setup 
% STANDARD GENES: tetR, lacI, deGFP, gamS

    %tetR
    kf_tetR_aTc = 1; kr_tetR_aTc = 0.1;
    kf_aTcdeg = 0.0001;
    addparameter(modelObj, 'TXTL_INDUCER_TETR_ATC_F', kf_tetR_aTc);
    addparameter(modelObj, 'TXTL_INDUCER_TETR_ATC_R', kr_tetR_aTc);
    addparameter(modelObj, 'TXTL_INDUCER_DEGRADATION_ATC', kf_aTcdeg);
    
    %lacI
    kf_lacI_IPTG = 0.1; kr__lacI_IPTG = 0.01;
    kf_IPTGdeg = 0.0001;
    addparameter(modelObj, 'TXTL_INDUCER_LACI_IPTG_F', kf_lacI_IPTG);
    addparameter(modelObj, 'TXTL_INDUCER_LACI_IPTG_R', kr__lacI_IPTG);
    addparameter(modelObj, 'TXTL_INDUCER_DEGRADATION_IPTG', kf_IPTGdeg);
    
    % deGFP
    Kmat_deGFP = log(2)/(15*60);
    addparameter(modelObj, 'TXTL_PROT_DEGFP_MATURATION', Kmat_deGFP);
    
    

% ADDITIONAL GENES (may be user defined)

%%
% OTHER REACTIONS
% extract promoter, utr and gene strings, along with their lengths. 
    

% tetR dimerization
kf_tetR_dimer = 8e-2;% originally 1e-4 1/(molecule*sec)
kr_tetR_dimer = 0.02; %0.00000001; % 1/sec
matchStrings = regexp(listOfSpecies,'(^protein tetR.*)', 'match');
listOftetRSpecies = vertcat(matchStrings{:});
countOftetRMonomers = 1;
for i = 1:length(listOftetRSpecies)
    if length(listOftetRSpecies{i}) > 4
        str1 = listOftetRSpecies{i}(end-4:end);
    end
    if length(listOftetRSpecies{i}) > 7
        str2 = listOftetRSpecies{i}(end-7:end);
    end
    if ~strcmp(str1, 'dimer') &&  ~strcmp(str2, 'tetramer')
        listOftetRMonomers{countOftetRMonomers} = listOftetRSpecies{i};
        countOftetRMonomers = countOftetRMonomers +1;
    end
end
listOftetRMonomers = listOftetRMonomers';
if ~isempty(listOftetRMonomers)
    tetRDimerization = true;
end
if tetRDimerization
    for i = 1:size(listOftetRMonomers,1)
        rN = regexprep(listOftetRMonomers{i}, {'( )'}, {''});
        uniqueNameF = sprintf('TXTL_PROT_DIMER_%s_F',rN);
        uniqueNameR = sprintf('TXTL_PROT_DIMER_%s_R',rN);
        addparameter(modelObj, uniqueNameF, kf_tetR_dimer);
        addparameter(modelObj, uniqueNameR, kr_tetR_dimer);
    end
end

% lacI dimerization, specify both forward and backward rr
kf_lacI_dimer = 8e-2;% originally 1e-4 1/(molecule*sec)
kr_lacI_dimer = 0.02; %0.00000001; % 1/sec


matchStrings = regexp(listOfSpecies,'(^protein lacI.*)', 'match');
listOflacISpecies = vertcat(matchStrings{:});
countOflacIMonomers = 1;
listOflacIMonomers = cell(0);
lacIDimerization = false;
for i = 1:length(listOflacISpecies)
    if length(listOflacISpecies{i}) > 4
        str1 = listOflacISpecies{i}(end-4:end);
    end
    if length(listOflacISpecies{i}) > 7
        str2 = listOflacISpecies{i}(end-7:end);
    end
    if length(listOflacISpecies{i}) > 5
        str3 = listOflacISpecies{i}(end-5:end);
    end    
    if ~strcmp(str1, 'dimer') &&  ~strcmp(str2, 'tetramer') && ~strcmp(str3, ':ClpXP') % i think putting ClpXP at the beginning solves this issue
        listOflacIMonomers{countOflacIMonomers} = listOflacISpecies{i};
        countOflacIMonomers = countOflacIMonomers +1;
    end
end
listOflacIMonomers = listOflacIMonomers';
if ~isempty(listOflacIMonomers)
    lacIDimerization = true;
end
if lacIDimerization
    for i = 1:size(listOflacIMonomers,1)
        rN = regexprep(listOflacIMonomers{i}, {'( )'}, {''});
        uniqueNameF = sprintf('TXTL_PROT_DIMER_%s_F',rN);
        uniqueNameR = sprintf('TXTL_PROT_DIMER_%s_R',rN);
        addparameter(modelObj, uniqueNameF, kf_lacI_dimer);
        addparameter(modelObj, uniqueNameR, kr_lacI_dimer);
    end
end
        
% lacI tetramerization

kf_lacI_tetramer = 0.00000602; % 1/(molecule*sec)
kr_lacI_tetramer = 0.000001; %0.000001; % 1/sec
matchStr = regexp(listOfSpecies,'(^protein lacI.*dimer$)','tokens','once');
listOflacIdimer = vertcat(matchStr{:});
lacItetramerization = false;
if ~isempty(listOflacIdimer)
    lacItetramerization = true;
end

if lacItetramerization
    for i = 1:size(listOflacIdimer,1)
        rN = regexprep(listOflacIdimer{i}, {'( )'}, {''});
        uniqueNameF = sprintf('TXTL_PROT_TETRAMER_%s_F',rN);
        uniqueNameR = sprintf('TXTL_PROT_TETRAMER_%s_R',rN);
        addparameter(modelObj, uniqueNameF, kf_lacI_tetramer);
        addparameter(modelObj, uniqueNameR, kr_lacI_tetramer);
    end
end

% Protein degradation
kf_degrad_complex = 0.005 ;
kr_degrad_complex = 0.001 ;
k_degrad = 0.001;
matchStrings = regexp(listOfSpecies,'(^protein .*)', 'match');
listOfProteinSpecies = vertcat(matchStrings{:});
countOfProteinMonomers = 1;
for i = 1:length(listOfProteinSpecies)
    if length(listOfProteinSpecies{i}) > 4
        str1 = listOfProteinSpecies{i}(end-4:end);
    end
    if length(listOfProteinSpecies{i}) > 7
        str2 = listOfProteinSpecies{i}(end-7:end);
    end
    if length(listOfProteinSpecies{i}) > 5
        str3 = listOfProteinSpecies{i}(end-5:end);
    end    
    if ~strcmp(str1, 'dimer') &&  ~strcmp(str2, 'tetramer') && ~strcmp(str3, ':ClpXP') % i think putting ClpXP at the beginning solves this issue
        listOfProteinMonomers{countOfProteinMonomers} = listOfProteinSpecies{i};
        countOfProteinMonomers = countOfProteinMonomers +1;
    end
end
listOfProteinMonomers = listOfProteinMonomers';
if ~isempty(listOfProteinMonomers)
    proteinDegradation = true;
end
if proteinDegradation
    for i = 1:size(listOfProteinMonomers,1)
        rN = regexprep(listOfProteinMonomers{i}, {'( )'}, {''});
        uniqueNameF = sprintf('TXTL_PROT_DEGRAD_COMPLEX_%s_F',rN);
        uniqueNameR = sprintf('TXTL_PROT_DEGRAD_COMPLEX_%s_R',rN);
        addparameter(modelObj, uniqueNameF, kf_degrad_complex);
        addparameter(modelObj, uniqueNameR, kr_degrad_complex);
        uniqueName = sprintf('TXTL_PROT_DEGRAD_%s',rN);
        addparameter(modelObj, uniqueName, k_degrad);
    end
end


end