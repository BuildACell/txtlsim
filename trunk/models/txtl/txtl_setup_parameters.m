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
kf_p70 = log(2)/0.1;			% 100 ms bind rate
kr_p70 = 10 * kf_p70;			% Km of 10 nM (from VN model)
addparameter(modelObj, 'TXTL_DNA_RNAP70_F', kf_p70);
addparameter(modelObj, 'TXTL_DNA_RNAP70_R', kr_p70);


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
matchStr = regexp(listOfSpecies,'(^RNA .*)','tokens','once');
listOfRNAs = vertcat(matchStr{:});
speciesIndex = findspecies(modelObj,listOfRNAs');
for k = 1:size(speciesIndex,2)
    ktx = log(2)/(modelObj.Species(speciesIndex(k)).UserData/30);		% 30 NTP/second transcription
    rN = regexprep(modelObj.Species(speciesIndex(k)).Name, {'( )'}, {''});
    uniqueName = sprintf('TXTL_TX_rate_%s',rN);
    addparameter(modelObj, uniqueName, ktx);
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
    end
end



end