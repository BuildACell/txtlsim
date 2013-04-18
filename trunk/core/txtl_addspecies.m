function simBioSpecies = txtl_addspecies(tube, name, amount, varargin)
%TXTL_ADDSPECIES   Add one or more molecular species to a tube
%
% Sobj = TXTL_ADDSPECIES(tube, name, amount) adds a molecule to a
% tube, in the gen amount (in nM).  The species can be a new species or
% one that already exists (in which case its concentration is added to
% what is already present).

% Written by Richard Murray, 11 Sep 2012
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

add_dna_mode = struct('add_dna_driver', {'Setup Species'});
%%%%%%%%%%%%%%%%%%% DRIVER MODE: USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    mode = 'User';
    
    
    % check
    if (size(amount,2) > 1)
        assert(size(name,2) == size(amount,2));
    end
    
    index = findspecies(tube, name);
    
    
%     % get all proteins already set up in txtl_add_dna:
%     speciesNames = get(tube.species, 'name');
%     matchStr = regexp(speciesNames,'(^DNA.*)','tokens','once');
%     listOfDNA = vertcat(matchStr{:});
%     DNAparts = regexp(listOfDNA, '--','split');
%     proteinsAlreadySetUp = cell(length(DNAparts),1);
%     for i = 1:length(DNAparts)
%         proteinsAlreadySetUp{i} = DNAparts{i}{3};
%     end
    
    
    if iscell(name) % if more than one species are to be added
        for k =1:size(index,2)
            if size(name{k},2) > 6 && strcmp(name{k}(1:7), 'protein') %is a protein and does not yet exist
                % the protein processing function
                setupProteins(tube, name{k}, amount{k}, index(k), add_dna_mode)
            else
                addOneSpecie(tube,name{k},amount{k},index(k));
            end
        end
    else
        if size(name,2) > 6 && strcmp(name(1:7), 'protein')
            setupProteins(tube, name, amount, index, add_dna_mode)
        else
            addOneSpecie(tube,name,amount,index); % not a protien. may or may not exist in the tube!
        end
    end
    
    indexPost = findspecies(tube, name);
    simBioSpecies = tube.Species(indexPost);
    
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: INTERNAL %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(varargin{1}, 'Internal')
    mode = varargin{1};
    
    % check
    if (size(amount,2) > 1)
        assert(size(name,2) == size(amount,2));
    end
    
    index = findspecies(tube, name);
    
    if iscell(name)
        for k =1:size(index,2)
            addOneSpecie(tube,name{k},amount{k},index(k));
        end
    else
        addOneSpecie(tube,name,amount,index);
    end
    indexPost = findspecies(tube, name);
    simBioSpecies = tube.Species(indexPost);
    
end


end


function varargout = addOneSpecie(tube,name,amount,index)
% if amount wasn't specified then make it zero
if isempty(amount)
    amount = 0;
end

if (index == 0)
    varargout{1} = addspecies(tube, name, amount);
else
    tube.Species(index).InitialAmount = ...
        tube.Species(index).InitialAmount + amount;
    tube.Species(index);
    varargout{1} = tube.Species(index) ;
end
end

function [parsedData, combinedStr] = txtl_parsespec(spec)

indivRegions = regexp(spec, '-','split'); %cell array of individual xyz(123) strings
namesAndLengths = regexp(indivRegions, '\w+','match'); %cell array of cells containing the names and lengths of dna regions
names = cell(1,length(namesAndLengths));
lengths = cell(1,length(namesAndLengths));
combinedStr = '';

%error checking followed by returning parsed strings
for k = 1:length(namesAndLengths)
    if isempty(namesAndLengths{k})
        error('txtl_add_dna:wrongStringFormat',...
            ['the string %s should be: name(length)-name2(length2)-' ...
            '...-nameN(lengthN), where the lengths are optional. eg: thio-junk(500)-ptet(50)'...
            'the name must start with an alphabet'], spec)
    else
        A = isstrprop(namesAndLengths{k}{1},'alpha');
        if ~A(1) % this happens when the name of the dna fragment does not start with an alphabet
            error('txtl_add_dna:wrongSpeciesName',...
                ['species named %s should start with an alphabet. Format is' ...
                ' name(length). Where the lengths are optional. eg: thio or junk(500)'],indivRegions{k})
        end
    end
    % return the parsed name and optional length
    names{k} = namesAndLengths{k}{1};
    if length(namesAndLengths{k}) == 1
        lengths{k} = [];
    else if length(namesAndLengths{k})>2
            error('txtl_add_dna:tooManyElements',...
                ['the string %s is not of the format name(length). '...
                'It has unwanted elements after '')'''],...
                indivRegions{k});
        else if length(namesAndLengths{k})==2
                lengths{k} = str2double(namesAndLengths{k}{2});
            end
        end
    end
    if k==1
        combinedStr = names{k};
    else
        combinedStr = [combinedStr '-' names{k}];
    end
    parsedData = [names;lengths];
end
% !TODO add error checking for numerical values for the lengths.

end

function setupProteins(tube, name, amount, index, add_dna_mode)

% get all proteins already set up in txtl_add_dna:
speciesNames = get(tube.species, 'name');
matchStr = regexp(speciesNames,'(^DNA.*)','tokens','once');
listOfDNA = vertcat(matchStr{:});
DNAparts = regexp(listOfDNA, '--','split');
proteinsAlreadySetUp = cell(length(DNAparts),1);
for i = 1:length(DNAparts)
    proteinsAlreadySetUp{i} = DNAparts{i}{3};
end
% here we want to add the protein, be it the dimer,
% tetramer, whatever:
proteinSpecie = addOneSpecie(tube,name,amount,index);
proteinName = name(9:end);
if  size(proteinName,2) > 1 && strcmp(proteinName(end), '*')
    proteinName = proteinName(1:end-1);
end

if size(proteinName,2) > 4 && strcmp(proteinName(end-4:end), 'dimer')
    proteinName = proteinName(1:end-5);
elseif  size(proteinName,2) > 7 && strcmp(proteinName(end-7:end), 'tetramer')
    proteinName = proteinName(1:end-8);
end
index2 = findspecies(tube, ['protein ' proteinName]);
%add the soecie without the *, dimer or tetramer.
proteinOriginal = addOneSpecie(tube, ['protein ' proteinName], 0, index2);
temp = strfind(proteinsAlreadySetUp, proteinName);
temp2 = vertcat(temp{:});
clear temp
scalarIndices = isscalar(temp2);
clear temp2
%            temp = strfind(proteinsAlreadySetUp, proteinName);
%            scalarIndices = cellfun(@isscalar, temp);
if scalarIndices == 0
    proteinParts = regexp(proteinName, '-','split');
    justProtein = proteinParts{1};
    if exist(['txtl_protein_' justProtein], 'file') == 2
        % Run the protein specific setup
        protData = txtl_parsespec(proteinName);
        eval(['txtl_protein_' justProtein '(add_dna_mode, tube, proteinOriginal, protData)']);
    elseif ~strcmp(justProtein, 'gamS') && ~strcmp(justProtein, 'sigma70')
        warning('Warning:ProteinFileNotFound','Protein %s file not defined.', justProtein)
    end
end

% NOTES
% lets say a protein dimer is added. even in that case, we should just set
% the that protein's reactions, since that dimer will dissociate to give
% that protein, and all of the proteins reactions will follow.

% no complexes should be allowed. so aTc:protein , and protein:ClpXP should
% not be allowed to be added. This is important since the current code does
% not have the capability to handle it. This is ok in practice
% since we cant add complexes in lab anyway. (I think)

end


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
