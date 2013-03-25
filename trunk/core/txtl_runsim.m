
% Written by Zoltan A Tuza, Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Documentation:
% there are 2 ways to call this function.
% [t_ode_output, x_ode_output, simData_output] = txtl_runsim(modelObj, configsetObj, t_ode_input, x_ode_input)
% [t_ode_output, x_ode_output, simData_output] = txtl_runsim(modelObj, configsetObj, t_ode_input, x_ode_input, simData_input)
% simData_output is optional

% when using this, it is Necessary to use the model object returned by this
% function in subsequent calls / plotting. This is because the first call to
% this function sets the reactions in the modelObject.


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



function [varargout] = txtl_runsim(varargin)
%%
switch nargin
    case 1
        % parameter estimation mode, runsim just assemble the
        % reactions, but won't run it
        modelObj = varargin{1};
        configsetObj = [];
        time_vector = [];
        data =  [];
    case 2
        modelObj = varargin{1};
        configsetObj = varargin{2};
        time_vector = [];
        data =  [];
    case 4
        modelObj = varargin{1};
        configsetObj = varargin{2};
        time_vector = varargin{3};
        data = varargin{4};
    case 5
        modelObj = varargin{1};
        configsetObj = varargin{2};
        time_vector = varargin{3};
        data = varargin{4};
        simData = varargin{5};
    otherwise
        error('txtl_runsim should be called either with 4 or 5 arguments.');
end
% check what proteins are present, but no corresponding DNA. this will mean
% that the protein must have been added. So now set up the reactions for
% that protein. and you are done! So basically, you are comparing two
% lists.
setupReactionsForNewProteinAdded(modelObj)

%! TODO zoltuz 2/4/13 review this part
m = get(modelObj, 'UserData');
datalen = size(m,1);

%% FIRST RUN
% set up reaction using txtl_add_dna if this is the first run of runsim.
% However, what about the case when we add dna after the first run?
if ~strcmp(m{datalen},'notFirstRun')
    % if it is first run, userdata only has dna spec data
    % if it is second run, the last element in the userdata is a string
    % 'notFirstRun'
    for i = 1:length(m)
        if ~isa(m{i},'txtl_reaction_config')
            for j= 1:length(m{i})
                txtl_add_dna(modelObj, m{i}{j}{1}, m{i}{j}{2}, m{i}{j}{3}, m{i}{j}{4}, m{i}{j}{5}, 'Setup Reactions');
            end
        end
    end
    
    
    % ADD: compare DNA protein list and protein list. If there are proteins
    % for which there is no dna, call that protein s setup reaction fine.
   
    
    newUserData = cat(1, m, 'notFirstRun');
    set(modelObj, 'UserData', newUserData)
    
%
% ATP degration is a first order reaction.    
ntp_deg = 0.00093;
txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
        'MassAction',{'NTPdeg_F',ntp_deg});
    

    
end % end of first run


%% SUBSEQUENT RUNS
if ~isempty(time_vector) && size(time_vector,1) > 1
    prevData = zeros(size(time_vector,1),size(modelObj.Species,1));
end

% Species-data pairs is needed
if iscell(data)
    SpName = findspecies(modelObj,data{:,1});
    for k=1:size(data,1)
        if size(data{k,2},1) == 1
            %first run initial amount provided
            modelObj.Species(SpName(k)).InitialAmount = normalizeSmallAmount(data{k,2});
        elseif size(data{k,2},1) > 1
            %setting up the initial amount to the latest simulation result
            modelObj.Species(k).InitialAmount = normalizeSmallAmount(data{k,2}(end));
            % reordering the data cell matrix to make it compatible with order in modelObj
            prevData(:,SpName(k)) = data{k,2}(:);
            
        else
            % if no data was given, we issue a warning
            warning('something went wrong on data(%d,2)',k);
        end
    end
    % we have as many colums of species in the model as in data, we set the
    % inital amount to that (here the order of data matters!)

elseif size(modelObj.Species,1) == size(data,2) && size(data,1) == 1
    for k=1:size(data,2)
        modelObj.Species(k).InitialAmount = normalizeSmallAmount(data(1,k));
    end
    % we have simulation data from a previous run
elseif size(modelObj.Species,1) == size(data,2) && size(data,1) > 1
    for k=1:size(data,2)
        modelObj.Species(k).InitialAmount = normalizeSmallAmount(data(end,k));
        prevData(:,k) = data(:,k);
    end
else
    % no data was provided, no action needed
end

% initial amounts set in modelObj.Species(k).InitialAmount.
% previousdata, if any, stored in prevData.
if ~isempty(configsetObj)
    
    simData = sbiosimulate(modelObj, configsetObj);
    
    if isempty(time_vector)
        x_ode = simData.Data;
        t_ode = simData.Time;
    else
        t_ode = [time_vector; simData.Time+time_vector(end)];
        x_ode = [prevData;simData.Data];
    end
    
    % TODO zoltuz 03/05/13 we need a new copyobject for simData -> merge
    % two simData object.
    varargout{1} = simData;

else

    % parameter estimation mode, no simulation
    x_ode = [];
    t_ode = [];
    simData = [];
end

switch nargout
    case 0
        varargout{1} = [];
    case 1
        varargout{1} = simData;
    case 2
        varargout{1} = t_ode;
        varargout{2} = x_ode;
    otherwise
        error('not supported operation mode');
        
end

end

% InitialAmount cannot be smaller than certain amount, therefore small
% amounts are converted to zero
function retValue = normalizeSmallAmount(inValue)

if (abs(inValue) < eps) || inValue<0 % treats values below eps as zero.
    retValue = 0;
else
    retValue = inValue;
end

end

function setupReactionsForNewProteinAdded(tube)

% Make a list of all the proteins added by DNAs:
speciesNames = get(tube.species, 'name');
matchStr = regexp(speciesNames,'(^DNA.*)','tokens','once');
listOfDNA = vertcat(matchStr{:});
DNAparts = regexp(listOfDNA, '--','split');
proteinsAlreadySetUp = cell(length(DNAparts),1);
for i = 1:length(DNAparts)
    proteinsAlreadySetUp{i} = DNAparts{i}{3};
end

% Make a list of protein in the model:
matchStr = regexp(speciesNames,'(^protein.*)','tokens','once');
listOfprotein = vertcat(matchStr{:});
for k = 1:size(listOfprotein, 1)
    colonIdx = strfind(listOfprotein{k},':');
    if isempty(colonIdx)
        proteinName = listOfprotein{k}(9:end);
        if  size(proteinName,2) > 1 && strcmp(proteinName(end), '*') %what about deGFP*-lva-terminator? what happens there?
            proteinName = proteinName(1:end-1);
        end
        if size(proteinName,2) > 4 && strcmp(proteinName(end-4:end), 'dimer')
            proteinName = proteinName(1:end-5);
        elseif  size(proteinName,2) > 7 && strcmp(proteinName(end-7:end), 'tetramer')
            proteinName = proteinName(1:end-8);
        end
        temp = strfind(proteinsAlreadySetUp, proteinName);
        %temp2 = vertcat(temp{:});
        %clear temp
        scalarIndices = isscalar(temp);
        %clear temp2
        %            temp = strfind(proteinsAlreadySetUp, proteinName);
        %            scalarIndices = cellfun(@isscalar, temp);
        if scalarIndices == 0
            proteinParts = regexp(proteinName, '-','split');
            justProtein = proteinParts{1};
            proteinIdx = findspecies(tube, ['protein ' proteinName]);
            protein = tube.Species(proteinIdx);
            if exist(['txtl_protein_' justProtein], 'file') == 2
                % Run the protein specific setup
                protData = txtl_parsespec(proteinName);
                eval(['txtl_protein_' justProtein '(''Setup Reactions'', tube, protein, protData)']);
            else
                warning('Warning:ProteinFileNotFound','Protein %s file not defined.', justProtein)
            end
        end
    end
end


% Note a protein like 'protein tetR-lva' is different from 'protein
% tetR'


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