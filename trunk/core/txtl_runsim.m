
% Written by Zoltan A Tuza, Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Documentation:
% 2 April 13
% 2 modes: if the last argument is not a string 'events', then this
% function runs in 'default' mode (See below for info). If the last
% argument to this function is a string 'events', then this function runs
% in 'events' mode. 

%%%%%% DEFAULT MODE %%%%%%
% [t_ode, x_ode, modelObj, simData] = txtl_runsim(modelObj, configsetObj, time_vector, data, simData)

% Input combinations: 
% modelObj (this runs a parameter estimation mode, no simulation)
% modelObj, configsetObj
% modelObj, configsetObj, time_vector, data
% modelObj, configsetObj, time_vector, data, simData

% Output combinations: 
% simData
% t_ode, x_ode
% t_ode, x_ode, modelObj
% t_ode, x_ode, modelObj, simData

%%%%%% EVENTS MODE %%%%%%
% [t_ode, x_ode, modelObj, simData] = txtl_runsim(modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, simData, 'events')

% Input combinations: 
% modelObj, configsetObj, eventTriggers, eventFcns, 'events'
% modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, 'events'
% modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, simData, 'events'

% Output combinations: 
% simData
% t_ode, x_ode
% t_ode, x_ode, modelObj
% t_ode, x_ode, modelObj, simData


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
if ischar(varargin{end})
    if strcmp(varargin{end},'events');
        mode = 'events';
    elseif strcmp(varargin{end},'sensitivity');
        mode = 'sensitivity';
    end
else
    mode = 'default';
end

if strcmp(mode,'default')
    %%%%%%%%%%%%%%%%%%%%%%% DEFAULT MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %! TODO zoltuz 2/4/13 review this part
    m = get(modelObj, 'UserData');
    
    % check what proteins are present, but no corresponding DNA. this will mean
    % that the protein must have been added. So now set up the reactions for
    % that protein. and you are done! So basically, you are comparing two
    % lists.
    
    add_dna_mode = struct('add_dna_driver', {'Setup Reactions'});
    setupReactionsForNewProteinAdded(modelObj, add_dna_mode)
    

    
    %% FIRST RUN
    % if m.FirstRun
    %     for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
    %         %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
    %         txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
    %             m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
    %     end
    %     m.FirstRun = False;
    % end
    
    
    for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
        %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
        if strcmp(m.DNAinfo{i}{6}, 'rxns_not_set_up')
            txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
                m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
            m.DNAinfo{i}{6} = 'rxns_already_set_up';
        end
        
        %
        % ATP degration is a first order reaction.
        ntp_deg = 0.00093;
        txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
            'MassAction',{'NTPdeg_F',ntp_deg});
        
        %is the line below necessary? in general is m a pointer or a struct?
        set(modelObj, 'UserData', m)
        
        
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
    
    %
    % After 3hours because of the ATP regeneration stops the remaining NTP
    % becomes unusable c.f. V Noireaux 2003.
    % for solver specific reason the we need some amount of "NTP_GOES_BAD",
    % otherwise the rapid transition of 0->1nM at 3hours stops the solver.
    ntp_deg = 0.00008;
    % txtl_addspecies(modelObj, 'NTP_REGEN_SUP',1, 'Internal');
    % txtl_addreaction(modelObj,'NTP_REGEN_SUP -> null',...
    %     'MassAction',{'NTP_F',0.00035});
    % txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP_UNUSE',...
    %     'MassAction',{'NTP_F',0.00035});
    
    txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
        'MassAction',{'NTPdeg_F',ntp_deg});
    %
    % txtl_addreaction(modelObj,'NTP_UNUSE + NTP_REGEN_SUP <-> NTP_UNUSE:NTP_REGEN_SUP',...
    %     'MassAction',{'NTPdeg_F',50; 'R',0.001});
    %
    % txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP + NTP_REGEN_SUP',...
    %     'MassAction',{'NTPdeg_F',30});
    
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
        case 3
            varargout{1} = t_ode;
            varargout{2} = x_ode;
            varargout{3} = modelObj;
        case 4
            varargout{1} = t_ode;
            varargout{2} = x_ode;
            varargout{3} = modelObj;
            varargout{4} = simData;
        otherwise
            error('not supported operation mode');
            
    end
elseif strcmp(mode, 'events')
    %%%%%%%%%%%%%%%%%%%%%%%%% EVENTS MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch nargin
        case 5
            modelObj = varargin{1};
            configsetObj = varargin{2};
            eventTriggers = varargin{3};
            eventFcns = varargin{4};
            time_vector = [];
            data =  [];
        case 7
            modelObj = varargin{1};
            configsetObj = varargin{2};
            eventTriggers = varargin{3};
            eventFcns = varargin{4};
            time_vector = varargin{5};
            data = varargin{6};
        case 8
            modelObj = varargin{1};
            configsetObj = varargin{2};
            eventTriggers = varargin{3};
            eventFcns = varargin{4};
            time_vector = varargin{5};
            data = varargin{6};
            simData = varargin{7};
        otherwise
            error('txtl_runsim should be called either with 2, 5, 7 or 8 arguments in events mode.');
    end
    evt = cell(1, length(eventTriggers));
    for i = 1:length(eventTriggers)
        evt{i} = addevent(modelObj, eventTriggers{i}(1:end), eventFcns{i}(1:end));
    end
    % check what proteins are present, but no corresponding DNA. this will mean
    % that the protein must have been added. So now set up the reactions for
    % that protein. and you are done! So basically, you are comparing two
    % lists.
    add_dna_mode = struct('add_dna_driver', {'Setup Reactions'});
    setupReactionsForNewProteinAdded(modelObj, add_dna_mode)
    
    %! TODO zoltuz 2/4/13 review this part
    m = get(modelObj, 'UserData');
    
    %% FIRST RUN
    % if m.FirstRun
    %     for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
    %         %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
    %         txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
    %             m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
    %     end
    %     m.FirstRun = False;
    % end
    
    
    for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
        %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
        if strcmp(m.DNAinfo{i}{6}, 'rxns_not_set_up')
            txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
                m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
            m.DNAinfo{i}{6} = 'rxns_already_set_up';
        end
        
        %
        % ATP degration is a first order reaction.
        ntp_deg = 0.00002;
        txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
            'MassAction',{'NTPdeg_F',ntp_deg});
        
        %is the line below necessary? in general is m a pointer or a struct?
        set(modelObj, 'UserData', m)
        
        
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
    
    %
    % After 3hours because of the ATP regeneration stops the remaining NTP
    % becomes unusable c.f. V Noireaux 2003.
    % for solver specific reason the we need some amount of "NTP_GOES_BAD",
    % otherwise the rapid transition of 0->1nM at 3hours stops the solver.
    ntp_deg = 0.00002;
    % txtl_addspecies(modelObj, 'NTP_REGEN_SUP',1, 'Internal');
    % txtl_addreaction(modelObj,'NTP_REGEN_SUP -> null',...
    %     'MassAction',{'NTP_F',0.00035});
    % txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP_UNUSE',...
    %     'MassAction',{'NTP_F',0.00035});
    
    txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
        'MassAction',{'NTPdeg_F',ntp_deg});
    %
    % txtl_addreaction(modelObj,'NTP_UNUSE + NTP_REGEN_SUP <-> NTP_UNUSE:NTP_REGEN_SUP',...
    %     'MassAction',{'NTPdeg_F',50; 'R',0.001});
    %
    % txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP + NTP_REGEN_SUP',...
    %     'MassAction',{'NTPdeg_F',30});
    
    % initial amounts set in modelObj.Species(k).InitialAmount.
    % previousdata, if any, stored in prevData.
    if ~isempty(configsetObj)
        set(configsetObj.RuntimeOptions, 'StatesToLog', 'all');
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
    end
    
    switch nargout
        case 0
            varargout{1} = [];
        case 1
            varargout{1} = simData;
        case 2
            varargout{1} = t_ode;
            varargout{2} = x_ode;
        case 3
            varargout{1} = t_ode;
            varargout{2} = x_ode;
            varargout{3} = modelObj;
        case 4
            varargout{1} = t_ode;
            varargout{2} = x_ode;
            varargout{3} = modelObj;
            varargout{4} = simData;
        otherwise
            error('not supported operation mode');
            
    end
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

function setupReactionsForNewProteinAdded(tube, add_dna_mode)

% Make a list of all the proteins added by DNAs:
speciesNames = get(tube.species, 'name');
matchStr = regexp(speciesNames,'(^DNA.*)','tokens','once');
listOfDNA = vertcat(matchStr{:});
DNAparts = regexp(listOfDNA, '--','split');
%proteinsAlreadySetUp = cell(length(DNAparts),1);
ii = 1;
for i = 1:length(DNAparts)
    if length(DNAparts{i})==3
        proteinsAlreadySetUp{ii} = DNAparts{i}{3};
        ii = ii+1;
    end
end
proteinsAlreadySetUp = vertcat(proteinsAlreadySetUp(:));
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
                eval(['txtl_protein_' justProtein '(add_dna_mode, tube, protein, protData)']);
            elseif ~strcmp(justProtein, 'gamS') && ~strcmp(justProtein, 'sigma70')
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