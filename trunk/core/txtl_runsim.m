
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


function [t_ode, x_ode, modelObj, varargout] = txtl_runsim(varargin)
    

    switch nargin
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
datalen = size(m,1);
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
    
    newUserData = cat(1, m, 'notFirstRun');
    set(modelObj, 'UserData', newUserData)
end
  
    
if ~isempty(time_vector) && size(time_vector,1) > 1
    prevData = zeros(size(time_vector,1),size(modelObj.Species,1));
end
% Species-data pairs is needed
if iscell(data) 
    SpName = findspecies(modelObj,data{:,1}); 
    for k=1:size(data,1)
       if size(data{k,2},1) == 1
           %first run initial amount provided
           modelObj.Species(k).InitialAmount = normalizeSmallAmount(data{k,2});
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
     txtl_addspecies(modelObj, 'NTP_REGEN_SUP',1);
      txtl_addreaction(modelObj,'NTP_REGEN_SUP -> null',...
          'MassAction',{'NTP_F',0.00035});                      
      txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP_UNUSE',...
          'MassAction',{'NTP_F',0.00035});  
    
      txtl_addreaction(modelObj,'NTP -> NTP_UNUSE',...
         'MassAction',{'NTPdeg_F',ntp_deg});
     
      txtl_addreaction(modelObj,'NTP_UNUSE + NTP_REGEN_SUP <-> NTP_UNUSE:NTP_REGEN_SUP',...
          'MassAction',{'NTPdeg_F',50; 'R',0.001});
      
      txtl_addreaction(modelObj,'NTP_UNUSE:NTP_REGEN_SUP -> NTP + NTP_REGEN_SUP',...
         'MassAction',{'NTPdeg_F',30});

% initial amounts set in modelObj.Species(k).InitialAmount. 
% previousdata, if any, stored in prevData. 
simData = sbiosimulate(modelObj, configsetObj);
   
if isempty(time_vector) 
    x_ode = simData.Data;
    t_ode = simData.Time;
else
    t_ode = [time_vector; simData.Time+time_vector(end)];
    x_ode = [prevData;simData.Data]; 
end

varargout{1} = simData;

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