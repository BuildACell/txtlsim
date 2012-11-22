
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


function [t_ode, x_ode, varargout] = txtl_runsim(varargin)

    if nargin ~= 4 && nargin~= 5
        error('txtl_runsim should be called either with 4 or 5 arguments.');
    else
        modelObj = varargin{1};
        configsetObj = varargin{2};
        data = varargin{4};
        time_vector = varargin{3};
        if nargin == 5
            simData = varargin{5};
        end
        
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

if (isempty(data))
    %this is the first run, therefore we have setup the parameters
    txtl_setup_parameters(modelObj);
end
 
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