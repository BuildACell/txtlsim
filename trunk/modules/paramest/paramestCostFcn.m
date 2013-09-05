% Written by Zoltan A Tuza, Mar 2013
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


function cost = paramestCostFcn(p_act,data)

% sampling time points
tspan = data.xdata;
% mapping external parameters
pVector = data.p0;

% deciding initial values are estimated or not
if size(data.x0_mapping,1) == 0
    pVector(data.p_mapping) = p_act;
else
    pVector(data.p_mapping) = p_act(1:size(data.p_mapping,1));
    
    % updating initial values with current estimate
    cases = size(data.x0,2);
    startId = size(data.p_mapping,1);
    initParams = size(data.x0_mapping,1);
    for k=1:cases
        data.x0(data.x0_mapping,k) = p_act(startId+k:(startId+k)+initParams-1);
    end
    
end

% evaluate parameters dependencies
if ~isempty(data.p_dependecies)
    cellfun(@(x) evalc(strrep(x,'p','pVector')),data.p_dependecies,'UniformOutput',false);
end

simTime =[];
%% start of simulations
try
    % run simulations with different initial values 
    for k=1:size(data.x0,2)
        tic
        [t(:,k),y{k}] = ode15s(@(t,x)data.modelFcn(t,x,pVector),tspan,data.x0(:,k));
        simTime(k) = toc;
        
        if sum(y{k}(:,data.outputCheck)) == 0
            
            cost = Inf;
            disp(sprintf('cost:%f\n',cost));
            return;
            
        end
    end
    
%% calculate cost
    cost = 0;
    
    
    for k=1:size(data.x0,2)
        cost = cost + sum((data.ydata(:,k)-sum(y{k}(:,data.targetSpecies),2)*data.speciesScaleFactor).^2);
    end
    
    
    
catch err
    cost = Inf;
    err.message
end

str1 = sprintf('costV:%6.2G\t',cost);
str2 = strtrim(sprintf('p:%f ',p_act));
str3 = sprintf('%6.2f sec',sum(simTime));
disp([str1 ' ' str2 ' ' str3]);

if data.debugMode ==1 && cost < Inf
    
    figure(1);
    hold on
    for k=1:size(data.x0,2)
        plot(data.xdata/60,[data.ydata(:,k) sum(y{k}(:,data.targetSpecies),2)*data.speciesScaleFactor])
    end
    hold off
    title(sprintf('current cost value:%f',cost))
    disp('push any key to continue')
    pause
    close all
    
    
end