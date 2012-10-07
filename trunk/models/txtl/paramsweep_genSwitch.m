% paramsweep_genSwitch.m
% Run the toggle switch function, gen_switch(), for different values of
% initial protein concentrations and plot the concentrations of one
% protein against the other. Show bistability and basins of convergence. 

% Written by Vipul Singhal, Sep 2012
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
%

%clean up
%close all
%clear all
clc

% define initial protein concentrations
lacI = [0 4]; 
% I do not yet know why an initial concentration of x ends up as x/5 when the simulation is run. 
% for instance 4 will come up at 0.8, 16 as 3.2.
% In 'gen_switch', right up until sbiosimulate, this initial concentration
% is correct. But something in sbiosimulate divides it by 5. 
tetR = [0 4];
activeInducer = {'both', 'IPTG', 'aTc'}; % aTc induces lacI, IPTG induces tetR
aTc_conc = 5;
IPTG_conc = 5;
timevects_tetR = cell(1,length(lacI)*length(tetR));
timevects_lacI = cell(1,length(lacI)*length(tetR));
xvects_tetR = cell(1,length(lacI)*length(tetR));
xvects_lacI = cell(1,length(lacI)*length(tetR));
% Run the simulation multiple times and store the results. 
for k = 1:length(activeInducer)
    for i = 1:length(lacI)
        for j = 1:length(tetR)
            
            [Mobj, t_ode, x_ode, names] = gen_switch(tetR(j),lacI(i), 60*60, activeInducer{k}, aTc_conc, IPTG_conc);
            ilacI = findspecies(Mobj, 'protein lacI');
            iTetR = findspecies(Mobj, 'protein tetR');
            timevects_tetR{1,(k-1)*length(lacI)*length(tetR) + (i-1)*length(tetR)+j} = t_ode/60;
            xvects_tetR{1,(k-1)*length(lacI)*length(tetR) + (i-1)*length(tetR)+j} = x_ode(:,iTetR);        
            timevects_lacI{1,(k-1)*length(lacI)*length(tetR) + (i-1)*length(tetR)+j} = t_ode/60;
            xvects_lacI{1,(k-1)*length(lacI)*length(tetR) + (i-1)*length(tetR)+j} = x_ode(:,ilacI);
        end
    end
end


close all

hold off
colorlist={'r.','r--','g.','g--','b.','b--','k.','k--','y-','y--','m-','m--','c-','c--'}; 
cllen=length(colorlist);

h = figure;

for k = 1:length(activeInducer)
    subplot(length(activeInducer),1,k);
    M=1;
    for i = 1:length(lacI)
        for j = 1:length(tetR)
            if M>cllen M=1; end;
            hold on
            plot(timevects_tetR{1,(k-1)*length(lacI)*length(tetR) +(i-1)*length(tetR)+j}, xvects_tetR{1,(k-1)*length(lacI)*length(tetR) +(i-1)*length(tetR)+j},colorlist{M}, timevects_lacI{1,(k-1)*length(lacI)*length(tetR) +(i-1)*length(tetR)+j}, xvects_lacI{1,(k-1)*length(lacI)*length(tetR) +(i-1)*length(tetR)+j},colorlist{M+1})
           legends{1,(i-1)*2*length(tetR)+2*j} = ['lacI; INIT: L=' num2str(lacI(i)) ',T=' num2str(tetR(j))];
           legends{1,(i-1)*2*length(tetR)+2*j-1} = ['tetR; INIT: L=' num2str(lacI(i)) ',T=' num2str(tetR(j))];
           M=M+2;
        end
    end
    legend(gca, legends, 'Location','NorthEastOutside')
    xlabel('time/min'); ylabel('lacI and tetR conc')
    title([activeInducer{k} ' as the inducer'])
    clear legends   
end


h2 = figure;
hold on
M = 1;
for k = 1:length(activeInducer)
    for i = 1:length(lacI)
        for j = 1:length(tetR)
             if M>cllen M=1; end;
            plot(xvects_tetR{1,(k-1)*length(lacI)*length(tetR)+(i-1)*length(tetR)+j}, xvects_lacI{1,(k-1)*length(lacI)*length(tetR)+(i-1)*length(tetR)+j},colorlist{M})
            legends{1,(k-1)*length(lacI)*length(tetR) +(i-1)*length(tetR)+j} = ['lacI = ' num2str(lacI(i)) ', tetR = ' num2str(tetR(j)) 'ind: ' activeInducer{k}];
            M = M+1;
        end
    end
end

xlabel('tetR'); ylabel('lacI')
legend(gca, legends, 'Location','NorthEastOutside')


