% txtl_protection.m 
% VS, Oct 2012
%%
% Run negautoreg.m without protection first, followed by with protection.
% Plot and compare the results. 
%%
% Written by Vipul Singhal, Oct 2012
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

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'ptet(50)', 'rbs(20)', 'tetR(647)', 5, 'linear');
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 5*60*60)
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end

[t_ode,x_ode] = txtl_runsim(Mobj,configsetObj,[],[]);


% Top row: protein and RNA levels

% DNA and mRNA plot
dataGroups{1,1} = 'DNA and mRNA';
dataGroups{1,2} = {'ALL_DNA'};
dataGroups{1,3} = {'b-','r-','b--','r--','y-','c-','g-','g--'};

% Gene Expression Plot
dataGroups{2,1} = 'Gene Expression';
%dataGroups{2,2} = {'protein deGFP-lva-terminator*'};
%dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'};

% Resource Plot
dataGroups{3,1} = 'Resource usage';

txtl_plot(t_ode,x_ode,Mobj,dataGroups);

% figure(1); clf(); subplot(2,1,1);
% iTetR = findspecies(Mobj, 'protein tetR');
% iGamS = findspecies(Mobj, 'protein gamS');
% iGFP = findspecies(Mobj, 'protein deGFP');
% iGFPs = findspecies(Mobj, 'protein deGFP*');
% plot(t_ode/60, x_ode(:, iTetR), 'b-', t_ode/60, x_ode(:, iGamS), 'r-', ...
%   t_ode/60, x_ode(:, iGFP) + x_ode(:, iGFPs), 'g--', ...
%   t_ode/60, x_ode(:, iGFPs), 'g-');
% 
% title('Gene Expression');
% lgh = legend({'TetR', 'GamS', 'GFPt', 'GFP*'}, 'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [nM]');
% xlabel('Time [min]');
% 
% % Second row, left: resource limits
% subplot(2,2,3);
% iNTP = findspecies(Mobj, 'NTP');
% iAA  = findspecies(Mobj, 'AA');
% iRNAP  = findspecies(Mobj, 'RNAP70');
% iRibo  = findspecies(Mobj, 'Ribo');
% mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
% plot(...
%   t_ode/60, x_ode(:, iAA)/x_ode(1, iAA), 'b-', ...
%   t_ode/60, x_ode(:, iNTP)/x_ode(1, iNTP), 'r-', ...
%   t_ode/60, x_ode(:, iRNAP)/max(x_ode(:, iRNAP)), 'b--', ...
%   t_ode/60, x_ode(:, iRibo)/x_ode(1, iRibo), 'r--');
% 
% title('Resource usage');
% lgh = legend(...
%   {'NTP [mM]', 'AA [mM]', 'RNAP70 [nM]', 'Ribo [nM]'}, ...
%   'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [normalized]');
% xlabel('Time [min]');
% 
% % Second row, right: DNA and mRNA
% subplot(2,2,4);
% iDNA_tetR = findspecies(Mobj, 'DNA ptet--rbs--tetR');
% iDNA_gamS = findspecies(Mobj, 'DNA p70--rbs--gamS');
% iRNA_tetR = findspecies(Mobj, 'RNA rbs--tetR');
% iRNA_gamS = findspecies(Mobj, 'RNA rbs--gamS');
% plot(t_ode/60, x_ode(:, iDNA_tetR), 'b-', ...
%   t_ode/60, x_ode(:, iDNA_gamS), 'r-', ...
%   t_ode/60, x_ode(:, iRNA_tetR), 'b--', ...
%   t_ode/60, x_ode(:, iRNA_gamS), 'r--');
% 
% title('DNA and mRNA');
% lgh = legend(...
%   names([iDNA_tetR, iDNA_gamS, iRNA_tetR, iRNA_gamS]), ...
%   'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [nM]');
% xlabel('Time [min]');

%% Run the simulation with protection
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E6');
tube2 = txtl_buffer('e1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'thio-junk(500)-ptet(50)', 'rbs(20)', 'tetR(647)', 5, 'linear');
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', 5, 'linear');
dna_gamS = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'gamS(1000)', 1, 'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3], [6, 2, 2]);

% Run a simulation
configsetObj = getconfigset(Mobj, 'active');
set(configsetObj, 'StopTime', 5*60*60)
if ~strcmp(version('-release'),'2012a')
 set(configsetObj, 'SolverType', 'ode23s');
end
[t_ode,x_ode] = txtl_runsim(Mobj,configsetObj,[],[]);


txtl_plot(t_ode,x_ode,Mobj,dataGroups);

% % Top row: protein and RNA levels
% figure(2); clf(); subplot(2,1,1);
% iTetR = findspecies(Mobj, 'protein tetR');
% iGamS = findspecies(Mobj, 'protein gamS');
% iGFP = findspecies(Mobj, 'protein deGFP');
% iGFPs = findspecies(Mobj, 'protein deGFP*');
% plot(t_ode/60, x_ode(:, iTetR), 'b-', t_ode/60, x_ode(:, iGamS), 'r-', ...
%   t_ode/60, x_ode(:, iGFP) + x_ode(:, iGFPs), 'g--', ...
%   t_ode/60, x_ode(:, iGFPs), 'g-');
% 
% title('Gene Expression');
% lgh = legend({'TetR', 'GamS', 'GFPt', 'GFP*'}, 'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [nM]');
% xlabel('Time [min]');
% 
% % Second row, left: resource limits
% subplot(2,2,3);
% iNTP = findspecies(Mobj, 'NTP');
% iAA  = findspecies(Mobj, 'AA');
% iRNAP  = findspecies(Mobj, 'RNAP70');
% iRibo  = findspecies(Mobj, 'Ribo');
% mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
% plot(...
%   t_ode/60, x_ode(:, iAA)/x_ode(1, iAA), 'b-', ...
%   t_ode/60, x_ode(:, iNTP)/x_ode(1, iNTP), 'r-', ...
%   t_ode/60, x_ode(:, iRNAP)/max(x_ode(:, iRNAP)), 'b--', ...
%   t_ode/60, x_ode(:, iRibo)/x_ode(1, iRibo), 'r--');
% 
% title('Resource usage');
% lgh = legend(...
%   {'NTP [mM]', 'AA [mM]', 'RNAP70 [nM]', 'Ribo [nM]'}, ...
%   'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [normalized]');
% xlabel('Time [min]');
% 
% % Second row, right: DNA and mRNA
% subplot(2,2,4);
% iDNA_tetR = findspecies(Mobj, 'DNA thio-junk-ptet--rbs--tetR');
% iDNA_gamS = findspecies(Mobj, 'DNA p70--rbs--gamS');
% iRNA_tetR = findspecies(Mobj, 'RNA rbs--tetR');
% iRNA_gamS = findspecies(Mobj, 'RNA rbs--gamS');
% plot(t_ode/60, x_ode(:, iDNA_tetR), 'b-', ...
%   t_ode/60, x_ode(:, iDNA_gamS), 'r-', ...
%   t_ode/60, x_ode(:, iRNA_tetR), 'b--', ...
%   t_ode/60, x_ode(:, iRNA_gamS), 'r--');
% 
% title('DNA and mRNA');
% lgh = legend(...
%   names([iDNA_tetR, iDNA_gamS, iRNA_tetR, iRNA_gamS]), ...
%   'Location', 'Northeast');
% legend(lgh, 'boxoff');
% ylabel('Species amounts [nM]');
% xlabel('Time [min]');
% 
% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
