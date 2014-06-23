% txtl_protection.m 
% VS, Jun 2014
% Testing linear DNA protection via gamS
% 
%
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
close all 
clear all
% Initialization
tube1 = txtl_extract('E31VNPRL');
tube2 = txtl_buffer('E31VNPRL');
simulationTime =  10*60*60;
% datalocation = '/Users/vipulsinghal/Dropbox/Research/TXTL Toolbox work/2014 June/Characterization/Zach ACS paper';

%% Linear DNA expression without protection

lin_nogamS_data = [30.679	0.191
66.516	0.389
99.597	0.389
132.677	0.389
163.001	0.389
196.081	0.389
229.054	0.389
261.930	0.389
295.010	0.389
328.091	0.389
363.928	0.389
394.251	0.389
427.332	0.389
460.412	0.389];
tube3 = txtl_newtube('linearDNA_deg');
txtl_add_dna(tube3, 'p70(40)', 'rbs(20)', 'deGFP(750)', 16, 'linear'); %810 bp as in zsun ACS
txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 0, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);

[t_noGam,x_noGam] = txtl_runsim(Mobj,simulationTime);
iGFP_noGam = findspecies(Mobj, 'protein deGFP*');
%txtl_plot(t_noGam,x_noGam,Mobj);

%% Linear DNA expression with protection

lin_gamS_data = ...
    [30.679	0.807
63.760	2.188
97.045	2.978
130.233	3.570
165.865	4.344
196.189	4.936
229.269	5.134
259.593	5.726
295.430	5.718
328.511	6.113
361.591	6.113
394.671	6.310
427.860	6.113
460.940	6.508];
tube3 = txtl_newtube('linearDNA_gamS');
txtl_add_dna(tube3, 'p70(40)', 'rbs(20)', 'deGFP(750)', 16, 'linear'); %810 bp as in zsun ACS
dna_gamS = txtl_add_dna(tube3, 'p70(50)', 'rbs(20)', 'gamS(1000)', 0, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);
txtl_addspecies(Mobj,'protein gamS', 3500);

[t_Gam,x_Gam] = txtl_runsim(Mobj,simulationTime);
iGFP_Gam = findspecies(Mobj, 'protein deGFP*');
%txtl_plot(t_Gam,x_Gam,Mobj);

%% Plasmid DNA expression
plasmid_data = ...
[29.904	1.488
63.436	3.949
93.760	6.713
129.705	8.884
162.785	10.858
198.623	12.437
228.946	13.424
264.783	14.214
295.107	14.806
330.944	15.398
364.025	15.990
394.553	16.385
427.741	16.583
460.617	16.780];

tube3 = txtl_newtube('plasmid_geneexpr');
dna_deGFP = txtl_add_dna(tube3, 'p70(40)', 'rbs(20)', 'deGFP(750)', 16, 'plasmid');
Mobj = txtl_combine([tube1, tube2, tube3]);

[t_plasmid,x_plasmid] = txtl_runsim(Mobj,simulationTime);
iGFP_plasmid = findspecies(Mobj, 'protein deGFP*');
%txtl_plot(t_plasmid,x_plasmid,Mobj);

%% compilation of the above three experiments into Fig2a of zsun's ACS paper. 

figure
p = plot(t_noGam/60, x_noGam(:,iGFP_noGam), 'b',...
    lin_nogamS_data(:,1), lin_nogamS_data(:,2)*1000, 'b*',...
    t_Gam/60, x_Gam(:,iGFP_Gam), 'r',...
    lin_gamS_data(:,1), lin_gamS_data(:,2)*1000, 'r*',...
    t_plasmid/60, x_plasmid(:,iGFP_plasmid), 'k',...
    plasmid_data(:,1), plasmid_data(:,2)*1000, 'k*');
set(p, 'LineWidth', 2)
legend(p, 'linear sim', 'linear DNA Exp',...
    'linear gamS@3.5uM sim', 'lin gamS@3.5uM exp',...
    'plasmid sim', 'plasmid exp', 'Location', 'NorthEastOutside')
title('Constitutive production: deGFP DNA @ 16nM','fontsize', 14)
xlabel('time/min','fontsize', 14)
ylabel('deGFP conc/nM','fontsize', 14)
set(gca,'fontsize', 14)

% 
% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
