% geneexpr.m - basic gene expression reaction
% R. M. Murray, 9 Sep 2012
% 
% This file contains a simple example of setting up a TXTL simulation
% for gene expression using the standard TXTL control plasmid.
% 
close all 
clear all

% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E30_1');
tube2 = txtl_buffer('E30_1');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('ptet_charac');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'tetR(1000)', ...	% promoter, rbs, gene
   1, ...					% concentration (nM)
  'plasmid');					% type
dna_deGFP = txtl_add_dna(tube3, ...
  'ptet(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   2, ...					% concentration (nM)
  'plasmid');

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);
   txtl_addspecies(Mobj, 'aTc', 100);

cs = getconfigset(Mobj);
set(cs.RuntimeOptions, 'StatesToLog', 'all');
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

 txtl_plot(simData,Mobj);
 
%   print('-djpeg','-r100',['case' num2str(caseNum) 'main.jpeg'])
%             saveas(gcf, ['case' num2str(caseNum) 'main.fig'])
% figure
% plot(t_ode/60, x_ode(:,1)+x_ode(:,5)+x_ode(:,19))
% title('RNAP+RNAP70+RNAP28')

%    1         contents        RNAP                                    100               
%    2         contents        protein sigma70                         35                
%    3         contents        protein sigma28                         20                
%    4         contents        Ribo                                    30                
%    5         contents        RNAP70                                  0                 
%    6         contents        RNase                                   100               
%    7         contents        AGTP                                    3.18005e+06       
%    8         contents        CUTP                                    1.90803e+06       
%    9         contents        AA                                      3.18005e+07       
%    10        contents        protein deGFP                           0                 
%    11        contents        protein deGFP*                          0                 
%    12        contents        RNA rbs--deGFP                          0                 
%    13        contents        Ribo:RNA rbs--deGFP                     0                 
%    14        contents        DNA p70--rbs--deGFP                     30                
%    15        contents        RNAP70:DNA p70--rbs--deGFP              0                 
%    16        contents        CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP    0                 
%    17        contents        term_RNAP70:DNA p70--rbs--deGFP         0                 
%    18        contents        AA:AGTP:Ribo:RNA rbs--deGFP             0                 
%    19        contents        RNAP28                                  0                 
%    20        contents        AGTP:RNAP70:DNA p70--rbs--deGFP         0                 
%    21        contents        CUTP:RNAP70:DNA p70--rbs--deGFP         0                 
%    22        contents        RNA rbs--deGFP:RNase                    0                 
%    23        contents        AA:AGTP:Ribo:RNA rbs--deGFP:RNase       0                 
%    24        contents        Ribo:RNA rbs--deGFP:RNase               0                 
%    25        contents        AGTP_UNUSE                              0   
% cellOfSpecies1 = { 'AGTP','CUTP','AGTP_UNUSE'
%     'DNA p70--rbs--deGFP', 'RNAP70:DNA p70--rbs--deGFP','term_RNAP70:DNA p70--rbs--deGFP'
%     'AGTP:RNAP70:DNA p70--rbs--deGFP','CUTP:RNAP70:DNA p70--rbs--deGFP','CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP'
%     'RNAP70', 'RNAP28', 'RNAP'  };
%              
%  plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies1)  
%  
% %   print('-djpeg','-r100',['case' num2str(caseNum) 'RNA.jpeg'])
% %             saveas(gcf, ['case' num2str(caseNum) 'RNA.fig'])
%             
% cellOfSpecies2 = {'Ribo','AA','RNase'
%                  'RNA rbs--deGFP','Ribo:RNA rbs--deGFP','AA:AGTP:Ribo:RNA rbs--deGFP'
%                  'RNA rbs--deGFP:RNase','Ribo:RNA rbs--deGFP:RNase', 'AA:AGTP:Ribo:RNA rbs--deGFP:RNase'};
% plotCustomSpecies2({Mobj}, {x_ode}, {t_ode}, cellOfSpecies2)
%             print('-djpeg','-r100',['case' num2str(caseNum) 'PROTEIN.jpeg'])
%             saveas(gcf, ['case' num2str(caseNum) 'PROTEIN.fig'])
           

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
