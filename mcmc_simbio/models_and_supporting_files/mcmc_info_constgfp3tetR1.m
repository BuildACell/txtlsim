function mcmc_info = mcmc_info_constgfp3tetR1(constmobj, tetRmobj)
% mcmc_info_constgfp3tetR1 - set up concurrent parameter inference for a
% set of common parameters informed by a constitutive gene expression
% circuit and a tetR repression circuit. More information on terminology
% and this file can be found by typing help mcmc_info in the command
% window. Before going through this file, we encourage you to read through
% and work through the tutorials on the constitutive gene expression: 
%
% proj_mcmc_tutorial.m
% proj_mcmc_tutorial_II.m
%
% (and optionally also look at the very similar files:
% proj_protein_constgfp3i.m
% proj_protein_constgfp3ii.m)
% 
% and most importantly the associated mcmc_info files 
%
% mcmc_info_constgfp3i.m
% mcmc_info_constgfp3ii.m
%
% to get a good idea of how mcmc_info files are set up, and set some
% context for the complexity of this file. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example demonstrates a slightly more complex example of the
% concurrence feature of the mcmc_simbio Bayesian Parameter inference
% toolbox. Here, we have two circuits: the constitutive gene expression
% circuit with model 
% 
% D + pol <-> D__pol  (kfdG, krdG)
% D__pol -> D + pol + protien (kcp)
%
% and the tetR repression circuit with model
%
% D_T + P <-> D_T:P -> D_T + P + T (kfdT, krdT; kcp)
% D_G + P <-> D_G:P -> D_G + P + G (kfdG, krdG; kcp)
% 2 T <-> T2 (kfdimTet, krdimTet)
% D_G + T2 <-> D_G:T2 (kfseqTet, krseqTet)
%
% Here each model is a different topology (thus there is only one
% geometry associated with each topology). The paramters that are shared 
% between the two topologies are: kfdG, krdG, kcp, pol. The remaining
% parameters are specific to the topology-geometry pair they appear in (in
% this case all the remaining parameters appear in the TetR repression
% circuit topology. Furthermore, we will set the forward rate parameters in
% all the reversible reaction to be fixed parameters, and therefore only
% estimate the reverse rate parameters. 
% 
% This file is used in proj_constgfp3tetR1.m. The data used by this file is
% artificially generated using the data_artificial_v2 function. 
% 

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% User readable description of the circuit. Will be used in the log file generated
% from the MCMC inference procedure.
circuitInfo1 = ...
    [' D + pol <-> D__pol  (k_f, k_r) \n'... 
    'D__pol -> D + pol + protien (kc)\n'];

circuitInfo2 = ...
    ['D_T + P <-> D_T:P -> D_T + P + T (kfdT, krdT, kcpT) \n'...
    'D_G + P <-> D_G:P -> D_G + P + G (kfdG, krdG, kcpG) \n'...
    '2 T <-> T2 (kfdim, krdim)\n'...
    'D_G + T2 <-> D_G:T2 (kfseq, krseq) \n'];

cpol = 100; % nM
rkfdG = 5; % nM-1s-1
rkrdG = 300; % s-1
rkfdT = 5;
rkrdT = 300;
rkcp = 0.012; %s-1
rkfdimTet = 20; % nM-1s-1
rkrdimTet = 10; % s-1
rkfseqTet = 20; % nM-1s-1
rkrseqTet = 10; % s-1

% Each topology has its own set of active names. These are the parameters
% and species that are set from values taken from the master vector (using
% indices drawn from paramMaps). 
activeNames1 = {'kfdG'
    'krdG'
    'kcp'
    'pol'};
activeNames2 = ...
    {'kfdG'
    'krdG'
    'kfdT'
    'krdT'
    'kfdimTet'
    'krdimTet'
    'kfseqTet'
    'krseqTet'
    'kcp'
    'pol'};

% This is the masterVector, each unique parameter only appears here once.
% The forward rate parameters (rkfdG, rkfdT, rkfdimTet, rkfseqTet) are not
% estimated, and their values are fixed to the values above. The remaining
% parameters are estimated during the MCMC procedure. 
masterVector = log([...
rkfdG 
rkrdG
rkfdT
rkrdT
rkfdimTet
rkrdimTet
rkfseqTet
rkrseqTet
rkcp
cpol]);

% fixedParams vector
fixedParams = [1 3 5 7];

% The indices of the elements of the master vector that are to be
% estimated. There are a total of six unique species and parameters that 
% are to be estimated. Thus, the MCMC algorithm explores a 6 dimensional
% space to perform the fitting. 
estParamsIx = setdiff((1:length(masterVector))', fixedParams);

% The names of the estiamted parameters and species in this case are: 
% {'krdG'} {'krdT'} {'krdimTet'} {'krseqTet'} {'kcp'} {'pol'}
estParams = activeNames2(estParamsIx);

% There are two topologies, and each topology has its own paramMaps matrix
% (composed of a single column vector)
% The first topology uses the 1st, 2nd, 9th and 10th elements of the master
% vector, which map to 'kfdG',  'krdG',  'kcp', and 'pol' in that topology.
paramMap1 = [1 2 9 10]'; 

% The second topology uses all the elements of the master vector, and maps
% them to the species and parameters in activeNames2. 
paramMap2 = (1:length(masterVector))';

% mastervec(estParamsIx) is the subset of elements of the master vector
% that are being estimated. These values are log transformed (to allow them 
% to take values on the entire real line), and we set the range of values
% to be within a 6D hypercupe of width 6 (in log space) centered on the
% nominal values given above. 
paramRanges =  [masterVector(estParamsIx)-3 masterVector(estParamsIx)+3];

%% Set data indices to map the models to corresponding data in the data_info struct. 
% data indices tell us which data set to use for each topology - geometry pair
% from the data_info struct array.
dataIndices1 = 1;
dataIndices2 = 2;

%% next we define the dosing strategy.
% each of the two topologies has its own dosing strategy. 

% topology 1 only has the initial DNA (dG) concentration as the species to be
% dosed. 
dosedNames1 = {'dG'};
dosedVals1 = [10 30 60];

% topology 2 has two species being dosed, the initial tetR DNA (dT) 
% concentration and the initial GFP DNA concentration (dG). dT takes three
% possible values: 0.1, 2 and 8, while dG takes two values: 10 and 30
% (arbitrary units), leading to a total of 3 x 2 = 6 dose combinations
% (reflected as the 6 columns of dosedVals2). 

dosedNames2 = {'dT'; 'dG'}; % this needs to be a column array

dosedVals2 = [0.1  0.1  2   2   8   8  ; 
              10   30   10  30  10  30];

%% Set the measured species in the model and data. 
% Both the topologies have the same measured species: pG. 
measuredSpecies1 = {{'pG'}};
measuredSpecies2 = {{'pG'}};

% the msIx1 and msIx2 are used to map the measuredSpecies (above) to the
% corresponding columns of the data_info.dataArray data matrix. 
msIx1 = 1; %
msIx2 = 1;

%% setup the MCMC simulation parameters

% 
% >> (sum(di(1).dataArray(:)) + sum(di(2).dataArray(:)))/10
% 
% ans =
% 
%    1.6505e+04
% so overall we have 1.6 x 10^5 as the full signal magnitude (this is assuming 
% no dose or topology reweighting). stdev of 1e2 makes
% a factor of 1000 or 0.1% as the noise model. 
% lets try this. I tired 1x10^4, this led to extremely broad distributions
% and loose fits. That is fine for a start -- high 'temperature'. 


stdev = 1e2; % try this out. this could be great! or bad...
tightening = 1; % default is 1. Type in help mcmc_info for more information 
nW = 3200*16;% number of walkers. good values: 200 - 400
stepsize = 1.3; % MCMC step size. try: 1.1 to 4. DO NOT USE 1.
niter = 10; % try: 2 - 50. Number of times to loop the MCMC. "help mcmc_info"
npoints = 3200*32*10*10; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of

%                        params is small)
thinning = 10; % good values: 10 to 40 ish. 
% Number of steps to skip before recording positions of the walkers. 

%% pull all this together into an output struct.
% the mcmc info struct now is an array struct, the way struct should be used!

runsim_info = struct('stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepSize', {stepsize}, ...
    'nIter', {niter}, ...
    'nPoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', true);

% for now we simply make the model_info have just one model (topology).
% But the code will be written in a way such that multiple models can be used.

model_info = struct(...
    'circuitInfo',{circuitInfo1, circuitInfo2},...
    'modelObj', {constmobj, tetRmobj},... 
    ... % array of model objects (different topologies)
    'modelName', {constmobj.name, tetRmobj.name},...; % model names.
    'namesUnord', {activeNames1, activeNames2}, ... 
    ... % names of parameters per model, unordered.
    'paramMaps', {paramMap1, paramMap2}, ...% paramMap: matrix mapping elements in the 
    ...                   % master vector to the parameters or species given 
    ...                   % by active names for a given topology. 
    ...                   % Type help mcmc_info for more information. 
    'dosedNames', {dosedNames1, dosedNames2},... 
    ...% cell arrays of species. cell array corresponds
    ...                               % to a model.
    'dosedVals', {dosedVals1, dosedVals2},...  % matrices of dose vals
    'measuredSpecies', {measuredSpecies1, measuredSpecies2}, ... 
    ...% cell array of cell arrays of
    ...                  % species names. the elements of the inner
    ...                  % cell array get summed.
    'measuredSpeciesIndex', {msIx1, msIx2},...
    'dataToMapTo', {dataIndices1, dataIndices2}); 
...% each dataToMapTo property within an 
% element of the model_info array is a vector of length # of geometries.


%% IGNORE this for now. It has a subtle use, and we will update the
% documentation to describe how this is used in a future release. 
semanticGroups = num2cell((1:length(estParamsIx))'); 

%% master parameter vector, param ranges,
master_info = struct(...
    'estNames', {estParams},...
    'masterVector', {masterVector},...
    'paramRanges', {paramRanges},...
    'fixedParams', {fixedParams},...
    'semanticGroups', {semanticGroups});

mcmc_info = struct('runsim_info', runsim_info, ...
    'model_info', model_info,...
    'master_info', master_info);

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
end