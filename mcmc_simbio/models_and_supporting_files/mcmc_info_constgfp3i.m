function mcmc_info = mcmc_info_constgfp3i(modelObj)
% mcmc_info_constgfp3i.m 
% This file contains information on how to set up a constitutive GFP
% experiment using the Simbiology model generated by the function 
% model_protein3. It is the simplest possible example of the parameter inference
% problem, and does not utilize the concurrent parameter inference capabilities 
% of mcmc_simbio. It has a single variant (geometry) of a single model
% (topology), and a simple set of parameters to be estimated. It is to be
% used with the project files proj_mcmc_tutorial.m and
% proj_protein_constgfp3i.m. 
% 
% The model to be used is set up by the file model_protein3.m in the 
% models_and_supporting_files directory, and has equations 
% 
% % dG + pol <-> dG_pol  (kf, kr) 
% dG_pol -> dG + pol + pG (kc) 
% 
% where dG is the GFP protein DNA, pol is a species signifying the
% transcription and translation machinery, and pG is the GFP protein output
% of the model. 
% 
% In the experiment, the initial dG concentration is set to different
% values (10, 30 and 60). I.e., the dosedNames in mcmc_info.model_info
% has a single species 'dG' in it, and the dosedVals is the row vector 
% [10 30 60]. For each dose, the corresponding pG trajectories are 
% measured. In the parameter inference problem, these simulated model 
% trajectories are compared to artificial 'data' trajectories obtained from 
% simulating the model. 
% 
% We fix kf to 5, and estimate the rest of the parameters. 
% 
% Type help mcmc_info to read more about the mcmc_info struct array. 

%% Start here. 
% User readable description of the circuit. Will be used in the log file generated
% from the MCMC inference procedure.
circuitInfo = ...
    [' dG + pol <-> dG_pol  (kf, kr \n'... )
    'dG_pol -> dG + pol + ppG (kc)\n'...
    'single topology, single geometry.'];

rkfdG = 5; % rate of kfdG, in nM-1 . s-1
rkrdG = 300; % rate of krdG, in s-1
rkc = 0.012; %rate of kc, in s-1
cpol = 100; % concentration of pol, in nM

% names of the parameters and species in the model whose values are to be
% set by the values in the master vector. Since there is only one topology
% and geometry, the master vector also only has four entries. 
activeNames = ...
    {'kfdG'
    'krdG'
    'kcp'
    'pol'}; 

% The value of kfdG will be fixed, so only three params and species get
% estimated. Thus, the MCMC algorithm only searches a 3D space. 
estParams = {'krdG'
    'kcp'
    'pol'};

% The master vector is a vector of parameter values that is to be
% distributed to the model-experiment pairs. Here, since there is only one
% topology-geometry pair, it matches the activeNames array above. In
% general the master vector will collect the full set of parameters and
% species that are to be set across all the topologies and geometries (both
% as fixed and to-be-estimated values)
masterVector = log([rkfdG 
                    rkrdG
                    rkc
                    cpol]);

% The fixedParams vector sets the index of the parameters in the
% masterVector that stay constant withing the MCMC estimation problem. This
% value within the master vector must be set to the actual value it will
% need to be in all the simulations. (All other values within the master
% vector can be arbitrary because they will be set to various proposed values
% during the estimation procedure)
fixedParams = [1];

% The estParamsIx array contains the indices of the parameters that are to
% be estimated in the model. The line below will always be the same, since
% the fixedParams and the estParamsIx form a partition of the indices of
% masterVector. 
estParamsIx = setdiff((1:length(masterVector))', fixedParams);

% The paramMap array here is simply a column vector of the same length as
% the masterVector. In general, for each unique model topology, there will
% be one paramMap. Please type help mcmc_info in the command line to learn
% more about specifying the paramMaps array. See, for example, 
% mcmc_info_constgfp3tetR1.m for an example of a multiple topology 
% estimation problem. 
paramMap = (1:length(masterVector))';

% The paramRanges array is an array of size length(estParamsIx) x 2, and
% contains the upper and lower bounds for each parameter to be estimated.
% We note that the parameters in masterVector are log transformed, and 
% therefore, so are the bounds. Here we pick a range of parameter values that are
% +- 3 away from the nominal parameter values in the masterVector. 
paramRanges =  [masterVector(estParamsIx)-3 masterVector(estParamsIx)+3];

% The dataIndices array specifies which dataset within the data_info struct
% each geometry within the given topology points to. Since we only have one
% topology-geometry pair, and we will be generating data specifically for
% this 
dataIndices = 1;

%% next we define the dosing strategy.
dosedNames = {'dG'}; % name string of the species to be dosed. 
dosedVals = [10 30 60]; % values that species is to take. 

measuredSpecies = {{'pG'}}; % species to be measured. Note the double cell.
msIx = 1; % column index of the data_info.dataArray array that 
% corresponds to the measured species. 


%% setup the MCMC simulation parameters

stdev = 10; % the standard deviation in the likelihood function. 

tightening = 1; % See mcmc_info for more info. 

nW =200; % number of walkers. good values: 200 - 400 ish

stepsize = 1.1; % MCMC step size. try: 1.1 to 4 ish. DO NOT USE 1. 

niter = 10; % actual: 2 - 50 ish. Number of times to loop the MCMC. See mcmc_info

npoints = 4e4; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of
%                        params is small)

thinning = 10; % good values: 10 to 40 ish. 
% Number of steps to skip before recording positions of the walkers. 

%% pull all this together into an output struct.

runsim_info = struct('stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepSize', {stepsize}, ...
    'nIter', {niter}, ...
    'nPoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', false);

model_info = struct(...
    'circuitInfo',{circuitInfo},...
    'modelObj', {modelObj},... % array of model objects (different topologies)
    'modelName', {modelObj.name},...; % model names.
    'namesUnord', {activeNames}, ... % names of parameters per model, unordered.
    'paramMaps', {paramMap}, ... % paramMap: matrix mapping elements in the 
    ...                   % master vector to the parameters or species given 
    ...                   % by active names for a given topology. 
    ...                   % Type help mcmc_info for more information.
    'dosedNames', {dosedNames},... % cell arrays of species. cell array corresponds
    ...                               % to a model.
    'dosedVals', {dosedVals},...  % matrices of dose vals
    'measuredSpecies', {measuredSpecies}, ... % cell array of cell arrays of
    ...                  % species names. the elements of the inner
    ...                  % cell array get summed.
    'measuredSpeciesIndex', {msIx},...
    'dataToMapTo', dataIndices); % each dataToMapTo property within an 
% element of the model_info array is a vector of length # of geometries.

%% ignore this for now. It has a subtle use, and we will update the
% documentation to describe how this is used in a future release. 
semanticGroups = num2cell((1:length(estParams))'); 
%arrayfun(@num2str, 1:10, 'UniformOutput', false);

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