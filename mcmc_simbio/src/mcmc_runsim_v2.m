function mi =  mcmc_runsim_v2(tstamp, projdir, data_info, mcmc_info, varargin)
% version 2 of the runsim file. 
% This version is for the multi modal version of the mcmc problem. 
% 
%
% mcmc_runsim: run the mcmc estimation using the affine invariant 
% ensemble sampler. 
% 
% INPUTS: 
%
% tstamp:       The time stamp string of the simulation generated by 
%               project_init. Format is 'yyyymmdd_HHMMSS'
% 
% projdir:      Directory where the generated data subdirectory
%               (simdata_yyyymmdd_HHMMSS, where yyyymmdd_HHMMSS is 
%               the time stamp) will be created. 
%
% tv:           time vector in the units of seconds. 
%
% da:           A matlab array containing experimental data. This array 
%               has dimensions nTimePoints x nMS x nReplicates x nDoses
%               where 
%               -nTimePoints:    length(tv), i.e., the number of time points. 
%               -nMS:            The number of measured species. Corresponds 
%                               to values given in the mcmc_info and data_info 
%                               structs. 
%               -nReplicates
%               -nDoses:         Number of dose combinations (initial conditions)
% 
% mobj:         Simbiology model object. 
%
% mi:           mcmc_info struct. Type 'help mcmc_info_dsg2014_mrna' or 
%               'help mcmc_info_template' into the MATLAB command line to learn 
%               more. 
%               
% Optional name-value pair arguments:
% InitialDistribution:  Initial distribution of walker points. This can be 
%                       Latin hyprcube sampled (Value: 'LHS'), gaussian 
%                       distributed (Value: 'gaussian') about the midpoint of 
%                       mi.paramranges or uniformly distributed (Value: 'unifrand'). 
% 
% Width:                Applies to the width of the gaussian or uniform random
%                       parameter distribution around the midpoint given by 
%                       mi.paramranges. 
%
% UserInitialize:       The user provides a matrix of initial walker positions. 
%                       When this input is specified, 'InitialDistribution' and 
%                       'Width' are ignored. 
%
% FitOption:            Allows for fitting to data in 3 modes:
%               'FitMedian':    This mode computes the curvewise median (Default).
%                               of the data over the replicates, and fits the model 
%                               to this. 
%               'FitMean':      Compute the mean of the replicates, and fit to this 
%                               mean
%               'FitAll'        Fit all the curves. 

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

p = inputParser;
p.addParameter('InitialDistribution', 'LHS', @ischar); % LHS, gaussian, unifrand
p.addParameter('Width', 0.1, @isnumeric)
p.addParameter('UserInitialize', [], @isnumeric)
p.addParameter('FitOption', 'FitMedian', @ischar); % 'FitMean', 'FitAll'
p.addParameter('pausemode', false, @islogical)
p.addParameter('DoseNormalization', false, @islogical)

p.addParameter('multiplier', 1, @isnumeric); % every 10 iterations, 
% the stepsize is multiplied by this number. Set it to something like 2 or 4 
% in the hope that the walkers
% will be able to explore more of the space during this iteration. 
p.parse(varargin{:});
p = p.Results;

% get the three objects. 
ri = mcmc_info.runsim_info;
mi = mcmc_info.model_info;
mai = mcmc_info.master_info;

[da, mi, tv] = exportmobj(mi, data_info, p.FitOption);
di = data_info;
% 
% nTopo = length(mi);
% nGeom = zeros(length(mi),1);
% for i = 1:nTopo
%     nGeom(i) = length(mi(i).dataToMapTo);
% end
% 
% % V2 data
% % for each topology geometry pair, compute the data to fit
% % ASSUME that the data dimensions across geometries is the same
% % only differs across topologies at most. 
% % Despite this, we allow for different geometries within a topology 
% % to point to different data info array elements, just that all of 
% % these elements must have the same number of timepoints, measured species, 
% % replicates and dosing combinations. (version 3 of this code can be even 
% % more general, with each topo-geom pair getting its own cell. this will be 
% % slower.)
% 
% da = cell(nTopo, 1);
% 
% for i = 1:nTopo % each topology
%     % each of the nGeom(i) geometries has a data info array element it points to. 
%     % since the dimensions of these is assumed to be equal, we just use the first
%     % one to set the empty array: 
%     data_info_element = mi(i).dataToMapTo(1);
%     currda = data_info(data_info_element).dataArray;
%     % Transform Experimental - compute mean or median or nothing
%     currda = computeFitOption(currda,  p.FitOption);
%     da{i} = currda;
%     tv{i} = data_info(data_info_element).timeVector;
%     for j = 2:nGeom(i) % each geometry
%         data_info_element = mi(i).dataToMapTo(j);
%         currda = data_info(data_info_element).dataArray;
%         % Transform Experimental - compute mean or median or nothing
%         currda = computeFitOption(currda,  p.FitOption);
%         % concatenate in the 5th dimension (the geometries dimension.)
%         da{i} = cat(5, da{i}, currda); 
%     end
%     % EXPORT MODEL object to get it ready for MCMC
%     % the resulting object is of class SimBiology.export.Model
%     % documentation: 
%     % https://www.mathworks.com/help/simbio/ref/simbiology.export.model-class.html
%     % Sven Mesecke's blog post on using the exported model class for 
%     % parameter inference applicaton. 
%     % http://sveme.org/how-to-use-global-optimization-toolbox-algorithms-for-
%     % simbiology-parameter-estimation-in-parallel-part-i.html
% 
%     mobj = mi(i).modelObj;
% 
%     enuo = mi(i).namesUnord;% estimated names unordered
%     
%     ep = sbioselect(mobj, 'Type', 'parameter', 'Name', ...
%         mi(i).namesUnord);% est parameters
% 
%     es = sbioselect(mobj, 'Type', 'species', 'Name', ...
%         mi(i).namesUnord);% est species
% 
%     aps = [ep; es]; % active parameters and species
% 
%     % reorder the parameter and species so they are in the same order as that
%     % in the model. 
%     eno = cell(length(aps), 1);% est names ordered
%     ds = sbioselect(mobj, 'Type', 'species', 'Name', mi(i).dosedNames);
%     emo{i} = export(mobj, [ep; es; ds]); % exported model object, dosed species names. 
%     SI = emo{i}.SimulationOptions;
% 
%     % each of the nGeom(i) geometries has a data info array element it points to. 
%     % since the dimensions of these is assumed to be equal, we just use the first
%     % one to set the empty array: 
%     data_info_element = mi(i).dataToMapTo(1);
%     SI.StopTime = data_info(data_info_element).timeVector(end);
%     accelerate(emo{i});
%     
%     mi(i).emo = emo{i}; % exported model object. 
%     orderingIx = zeros(length(aps),1);
%     orderingIx2 = orderingIx;
%     for k = 1:length(aps)
%         eno{k} = aps(k).Name;
%         for kk = 1:length(enuo)
%             if strcmp(eno{k}, enuo{kk} )
%                 orderingIx(k) = kk; % eno = enuo(orderingIx);
%                 % the kth element of orderingIx is kk. so the kth element of 
%                 % enuo(orderingIx) is enuo(kk). But this is just eno(k). And eno 
%                 % has the property of the kth element being eno(k). (as seen 
%                 % from "if eno{k} == enuo{kk} ")
% 
%                 orderingIx2(kk) = k; %i.e., enuo = eno(orderingIx2); 
%                 % the kkth element of orderingIx2 is k. so the kk th element of
%                 % eno(orderingIx2) is eno(k). But the vector with this property is 
%                 % simply enuo. (as seen from "if eno{k} == enuo{kk} ")
%             end
%         end
%     end
% 
%     mi(i).orderingIx = orderingIx; % these two arrays will be VERY useful. 
%     mi(i).orderingIx2 = orderingIx2; % this one being the second. 
%     mi(i).namesOrd = eno; % est names ordered. 
% end

% V2 model export - export with all parameters, and set fixed and estimated 
% parameters per iteration of mcmc. One exported model for each topology
% make sure the reordering of the parameters when exporting is carefully 
% taken care of. 

%% COMPUTE INITIAL WALKER POSITIONS
% a very cool idea: if a parameter with the same semantic meaning 
% appears as individual parameters across different topologies and 
% geometries, then we expect that generally the final estimated values 
% should be close between the different verions of the parameter 
% (the number of RNAP in extract 1, extract 2 and so on should be similar
% to an order of magnitude, even if they are differet.) 
% To be clear, this is not saying that parameters that get shared across 
% topologies and geometries need to be close. Those are EQUAL BY DEFINITION. 
% All it is saying is that within a master vector if a parameter is semantically
% similar to another, they should have the same STARTING values. 

if isempty(p.UserInitialize)
    minit = integrableLHS_v2(mi, mai, ri, ...
        'distribution', p.InitialDistribution, ...
        'width', p.Width);

else
    minit = p.UserInitialize;
    % assume all the user defined points are integrable.
end

%% SETUP FUNCTIONS
% setup the log prior, log likelihood function and lognormvec functions
lognormvec=@(res,sig) -(res./sig).^2 -log(sqrt(2*pi)).*sig;

logprior = @(logp) all(mai.paramRanges(:, 1) < logp) &&...
    all(logp < mai.paramRanges(:,2));

sigg = ri.stdev/ri.tightening;


% need to transform the data array to summary stats before 
% sending it into the gen residuals function. 
mv = mai.masterVector;
estParamIx = setdiff((1:length(mv))', mai.fixedParams);

loglike = @(logp) gen_residuals_v2(logp, estParamIx, ...
                        mv, da,...
                        tv, mi, lognormvec, sigg);

% BURN IN: run the burn in simulation
if isempty(p.UserInitialize)
    tic
    [m] =gwmcmc_vse(minit,{logprior loglike},...
        ri.nPoints,...
        'StepSize',ri.stepSize , ...
        'ThinChain',ri.thinning,...
        'Parallel', ri.parallel);
    toc
    
    minit = m(:,:,end);
    clear m
else
    disp('User initialized intitial walker positions, skipping burn in phase')
end

% specify where to save things 
cfname = cell(ri.nIter, 1);
specificproj = [projdir '/simdata_' tstamp];

%% we save useful variables in a one off manner here (ie, outside the loops)
fname = ['full_variable_set_' tstamp]; % filename
save([specificproj '/' fname]);

% run the actual simuations, saving the data every iteration
for i = 1:ri.nIter %
    if p.pausemode
        if ~mod(i, 1)
            fprintf('Pausing for 2 minutes before starting run number %d. \n', i);
            pause(120)        
        end
    end
    
    
    if  i == 3 || ~mod(i, 3)
        ssize = p.multiplier*ri.stepSize;
        fprintf('Mixup round! The step size for this iteration is set to \n %d * %d = %d.\n',...
            p.multiplier, ri.stepSize, ssize);
    else
        ssize = ri.stepSize;
        
        
    end
    
    
    tic
    fprintf('starting mcmc %d\n', i);
    [m] = gwmcmc_vse(minit,{logprior loglike},...
        ri.nPoints, ...
        'StepSize',ssize , ...
        'ThinChain',ri.thinning, 'Parallel', ri.parallel);
    
    fprintf('ending mcmc %d\n', i);
    toc
    fname = ['mcmc' tstamp '_ID' num2str(i)] ;
    cfname{i} = fname;
    save([specificproj '/' fname], 'm');
    % the only thing that is different in each run are the above
    %             pause(1)
    minit = m(:,:,end);% + 0.1*randn(size(m(:,:,end-1)));
    
    clear m
    
end

% generate log file
    if isempty(p.UserInitialize)
        initialization_used = p.InitialDistribution;
    else
        initialization_used = 'User_initialized';
    end
mcmc_log_v2(tstamp,projdir, specificproj, mcmc_info, data_info, initialization_used);
end



