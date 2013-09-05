%% load datasets

e15_pr_gfp_mg_0821 = wholeExpfileReader('ACS_paper_data/e15_pr_gfp_pr_gfp_mg_grad_0821_data.csv',';',2);
e15_pr_gfp_mg_0821 = processMGData(e15_pr_gfp_mg_0821);
e15_pr_gfp_mg_0821.valveNames = {'pos cont','neg cont','pr-gfp 1nM','pr-gfp 2nM','pr-gfp 3nM','pr-gfp 5nM','pr-gfp 10nM','pr-gfp 15nM'...
                                 'pr-gfp-s15-mg 1nm','pr-gfp-s15-mg 2nM','pr-gfp-s15-mg 3nM','pr-gfp-s15-mg 5nM','pr-gfp-s15-mg 10nM','pr-gfp-s15-mg 15nM'};
 

data.xdata = e15_pr_gfp_mg_0821.t_vec;
data.ydata = e15_pr_gfp_mg_0821.noBg(:,3:8,2)./2683.9;

%% initial values and parameters

geneexpr
close all;
 dataOut = parseGetEqOutput(Mobj);

data.x0 = dataOut.initialValues;
data.p0 = dataOut.parameters;

data.modelFcn = str2func(dataOut.modelFcn);
data.costFcn = @paramestCostFcn;

% building initial cases
%
% select a species (e.g. DNA or induder, which are subject to variations in the experiments)
speciesInd = findStringInAList(dataOut.speciesNames,'[DNA p70--rbs--deGFP]');
% specify the initial values for the selected species  
speciesInitialValues = [1 2 3 4]; % 1-4nM
% transform x0 vector to the right format (each column -> one different initial values vector)
data.x0 = repmat(data.x0,1,size(speciesInitialValues,2)); 
data.x0(speciesInd,:) = speciesInitialValues;

% which species are subject to initial value estimation
%data.x0_mapping = speciesInd;
data.x0_mapping = [];
% adjusting parameter vector size based on the initial values
if ~isempty(data.x0_mapping)
    data.p0 = [data.p0; data.x0(speciesInd,:)'];
end

% select parameters for estimation, list is given by parseGetEqOutput (parameterNames field)
ParametersToEstimate  = {'TXTL_TL_rate',...
                         'TXTL_UTR_RBS_F',...
                         'ATPdeg_F'
                         };
                     
pSelect = cellfun(@(x) findStringInAList(dataOut.parameterNames,x),ParametersToEstimate);

% selecting parameters of interest 
data.p_mapping = pSelect;

gfp = findStringInAList(dataOut.speciesNames,'[protein deGFP*]');

% cost function target species 
data.targetSpecies = gfp;
% conversion btw nM vs uM
data.speciesScaleFactor = 0.001;
% ploting simulation results with the origianl data in each step
data.debugMode =0;
% output check 
data.outputCheck = gfp;

% parameter perturbation
% select parameters to perturb:
parametersToPerturb = data.p_mapping;
data.parameterPerturbNum = 10;
% perturbation range in percentage
data.parameterPerturbRange = 0.5;
% give the indexes of parameters are subject to parameter estimation

data.parameterPerturbList = find(ismember(data.p_mapping,parametersToPerturb) == 1);
data.notParameterPerturbList = find(ismember(data.p_mapping,parametersToPerturb) == 0);

% parameter dependency rules
ntp_consumption_rate = findStringInAList(dataOut.parameterNames,'TXTL_NTP_consumption');
aa_consumption_rate = findStringInAList(dataOut.parameterNames,'TXTL_TL_AA_consumption');
K_tx = findStringInAList(dataOut.parameterNames,'TXTL_transcription_rate1');
K_tl = findStringInAList(dataOut.parameterNames,'TXTL_TL_rate');
data.p_dependecies{1} = sprintf('p(%d)=9*p(%d)',ntp_consumption_rate,K_tx);
data.p_dependecies{2} = sprintf('p(%d)=2*p(%d)',aa_consumption_rate,K_tl);

options = psoptimset('TolMesh',1e-8,'MaxIter',2000);

% parameters subject to parameter estimation
pVec = data.p0(data.p_mapping);

% default generate lower and upper bounds
default_lb = 0;
default_ub = 100;

LB = repmat(default_lb,size(pVec),1);
UB = repmat(default_ub,size(pVec),1);

% handle perturbed initial parameter cases
% Latin Hypercube sampling from a uniform distribution
p0Array = zeros(size(data.parameterPerturbList,2)+size(data.notParameterPerturbList,2),data.parameterPerturbNum);
p0Array(data.parameterPerturbList,:) = lhsu(pVec(data.parameterPerturbList)-data.parameterPerturbRange.*pVec(data.parameterPerturbList),(1+data.parameterPerturbRange).*pVec(data.parameterPerturbList),data.parameterPerturbNum)';
p0Array(data.notParameterPerturbList,:) = repmat(pVec(data.notParameterPerturbList),1,data.parameterPerturbNum);


for k = 1:data.parameterPerturbNum
    
[x(:,k),fval(k),exitflag(k),output(k)] = ...
    patternsearch(@(y)data.costFcn(y,data),p0Array(:,k),[],[],[],[],LB,UB,[],options);
end






