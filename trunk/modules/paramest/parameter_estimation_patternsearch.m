%% load datasets

% e15_pr_gfp_mg_0821 = wholeExpfileReader('ACS_paper_data/e15_pr_gfp_pr_gfp_mg_grad_0821_data.csv',';',2);
% e15_pr_gfp_mg_0821 = processMGData(e15_pr_gfp_mg_0821);
% e15_pr_gfp_mg_0821.valveNames = {'pos cont','neg cont','pr-gfp 1nM','pr-gfp 2nM','pr-gfp 3nM','pr-gfp 5nM','pr-gfp 10nM','pr-gfp 15nM'...
%                                  'pr-gfp-s15-mg 1nm','pr-gfp-s15-mg 2nM','pr-gfp-s15-mg 3nM','pr-gfp-s15-mg 5nM','pr-gfp-s15-mg 10nM','pr-gfp-s15-mg 15nM'};
%  

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
% select a species
DNAind = findStringInAList(dataOut.speciesNames,'[DNA p70--rbs--deGFP]');
% specify the initial values for the selected species  
DNAinitialValues = [1 2 3 4]; % 1-4nM
% transform x0 vector to the right format (each column -> one different initial values vector)
data.x0 = repmat(data.x0,1,size(DNAinitialValues,2)); 
data.x0(DNAind,:) = DNAinitialValues;

% which species are subject to initial value estimation
%data.x0_mapping = DNAind;
data.x0_mapping = [];
% adjusting parameter vector size based on the initial values
if ~isempty(data.x0_mapping)
    data.p0 = [data.p0; data.x0(DNAind,:)'];
end

% select parameters for estimation

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
data.debugMode =1;
% output check 
data.outputCheck = gfp;

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


[x,fval,exitflag,output] = patternsearch(@(y)data.costFcn(y,data),pVec,[],[],[],[],LB,UB,[],options);






