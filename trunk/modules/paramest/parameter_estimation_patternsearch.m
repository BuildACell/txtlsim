%% load datasets

pr_pr1_pr2_27 = data_importer('data/pr_pr1_pr2_02_27.csv','victor');
pr_pr1_pr2_28 = data_importer('data/pr_pr1_pr2_02_28.csv','victor');
pr_pr1_pr2_04 = data_importer('data/pr_pr1_pr2_03_04.csv','victor');


% pr
[data_1nm_mean data_1nm_std]  = getStdMean([pr_pr1_pr2_27.noBg(:,1) pr_pr1_pr2_28.noBg(:,1) pr_pr1_pr2_04.noBg(:,1)]);
[data_2nm_mean data_2nm_std]  = getStdMean([ pr_pr1_pr2_28.noBg(:,3) pr_pr1_pr2_04.noBg(:,3)]);
[data_3nm_mean data_3nm_std]  = getStdMean([pr_pr1_pr2_27.noBg(:,4) pr_pr1_pr2_04.noBg(:,4)]);
[data_4nm_mean data_4nm_std]  = getStdMean([ pr_pr1_pr2_28.noBg(:,5) pr_pr1_pr2_04.noBg(:,5)]);


data.xdata = pr_pr1_pr2_27.t_vec ;
data.ydata = [data_1nm_mean./(225699.04) data_2nm_mean./(225699.04) data_3nm_mean./(225699.04) data_4nm_mean./(225699.04)];

%% initial values and parameters

dataOut = parseGetEqOutput(Mobj);

data.x0 = dataOut.initialValues;
data.p0 = dataOut.parameters;

data.modelFcn = str2func(dataOut.modelFcn);
data.costFcn = @extract_calibration_cost_fcn;

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

K_tx = findStringInAList(dataOut.parameterNames,'TXTL_transcription_rate1');
ntpdeg = findStringInAList(dataOut.parameterNames,'NTPdeg_F');
rbs_f = findStringInAList(dataOut.parameterNames,'TXTL_UTR_RBS_F');

% selecting parameters of interest 
data.p_mapping = [K_tx,ntpdeg,rbs_f];

gfp = findStringInAList(dataOut.speciesNames,'[protein deGFP*]');

% cost function target species 
data.targetSpecies = gfp;
% conversion btw nM vs uM
data.speciesScaleFactor = 0.001;
% ploting simulation results with the origianl data in each step
data.debugMode =0;
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






