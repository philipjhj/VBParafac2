warning off MATLAB:nearlySingularMatrix
addpath(genpath('misc'))

DATA_SET_PATH='example_data.mat';
M_ESTIMATE=5;

% For example data
if ~isfile(DATA_SET_PATH);
    simulate_example_data
end

load(DATA_SET_PATH) % Should contain X, samples as third mode

% DATA PREPROCESSING
% Not needed for the example data
% If you load an data matrix "X" (I x J x K with K = number of samples)
% use the code below to prepare the data
% data = dataClass;
% data.Xunfolded = X;

% INITIALIZE MODEL
myModel=varBayesModelParafac2(data, M_ESTIMATE);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 0;
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 500;
myModel.opts.maxTime = 60*60*9;
myModel.opts.maxIter = 1e4;
myModel.opts.initMethod = 'mle';
myModel.opts.noiseLearningDelay = 50;
myModel.opts.scaleLearningDelay = 0;
myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'max';
myModel.opts.estimationNoise = 'avg';
myModel.opts.rngInput = 1;
myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};

rng('default'); rng('shuffle');

myModel.partitionData(myModel.fullData.X)
trainTic=tic;
myModel.fitTrainingData;
toc(trainTic)

% SAVE TRAINED MODEL
% save(OUTPUT_PATH);

print('success');