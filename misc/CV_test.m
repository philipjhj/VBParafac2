% addpath(genpath('VBParafac2'))

test_data='real';
if strcmp(test_data,'synthetic')
    I=150;
    J=I;
    K=10;
    M=4;
    Mesti = M;
    
    options.dimensions = [I J K M];
    options.initMethod = 'kiers';
    % options.initMethod = 'generative';
    options.congruence = 0.4;
    % 1e4 1e-3 i ARD tests
    options.precision = [1e2 1e-6];
    options.SNR = 5;
    options.noiseType = 'homo';
    % [1e4 1e-8] creates problems for qC
    
    % sumSNR = 0;
    % for k = 1:100
    rng('shuffle')
    % rng(1)
    data = varBayesModelParafac2.generateDataFromModel(options);
    
elseif strcmp(test_data,'real')
    load('/media/data/DataAndResults/VBParafac2paper/amino.mat')
    data = dataClass;
    data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);
    Mesti=3;
    
end


myModel=varBayesModelParafac2(data,Mesti);

myModel.opts.verbose = 0;
myModel.opts.debugFlag = 0;
myModel.opts.estimationP= 'parafac2svd';
% myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'maxNoARD';
myModel.opts.estimationNoise = 'max';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 10;
% myModel.opts.rngInput = 31; % bad ones: 8, good ones: 15, 16, 19, 21, CORRECT: 31
myModel.opts.maxIter = 1000;
% myModel.opts.maxTime = 4;

myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};

rng('default')
Ms = 1:7
myModel.crossValidateM(Ms)