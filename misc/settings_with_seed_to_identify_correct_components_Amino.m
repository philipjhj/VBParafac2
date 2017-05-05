


load('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFit_real_data__AminoAcid_amino__mEsti_3.mat')
normalModel=myModel{2}

load('/media/data/DataAndResults/VBParafac2paper/amino.mat')
data = dataClass;
data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);
Mesti=3;





myModel=varBayesModelParafac2(data,Mesti);
%
% size(myModel.data.X)
%
% myModel=varBayesModelParafac2(Y,100);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 0;
myModel.opts.estimationP= 'parafac2svd';
% myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'maxNoARD';
myModel.opts.estimationNoise = 'max';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 10;

% ELBO with rng=16: -4.16e+08, ELBO with rng=31: -4.17e+08

% rng=31 identifies the correct components, but it has a lower ELBO than
% many of the bad solutions. Also, its iterations contains errors with a
% relative magnitude around 1e-12

% bad ones: 8, good ones: 15, 16, 19, 21, 56 CORRECT: 31, 986461368, 986471051
% myModel.opts.rngInput = 986471051; 
myModel.opts.maxIter = 1000;
% myModel.opts.maxTime = 4;

myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};
%
rng('default')
rng(1)
myModel.partitionData(myModel.fullData.X)
% tic
myModel.fitTrainingData;
% toc



bestVBmodel=myModel;
% bestVBmodel=myModel{m}{idx};

CC=corrcoef([bestVBmodel.qDistTrain.qC.mean y]);
CC=abs(CC(1:3,4:end));

[~,IDXs]=max(CC);
CCsorted=CC(IDXs,:);

CCC=corrcoef([normalModel.C y]);
CCC=abs(CCC(1:3,4:end));

[~,IDXs]=max(CCC);
CCCsorted=CCC(IDXs,:);

% figure(1);
subplot(2,1,1);
imagesc(CC); 
colorbar
% figure(2);
subplot(2,1,2);
imagesc(CCC);
colorbar






