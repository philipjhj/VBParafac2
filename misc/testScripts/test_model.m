%load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
% set_wd(2)
% set(0,'DefaultFigureWindowStyle','docked')
load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Int2.mat')
%%

warning on MATLAB:nearlySingularMatrix

%%

load('/media/data/DataAndResults/VBParafac2paper/data/aminoAcid/amino.mat')
data = dataClass;
data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);
Mesti=3;
%%
% myModel=varBayesModelParafac2;

I=20;
J=I;
K=4;
M=2;
Mesti = M;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
% options.initMethod = 'generative';
options.congruence = 0.4;
% 1e4 1e-3 i ARD tests
options.precision = [1e2 1e-6];
options.SNR = -4;
options.noiseType = 'homo';
% [1e4 1e-8] creates problems for qC

% sumSNR = 0;
% for k = 1:100
rng('shuffle')
% rng(1)
data = varBayesModelParafac2.generateDataFromModel(options);

%
% sumSNR = sumSNR+10*log10(norm(data.Xtrue(:),'fro')^2/norm(data.Etrue(:),'fro')^2);
% end
% sumSNR/100
% I2 = I2/(norm(I2(:),'fro')/numel(I2));
% data.X = permute(I2,[2 1 3]);
%
% %
% %
% data=dataClass;
% data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);

% normalModel = normalParafac2(data.X);
% normalModel = normalParafac2(permute(I2,[2 1 3]));

%%
% for m = 2:5
    rng('default')
%normalModel.fitParafac2(3)
%
% normalModel.CCDParafac2
%  end
%
%
% parafac2(data.Xtrue,Mesti,[0 0],[0 0 0 0 0]);
% normalModel.Parafac2Fit(data.Xtrue)
%
% rng('default')
myModel=varBayesModelParafac2(data,Mesti);
%
% size(myModel.data.X)
%
% myModel=varBayesModelParafac2(Y,100);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 0;
myModel.opts.estimationP= 'parafac2svd';
% myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'max';
myModel.opts.estimationNoise = 'maxShared';
myModel.opts.initMethod = 'mle';
myModel.opts.noiseLearningDelay=20;
myModel.opts.scaleLearningDelay=10;

myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 1;
% myModel.opts.rngInput = 15; % bad ones: 8, good ones: 15, 16, 19, 21, CORRECT: 31
myModel.opts.maxIter = 5000;
% myModel.opts.maxTime = 4;

myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};
%%
rng('default')
% rng(1)
myModel.partitionData(myModel.fullData.X)
% tic
myModel.fitTrainingData;
% toc 

% myModel.opts.activeParams = {'qA','qC','qP','qSigma','qF'};
%%
Ms = 1:5
rng('default')
myModel.crossValidateM(Ms)
%%

cellfun( @(S) S.Data.stopReason, myModel.cvRunsTrain, 'uni', false )

ELBOS = cell2mat(cellfun( @(S) S.Data.ELBO, myModel.cvRunsTrain, 'uni', false ));
ELBOS-max(ELBOS,[],2)

%%
best_ELBO_true = squeeze(mean(max(myModel.CV_ELBOS(:,:,3),[],2)));
squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true
(squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true)/best_ELBO_true

%%
best_ELBO_all = squeeze(mean(max(myModel.CV_ELBOS_test,[],2)));
plot(Ms,best_ELBO_all,'o-')
xlabel('# of components')
ylabel('ELBO')
grid on
% figname=input('Write name of figure: ','s');
% saveas(gcf,strcat('output/paper/temps/',figname,'.jpg'))

%%
tic
myModel.fitTrainingData;
toc
%myModel.qDist.SNR
% myModel.Parafac2Fit
%%
% return
for k=1:myModel{3}{1}.fullData.K
    
clf

myModel{2}{1}.plotSolutionSynthK(k,0)
% myModel{3}{1}.Parafac2Fit(myModel{3}{1}.qDistTrain)
pause
end
%%
clf
myModel.plotELBO([100 myModel.data.iter-1])


%%
for k = 1:4
    
    clf
    myModel.plotSolutionRealMatrixK(k)
   keyboard
   
end
%%
figure
subplot(2,1,1)
imagesc(myModel.qDist.qC.mean)
colorbar
subplot(2,1,2)
imagesc(1./myModel.qDist.qAlpha.mean) % Lav som procent varians i stedet
colorbar
%%
m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
mask = m.mask;
set(0,'DefaultFigureWindowStyle','docked')
%%
close all
k = 1;
m = 1:25;

% close all
myModel.plotSolutionReal3D(k,m,mask)


%%

figure(2)
plot(nonzeros(myModel.ELBO_chain));

%%
myModel.restartqDist
myModel.computeVarDistribution


%%
k=2;
U = myModel.qDist.qA.mean*diag(myModel.qDist.qC.mean(k,:));
V = myModel.qDist.qP.mean(:,:,k)*myModel.qDist.qF.mean;

m=5;
plotComponent(U(:,1:m),V(:,1:m)', mask, [ 53    63    46])



%%

[A,F,C,P,fit]=parafac2(myModel.data.X,myModel.data.M,[0 0]);
%%
clf
plotParafac2SolutionK(1,myModel.data.X,A,C,F,cat(3,P{:}),myModel.data.Atrue,...
    myModel.data.Ctrue,myModel.data.Ftrue,myModel.data.Ptrue)
