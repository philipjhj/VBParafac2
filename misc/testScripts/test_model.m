%load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
% set_wd(2)
% set(0,'DefaultFigureWindowStyle','docked')
%load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Int2.mat')
%%

warning off MATLAB:nearlySingularMatrix
%%
% myModel=varBayesModelParafac2;
rng('default')
I=150;
J=I;
K=30;
M=4;
Mesti = 10;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
% options.initMethod = 'generative';
options.congruence = 0.4;
% 1e4 1e-3 i ARD tests
options.precision = [1e2 1e-6];
options.SNR = -12;
options.noiseType = 'homo';
% [1e4 1e-8] creates problems for qC

% sumSNR = 0;
% for k = 1:100
% rng('shuffle')
rng(1)
data = varBayesModelParafac2.generateDataFromModel(options);

%
% sumSNR = sumSNR+10*log10(norm(data.Xtrue(:),'fro')^2/norm(data.Etrue(:),'fro')^2);
% end
% sumSNR/100
% I2 = I2/(norm(I2(:),'fro')/numel(I2));
% data.X = permute(I2,[2 1 3]);
%
%
%
%normalModel = normalParafac2(data.X);
% normalModel = normalParafac2(permute(I1,[2 1 3]));

%normalModel.fitParafac2(4)
%
%
%normalModel.Parafac2Fit(data.Xtrue)
%
rng('default')
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
myModel.opts.estimationNoise = 'avg';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 1;
myModel.opts.rngInput = 7;
% myModel.opts.maxIter = 50;

myModel.opts.activeParams = {'qA','qF','qP','qC','qSigma','qAlpha'};
%
myModel.partitionData(myModel.fullData.X)
% tic
myModel.fitTrainingData;
% toc

% myModel.opts.activeParams = {'qA','qC','qP','qSigma','qF'};
%%
Ms = 2:2;
myModel.crossValidateM(Ms)
%%

%%
best_ELBO_true = squeeze(mean(max(myModel.CV_ELBOS(:,:,3),[],2)));
squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true
(squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true)/best_ELBO_true

%%
best_ELBO_all = squeeze(mean(max(myModel.CV_ELBOS,[],2)));
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
return
for k=1:myModel.data.K
    
clf

myModel.plotSolutionSynthK(k,0)
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
