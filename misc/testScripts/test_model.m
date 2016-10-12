load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%

% warning on MATLAB:nearlySingularMatrix

% myModel=varBayesModelParafac2;

I=50;
J=50;
K=10;
M=2;
Mesti = 8;

dimensions = [I J K M];
% rng(3)
data = varBayesModelParafac2.generateDataFromModel([I J K M],[1e12 1e-6]);


myModel=varBayesModelParafac2(data,Mesti);

% myModel=varBayesModelParafac2(Y,100);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 1;
myModel.opts.estimationP= 'vonmises';
% myModel.opts.estimationP = 'parafac2svd';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.showIter = 1;
% myModel.maxTime = realmax;

%myModel.qDist.SNR
clc
myModel.qDist.opts.activeParams = {'qA','qF','qC','qP','qAlpha','qSigma'};
% myModel.qDist.activeParams_opt = {'qC','qAlpha'};


% clc

% myModel.data.iter = myModel.data.iter-1;

tic
% rng(3)
myModel.computeVarDistribution(400);
toc
%myModel.qDist.SNR

%%
figure
for k=1:myModel.data.K
    
clf

myModel.plotSolutionSynthK(k,0)
pause
end
%%
clf
myModel.plotELBO([100 myModel.data.iter-1])


%%
for k = 1:1
    
    clf
    myModel.plotSolutionRealMatrixK(k)
%    pause
   
end

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
