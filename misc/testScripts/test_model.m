load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%

% warning off MATLAB:nearlySingularMatrix

% myModel=varBayesModelParafac2;

I=20;
J=I;
K=6;
M=3;
Mesti = 3;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
options.congruence = 0.4;
options.precision = [1e12 1e-3];

rng(3)
data = varBayesModelParafac2.generateDataFromModel(options);


myModel=varBayesModelParafac2(data,Mesti);

% myModel=varBayesModelParafac2(Y,100);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 2;
myModel.opts.estimationP= 'vonmises';
% myModel.opts.estimationP = 'parafac2svd';
myModel.opts.estimationARD = 'avg';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 1;
myModel.opts.rngInput = 1463436224;%'shuffle';

% data set; rng(3)
% seed; 1461191309
% iter (avg); 441
% iter (max); 44

% myModel.opts.maxTime = 1;




%myModel.qDist.SNR
% clc

myModel.qDist.opts.activeParams = {'qA','qF','qC','qP','qAlpha','qSigma'};
% myModel.qDist.activeParams_opt = {'qC','qAlpha'};


% clc

% myModel.data.iter = myModel.data.iter-1;
% myModel.restartqDist;
% myModel.opts.maxTime = 5;
tic
myModel.computeVarDistribution;
toc
%myModel.qDist.SNR
myModel.Parafac2Fit
%%

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
