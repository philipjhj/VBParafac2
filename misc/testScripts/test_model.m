load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%

% warning on MATLAB:nearlySingularMatrix

% myModel=varBayesModelParafac2;

I=50;
J=I;
K=5;
M=2;
Mesti = 20;

dimensions = [I J K M];
data = varBayesModelParafac2.generateDataFromModel([I J K M]);


myModel=varBayesModelParafac2(data,Mesti);


% myModel=varBayesModelParafac2(Y,100);


myModel.qDist.debugflag = 0;
myModel.verbose = 1;
myModel.qDist.method = 'vonmises';
% myModel.qDist.method = 'parafac2svd';


%myModel.qDist.SNR

myModel.qDist.activeParams_opt = {'qA','qC','qF','qP','qAlpha','qSigma'};
% myModel.qDist.activeParams_opt = {'qC','qAlpha'};

% clc
tic
myModel.computeVarDistribution(6500);
toc
%myModel.qDist.SNR
%%

for k=1:1%myModel.data.K
clf
myModel.plotSolutionSynthK(k,0)
% pause
end

%%
for k = 1:1
    
    clf
    myModel.plotSolutionRealMatrixK(k)
%    pause
   
end

%%
m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
mask = m.mask;

k = 1;
m = 1:5;

clf
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
