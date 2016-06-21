load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%

warning on MATLAB:nearlySingularMatrix

% myModel=varBayesModelParafac2;

I=10;
J=I;
K=10;
M=2;
dimensions = [I J K M];
data = varBayesModelParafac2.generateDataFromModel([I J K M]);
Mesti = 2;
myModel=varBayesModelParafac2(data,Mesti);


% myModel=varBayesModelParafac2(Y(:,:,1:3),5);


myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.SNR

myModel.qDist.activeParams_opt = {'qA','qC','qP','qF','qSigma','qAlpha'};
tic
myModel.computeVarDistribution(500);
toc
myModel.qDist.SNR


%%
clf
myModel.plotSolutionSynthK(2,1)


%%
for k = 1:3
    
    clf
    myModel.plotSolutionRealK(k)
   pause
   
end



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
