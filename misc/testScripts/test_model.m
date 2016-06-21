load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%

warning on MATLAB:nearlySingularMatrix

myModel=varBayesModelParafac2;%(Y(:,:,:),30);
myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.SNR

myModel.qDist.activeParams_opt = {'qA','qC','qP','qF','qSigma','qAlpha'};

myModel.qDist.activeParams = {'qA','qF','qP','qAlpha','qSigma','qC'};
% myModel.qDist.activeParams = {'qF'};
tic
myModel.computeVarDistribution;
toc
myModel.qDist.SNR
%%

figure(2)
plot(nonzeros(myModel.ELBO_chain));

%%
myModel.restartqDist
myModel.computeVarDistribution

%%

for k = 1:1
    clf
    
    myModel.plotSolution(k,0)
%    pause
end



%%

[A,F,C,P,fit]=parafac2(myModel.data.X,myModel.data.M,[0 0]);
%%
clf
plotParafac2SolutionK(1,myModel.data.X,A,C,F,cat(3,P{:}),myModel.data.Atrue,...
    myModel.data.Ctrue,myModel.data.Ftrue,myModel.data.Ptrue)
