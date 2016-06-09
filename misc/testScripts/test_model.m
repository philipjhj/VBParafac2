myModel=varBayesModelParafac2;

myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.SNR

% qA qC qF qP qSigma qAlpha
myModel.qDist.activeParams = {'qC','qF','qA','qAlpha','qSigma','qP'};

% Kandidater til fejl:
% Problemer med VB:
% Sigma og Alpha skal korrekt parameteriseres som precisionen (done)
% qF har problemer når Mtrue > M
% Entropy af multinormals bliver negative?


% Scaling
% Vonmises: 
% - only have problems with scaling in I and J
% SVD:
% - Problemer med updates for P ved scaling af M 


myModel.computeVarDistribution;
myModel.qDist.SNR


%%
for k=1:myModel.data.K;
    myModel.qDist.qP.mean(:,:,k)'*myModel.qDist.qP.mean(:,:,k) %+ sum(myModel.qDist.qP.variance(:,:,:,k),3)
end













% Plan:
% Forsøg at lav MLE på F i stedet for et e-step.
% Optimer program:
% ## Beregn ens værdier en gang (f.eks. variance matricer)
% ## Restrukturer kode (find ens led og skriv dem i en funktion)
% ## Profile og optimer beregninger
% ## GPU?
