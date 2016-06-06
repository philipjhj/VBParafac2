myModel=varBayesModelParafac2;

% myModel.qDist.method = 'vonmises';

myModel.qDist.SNR

% qA qC qF qP qSigma qAlpha
myModel.qDist.activeParams = {'qA','qC','qP','qF','qSigma','qAlpha'};

% Kandidater til fejl:
% F mean opdateres ikke korrekt (er udkommenteret), check udledning
%

% Plan:
% Overvej beregning af gamma mean fra inverse gamma parameters
% Optimer program:
% ## Omskriv manopt funktioner til konstante værdier
% ## Beregn ens værdier en gang (f.eks. variance matricer)
% ## Restrukturer kode (find ens led og skriv dem i en funktion)
% ## Profile og optimer beregninger
% ## GPU?


myModel.computeVarDistribution;
myModel.qDist.SNR


%%
for k=1:myModel.data.K;
    myModel.qDist.qP.mean(:,:,k)'*myModel.qDist.qP.mean(:,:,k) %+ sum(myModel.qDist.qP.variance(:,:,:,k),3)
end
