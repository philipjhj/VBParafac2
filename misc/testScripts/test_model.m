load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%
% warning off MATLAB:nearlySingularMatrix

myModel=varBayesModelParafac2;%(Y(:,:,:),10);

myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.SNR

% qA qC qF qP qSigma qAlpha
myModel.qDist.activeParams_opt = {'qF','qA','qC','qP','qAlpha','qSigma'};

% Kandidater til fejl:
% Problemer med VB:
% qP variance bliver tæt på singular når M sættes højere end Mtrue. Kan
% også ske andre gange (på syntetisk data), men skaber ikke fejl.

% Scaling
% Parfor loops
% GPU

% Tests Kode
% Kør tests for at verificere at koden kører korrektbb

% Tests Data
% Synthetic:
% - Vis læring af komponenter
% - Vis ARD
% - Variere SNR
% - Forskellig støj pr. subjekt
% Rasmus bro's data
% Hjernedata
% Sammenlig med N-way toolbox'ens løsning

myModel.computeVarDistribution;
myModel.qDist.SNR
%%
figure(2)
plot(nonzeros(myModel.ELBO_chain));

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
