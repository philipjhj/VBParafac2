

myAnalysis = analysisVBParafac2([]);
%%


load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Fmax_it_models.mat')
load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Int1.mat')

normalModel = normalParafac2(I1,A{1,5},C{1,5},H{1,5},P{1,5});
%%
saveFlag = 1;

% myAnalysis.plotReconElutionProfiles(myModel,'test',saveFlag)
myAnalysis.plotReconElutionProfiles(normalModel,'NormalParafac2_M_10',saveFlag)

%%
close all


%%


testDir='/media/data/DataAndResults/Thesis/output/results/results_ARD_tests/';

myAnalysis.computeTestResults(testDir)


%%

subplot(2,2,1)
imagesc(normalModel.A)
colorbar
title('Normal Parafac2 - A')
subplot(2,2,2); 
imagesc(myModel.qDist.qA.mean)
title('VB Parafac2 - A')
colorbar

subplot(2,2,3)
imagesc(normalModel.C)
colorbar
title('Normal Parafac2 - C')
subplot(2,2,4); 
imagesc(myModel.qDist.qC.mean)
colorbar
title('VB Parafac2 - C')






