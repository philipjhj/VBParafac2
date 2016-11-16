

myAnalysis = analysisVBParafac2([]);
%%


load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Fmax_it_models.mat')
load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Int1.mat')

normalModel = normalParafac2(I1,A{1,5},C{1,5},H{1,5},P{1,5});
%%
saveFlag = 1;

% myAnalysis.plotReconElutionProfiles(myModel,fname,saveFlag)
myAnalysis.plotReconElutionProfiles(normalModel,'Normal_Parafac2_optimal_Apple_Int1',saveFlag)

%%
close all


%%


testDir='/media/data/DataAndResults/Thesis/output/results/results_ARD_tests/';

myAnalysis.computeTestResults(testDir)


