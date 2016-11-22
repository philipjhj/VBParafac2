

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
myAnalysis = analysisVBParafac2([]);

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

% testTitle = 'results_RealData_tests';
testTitle= 'results_ARD_tests';
% testTitle= 'test_folder';

myAnalysis = myAnalysis.computeTableELBOAll(resultsPath,testTitle);

load gong.mat;
sound(y)

%%

myAnalysis.find_max_ELBO_lin_idx;
myAnalysis.find_max_ELBO_sub_idx;

myAnalysis.find_max_ELBO_filenames;
%%
myAnalysis.dispBestRuns
%%
myAnalysis=myAnalysis.computeTableResultsFull;
%%
myAnalysis.computeTableFoundComponents;


imagesc((myAnalysis.nFoundComponents))
colorbar
axis image
yticks(1:36)


testConfig = sortrows(...
    unique(myAnalysis.max_ELBO_sub_idx(:,1:3),'rows'),myAnalysis.sortOrder);
ynames = cell(1,size(testConfig,1));

for i = 1:size(testConfig,1)
    
    myString = '';
    for j = 1:3 
        myString = sprintf('%s %s',myString,...
            myAnalysis.testOpts.(myAnalysis.testOpts_names{j}){...
            testConfig(i,j)});
    end
    ynames{i} = myString;
end

% ynames = unique(cat(1,ynames{:}));

set(gca,'YTickLabel',ynames)
xlabel('Data Sets')
%%
close all
%%
myAnalysis.fontsize = 12;
myAnalysis.plotTableResults
close all
%%
%%
myAnalysis.plotARDtest

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






