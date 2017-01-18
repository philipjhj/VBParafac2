


%% ARD TESTS (SYNTHETIC) NO NOISE


myAnalysis = analysisVBParafac2([]);

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

% testTitle = 'results_RealData_tests';
% testTitle = 'results_RealData_tests2';
testTitle='results_ARD_final';


myAnalysis = myAnalysis.computeTableELBOAll(resultsPath,testTitle);

load gong.mat;
sound(y)

%

% hvor mange parameter kan variere
myAnalysis.n_uniq_parameters = 1; %ARD_final=1,
myAnalysis.sortOrder = [1];

myAnalysis.find_max_ELBO_lin_idx;
myAnalysis.find_max_ELBO_sub_idx;

myAnalysis.find_max_ELBO_filenames;
%
myAnalysis.dispBestRuns
%
myAnalysis=myAnalysis.computeTableResultsFull;

load gong.mat;
sound(y)
%%

myAnalysis.computeTableFoundComponents;

%%


myAnalysis.selected_rows = [1:2];
myAnalysis.fontsize = 12;   
myAnalysis.plotTableResults


%% Normal Parafac2

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

testTitle = 'results_ARD_final_normal/';

normalParafac2congruence = zeros(1,myAnalysis.data_count);


nDataConfigs = myAnalysis.data_count;
nDataRepeats = 5;
AllResults = cell(nDataConfigs,nDataRepeats); %ARD_tests2_normalParafac2

table_high_cong = zeros(1,nDataConfigs);
table_low_cong = zeros(1,nDataConfigs);
table_fit = zeros(1,nDataConfigs);
table_nActive = zeros(1,nDataConfigs);

for dataconfig = 1:nDataConfigs
    for datasetRNG = 1:nDataRepeats
        myFile = strcat(resultsPath,testTitle,myAnalysis.data_names{dataconfig},...
            '__Normal_Parafac2_mEsti_10_datasetRNG_',num2str(datasetRNG),'.mat');
%         disp(myFile)
        load(myFile)
        
        AllResults{dataconfig,datasetRNG}.congruenceCompare = zeros(myModel.M,myModel.M); % mEsti = mTrue
        
        myModel.compute_reconstruction(data);
        
        for mEsti = 1:myModel.M
            for mTrue = 1:myModel.M
        
                Xrecon_m = myModel.Xrecon_m(:,:,:,mEsti);
                Xtrue_m = myModel.Xtrue_m(:,:,:,mTrue);
        
                congruence_m = congruenceScore(Xrecon_m(:),Xtrue_m(:));
                AllResults{dataconfig,datasetRNG}.congruenceCompare(mEsti,mTrue) = ...
                    congruence_m;
            end
        end
        
        AllResults{dataconfig,datasetRNG}.fit = myModel.Parafac2Fit(data.Xtrue);
    end
    
        % compute tables
        
    for datasetRNG = 1:nDataRepeats
        table_high_cong(dataconfig) = table_high_cong(dataconfig)...
            +sum(sum(AllResults{dataconfig,datasetRNG}.congruenceCompare>0.95));
        
        table_low_cong(dataconfig) = table_low_cong(dataconfig)...
            +sum(sum(AllResults{dataconfig,datasetRNG}.congruenceCompare>0.85));
        
        table_fit(dataconfig) = table_fit(dataconfig)+AllResults{dataconfig,datasetRNG}.fit;
    end
    
end

% same order as plot function
parafac2tables{1} = table_low_cong/nDataRepeats;
parafac2tables{2} = table_high_cong/nDataRepeats;
parafac2tables{5} = table_fit/nDataRepeats;
%%

myAnalysis.selected_rows = [1:2];
myAnalysis.fontsize = 12;   
myAnalysis.plotTableResults(parafac2tables)


%%









%% REAL DATA TESTS


% Only Apple data set to begin with!

myAnalysis = analysisVBParafac2([]);

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

testTitle = 'results_RealData_tests3';


myAnalysis = myAnalysis.computeTableELBOAll(resultsPath,testTitle);

load gong.mat;
sound(y)

%%

% hvor mange parameter kan variere
myAnalysis.n_uniq_parameters = 1; %ARD_final=1,
myAnalysis.sortOrder = [1];


myAnalysis.find_max_ELBO_lin_idx;
myAnalysis.find_max_ELBO_sub_idx;

myAnalysis.find_max_ELBO_filenames;
%%
myAnalysis.dispBestRuns
%%
myAnalysis=myAnalysis.computeTableResultsFull;

load gong.mat;
sound(y)
%%
myAnalysis.sortOrder = [1];
myAnalysis.computeTableFoundComponents;

%%


myAnalysis.selected_rows = [1:2];
myAnalysis.fontsize = 12;   
myAnalysis.plotTableResults


%% Normal Parafac2

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

testTitle = 'results_ARD_tests2_normalParafac2/';

normalParafac2congruence = zeros(1,myAnalysis.data_count);


nDataConfigs = myAnalysis.data_count;
nDataRepeats = 5;
AllResults = cell(nDataConfigs,nDataRepeats); %ARD_tests2_normalParafac2

table_high_cong = zeros(1,nDataConfigs);
table_low_cong = zeros(1,nDataConfigs);
table_fit = zeros(1,nDataConfigs);
table_nActive = zeros(1,nDataConfigs);

for dataconfig = 1:nDataConfigs
    for datasetRNG = 1:nDataRepeats
        myFile = strcat(resultsPath,testTitle,myAnalysis.data_names{dataconfig},...
            '__Normal_Parafac2_mEsti_10_datasetRNG_',num2str(datasetRNG),'.mat');
%         disp(myFile)
        load(myFile)
        
        AllResults{dataconfig,datasetRNG}.congruenceCompare = zeros(myModel.M,myModel.M); % mEsti = mTrue
        
        myModel.compute_reconstruction(data);
        
        for mEsti = 1:myModel.M
            for mTrue = 1:myModel.M
        
                Xrecon_m = myModel.Xrecon_m(:,:,:,mEsti);
                Xtrue_m = myModel.Xtrue_m(:,:,:,mTrue);
        
                congruence_m = congruenceScore(Xrecon_m(:),Xtrue_m(:));
                AllResults{dataconfig,datasetRNG}.congruenceCompare(mEsti,mTrue) = ...
                    congruence_m;
            end
        end
        
        AllResults{dataconfig,datasetRNG}.fit = myModel.Parafac2Fit(data.Xtrue);
    end
    
        % compute tables
        
    for datasetRNG = 1:nDataRepeats
        table_high_cong(dataconfig) = table_high_cong(dataconfig)...
            +sum(sum(AllResults{dataconfig,datasetRNG}.congruenceCompare>0.95));
        
        table_low_cong(dataconfig) = table_low_cong(dataconfig)...
            +sum(sum(AllResults{dataconfig,datasetRNG}.congruenceCompare>0.85));
        
        table_fit(dataconfig) = table_fit(dataconfig)+AllResults{dataconfig,datasetRNG}.fit;
    end
    
end

% same order as plot function
parafac2tables{1} = table_low_cong/nDataRepeats;
parafac2tables{2} = table_high_cong/nDataRepeats;
parafac2tables{4} = table_fit/nDataRepeats;
%%

myAnalysis.selected_rows = [1:2];
myAnalysis.fontsize = 12;   
myAnalysis.plotTableResults(parafac2tables)







