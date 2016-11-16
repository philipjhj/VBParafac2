

myDir = '/media/data/DataAndResults/Thesis/output/results/results_ARD_tests/';


files=dir(myDir);
files=cat(1,{files.name})';

% filename = 'ARD_tests_dim_20_20_10_4_pMethod_parafac2svd_mEsti_2_ARD_avg_datasetRNG_1_initRNG_1*.mat';


% load(filename);

pMethods={'parafac2svd','vonmises'};

ARDMethods={'max','avg','off'};

filepath = strcat(myDir,filename);

%  mEsti x Fit score x ARDmethod x pMethod x dataset x init
fitarray = zeros(7,2,3,2,10,10);
countarray = zeros(7,1,3,2,10,10);

for mEsti = 1:7
    for ARD = 1
        for P = 1:2
            for data = 1:10
                for init = 1:10
                    
%                     filename = 'ARD_tests_dim_20_20_10_4_pMethod_parafac2svd_mEsti_2_ARD_avg_datasetRNG_1_initRNG_1_29_10_2016_02_52_36.mat';
                    filename = strcat('ARD_tests_dim_20_20_10_4_pMethod_',pMethods{P},'_mEsti_',num2str(mEsti),'_ARD_',ARDMethods{ARD},'_datasetRNG_',num2str(data),'_initRNG_',...
                        num2str(init),'*.mat');
                        
                        
                        
                    filematches = regexp(files, regexptranslate('wildcard',filename),'match');
%                     disp(filematches)
                    filematches = filematches(~cellfun('isempty',filematches));
                    
                    
                    
                    if ~isempty(filematches)
%                         disp(filematches{1})
                        load(strcat(myDir,filematches{1}{1}));
                        
                        [fit, fit_true]=myModel.Parafac2Fit;
%                         disp(fit)
%                         disp(fit_true)
%                         
                        fitarray(mEsti,1,ARD,P,data,init) = fit;
                        fitarray(mEsti,2,ARD,P,data,init) = fit_true;
                        countarray(mEsti,2,ARD,P,data,init) = countarray(mEsti,1,ARD,P,data,init)+1;
                        
                    end

                end
            end
        end
    end
end



%%


% Needs to check if results were generated

avgmax = mean(max(fitarray,[],6),5);

%%
k=1;
figure(k)
% for i = 1:3
    for j = 1:2
%         subplot(3,2,j+2*(i-1))
        dataMethod = avgmax(2:end,:,2:3,j);
        dataMethod = reshape(dataMethod,size(dataMethod,1),size(dataMethod,2)*size(dataMethod,3));
%         bar(2:7,dataMethod)
        plot(2:7,dataMethod)
        axis tight
        plt = Plot();
        plt.BoxDim = [8, 8];
        plt.LineStyle = {'-.','-.'};
        plt.Colors = {'r','r','b','b','g','g'};
        plt.YLim = [70 100];
        plt.Legend=ARDMethods;
%         keyboard
        
    end
% end

%%

 % create a Plot object and grab the current figure
% plt.XLabel = 'x'; % xlabel
% plt.YLabel = 'y'; %ylabel

% plt.LineWidth = [0.001, 6];
%     plt.Colors = {'k','b'};
% plt.Markers = {'.','.'};
% plt.MarkerSpacing = [10,10];
% plt.Title = 'PCA on simple data'; % plot title
% plt.XGrid = 'on';
% plt.YGrid = 'on';

% plt.ShowBox = 'off';
% plt.XTick = [];
% plt.YTick = [];
% plt.FontSize = 35;
% plt.FileName = 'pca_example';

%%








