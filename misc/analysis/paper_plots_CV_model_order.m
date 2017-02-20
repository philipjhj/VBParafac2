



%%

pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;

Ms = 2:10;
M = numel(Ms);
T = 5;
allTestValues = zeros(M,numel(SNR),numel(noiseType),numel(pMethod));

for pMethodIdx=1:2
    for noiseIdx = 1:2
        for SNRIdx = 1:numel(SNR)
            filename = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/EarlyResultsToBeDeleted/par_CV_sim_data__dim_150_150_30_4_04__pMethod_',...
                pMethod{pMethodIdx},'_SNR_',num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat');
            if exist(filename, 'file') == 2
                load(filename);
            else
                continue
            end
            for m = 1:M
                allTestValues(m,SNRIdx,noiseIdx,pMethodIdx) = mean(max(myModel{m}.CV_ELBOS_test,[],2));
            end
            
            figure(pMethodIdx)
            subplot(1,2,noiseIdx)
            hold on
            plot(Ms,allTestValues(:,SNRIdx,noiseIdx,pMethodIdx),'o-','Color',colors(SNRIdx))
            grid on
            
            xlabel('Number of components')
            ylabel('ELBO')
            title(noiseType{noiseIdx})
        end
        legend(strtrim(cellstr(num2str(SNR','SNR=%d'))))
    end
end

subplot(1,2,1)
ylimits(1,1:2)=ylim;
subplot(1,2,2)
ylimits(2,1:2)=ylim;

ybot = min(ylimits(:,1));
ytop = max(ylimits(:,2));

subplot(1,2,1)
ylim([ybot ytop]);
subplot(1,2,2)
ylim([ybot ytop]);


% plot(myModel{m}{t}.dataTrain.n_components_hard)