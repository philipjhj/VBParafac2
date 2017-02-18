

pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;

Ms = 2:10;
M = numel(Ms);
T = 5;
evaltimeAvg = zeros(M,T,numel(SNR),numel(noiseType),numel(pMethod));

for pMethodIdx=1:2
    for noiseIdx = 1:2
        for SNRIdx = 1:numel(SNR)
            if strcmp(noiseType{noiseIdx},'hetero') && SNR(SNRIdx)==0 && strcmp(pMethod{pMethodIdx},'parafac2svd')
                continue
            end
            load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_sim_data__dim_150_150_30_4_04__pMethod_',...
                pMethod{pMethodIdx},'_SNR_',num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat'))
            
            for m = 1:M
                for t = 1:T
                    evaltimeAvg(m,t,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{t}.dataTrain.n_components_hard(end);
%                     evaltimeAvg(m,t,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{t}.dataTrain.evaltime(end);
                end
            end  
        end
    end
end
%%
evaltimeAvgMean = squeeze(mean(evaltimeAvg,2));

for p = 1:2
    for n = 1:2
        for s = 1:4
            figure(p)
            subplot(4,2,sub2ind([2 4],n,s))
            bar(evaltimeAvgMean(:,s,n,p))
            ylabel('s');
            title(strcat(pMethod{p},'~',noiseType{n},'~',num2str(SNR(s),'SNR=%d')))
        end
        
    end
    
end










%%
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