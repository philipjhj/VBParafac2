

pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;

Ms = 2:10;
M = numel(Ms);
T = 5;
% evaltimeAvg = zeros(M,T,numel(SNR),numel(noiseType),numel(pMethod));
alphaMeans = zeros(10,M,numel(SNR),numel(noiseType),numel(pMethod));

for pMethodIdx=1:2
    for noiseIdx = 1:2
        for SNRIdx = 1:numel(SNR)
            if strcmp(noiseType{noiseIdx},'hetero') && SNR(SNRIdx)==0 && strcmp(pMethod{pMethodIdx},'parafac2svd')
                continue
            end
            load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/EarlyResultsToBeDeleted/ARD_sim_data__dim_150_150_30_4_04__pMethod_',...
                pMethod{pMethodIdx},'_SNR_',num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat'))
            
            for m = 1:M
                maxELBO = -realmax;
                bestT = -1;
                for t = 1:T
%                     keyboard
                    if myModel{m}{t}.dataTrain.ELBO > maxELBO
                        bestT = t;
                    end
%                     evaltimeAvg(m,t,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{t}.dataTrain.n_components_hard(end);
%                     evaltimeAvg(m,t,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{t}.dataTrain.evaltime(end);
                end
                alphaMeans(1:(m+1),m,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{bestT}.qDistTrain.qAlpha.mean';
            end  
        end
    end
end
%% ELBO curves
evaltimeAvgMean = squeeze(mean(evaltimeAvg,2));

for p = 1:2
    for n = 1:2
        for s = 1:4
            figure(p)
            subplot(4,2,sub2ind([2 4],n,s))
            bar(evaltimeAvgMean(:,s,n,p))
            ylabel('s');
            title(strcat(pMethod{p},'~',noiseType{n},'~',num2str(SNR(s),'SNR=%d')))
        endmyModel{m}{t}.qDistTrain.qAlpha.mean
        
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

%% Hinton diagram

varianceValues = 1./alphaMeans;
varianceValues(varianceValues==Inf) = 0;
varianceValues = sort(varianceValues,1,'ascend');

for p = 1:2
    for n = 1:2
        for s = 1:4
            hinton(varianceValues(:,:,s,n,p),sub2ind([2, 4],n,s),[4 2]);
        end
    end
end
%% Hinton 2
varianceValues = 1./alphaMeans;
varianceValues(varianceValues==Inf) = 0;
varianceValues = sort(varianceValues,1,'ascend');
varianceValues = permute(varianceValues,[1 3 2 4 5]);
varDim = size(varianceValues);

varianceValues = reshape(varianceValues,prod(varDim([1 2])),prod(varDim([3 4])),2);
varianceValues = [varianceValues(:,1:9,1) zeros(size(varianceValues,1),1) varianceValues(:,10:18,1)];

hinton(varianceValues(:,:,1))

