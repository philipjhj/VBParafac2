
colors = 'rgbm';
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;
for noiseIdx = 1:2
for SNRIdx = 1:numel(SNR)
load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_sim_data__dim_150_150_30_4_04__pMethod_vonmises_SNR_',...
    num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat'))

% FIND BEST INITIALIZATION

Ms = 2:10;
M = numel(Ms);
ELBOs = zeros(M,T);
for m = 1:M
    T = 5;
    for t = 1:T
        ELBOs(m,t) = myModel{m}{t}.dataTrain.ELBO;
%         disp(ELBOs(m,t))
%         disp(t)
%         disp(m)
    end
    
end

[v,i]=max(ELBOs,[],2);


subplot(1,2,noiseIdx)
hold on
plot(Ms,v,'o-','Color',colors(SNRIdx))
grid on
legend(strtrim(cellstr(num2str(SNR'))))


end
title(noiseType{noiseIdx})
end
% plot(myModel{m}{t}.dataTrain.n_components_hard)