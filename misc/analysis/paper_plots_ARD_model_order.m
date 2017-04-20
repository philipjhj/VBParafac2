
pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;

for pMethodIdx = 1:2
    for noiseIdx = 1:2
        for SNRIdx = 1:numel(SNR)
            Ms = 2:10;
            M = numel(Ms);
            T = 5;
            ELBOs = zeros(M,T);
            for m = 1:M
                
                load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_sim_data__dim_150_150_30_4_04__pMethod_',pMethod{pMethodIdx},'_SNR_',...
                    num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_mEstimate_',num2str(Ms(m)),'_datasetRNG_1.mat'))
                
                
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
            xlim([2 10])
            
        end
        title(noiseType{noiseIdx})
    end
end

% plot(myModel{m}{t}.dataTrain.n_components_hard)