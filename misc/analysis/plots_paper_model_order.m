

testName = 'ARD';
pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
markers = 'ox';
linetypes={'--','-'};
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;
for pMethodIdx = 1:2
for noiseIdx = 1:2
    for SNRIdx = 1:numel(SNR)
        filename = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/',testName,'_sim_data__dim_150_150_30_4_04__pMethod_',...
            pMethod{pMethodIdx},'_SNR_',num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat');
        if exist(filename, 'file') == 2
            load(filename);
        else
            continue
        end
%         load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_sim_data__dim_150_150_30_4_04__pMethod_parafac2svd_SNR_',...
%             num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_datasetRNG_1.mat'))
        
        % FIND BEST INITIALIZATION
        
        Ms = 2:10;
        M = numel(Ms);
        T = 5;
        ELBOs = zeros(M,T);
        for m = 1:M
            
            for t = 1:T
                ELBOs(m,t) = myModel{m}{t}.dataTrain.ELBO;
                %         disp(ELBOs(m,t))
                %         disp(t)
                %         disp(m)
            end
            
        end
        
        [v,i]=max(ELBOs,[],2);
        
%         figure(pMethodIdx)
%         subplot(1,2,noiseIdx)
        hold on
        plot(Ms,v,linetypes{noiseIdx},'Color',colors(SNRIdx),'Marker',markers(pMethodIdx))
        grid on
        
        xlabel('Number of components')
        ylabel('ELBO')
    end
    title(strcat(noiseType{noiseIdx},'scedastic'))
end
legend(strtrim(cellstr(num2str(SNR','SNR=%d'))))
% subplot(1,2,1)
% ylimits(1,1:2)=ylim;
% subplot(1,2,2)
% ylimits(2,1:2)=ylim;
% 
% ybot = min(ylimits(:,1));
% ytop = max(ylimits(:,2));
% 
% subplot(1,2,1)
% ylim([ybot ytop]);
% subplot(1,2,2)
% ylim([ybot ytop]);
end



% plot(myModel{m}{t}.dataTrain.n_components_hard)