

noiseType = {'homo','hetero'};
SNR=0:-4:-12;

allCCD = zeros(9,4,numel(noiseType),2);
for noiseTypeIdx = 1:2
    for mEsti=2:10
        for snrIdx=1:4
            load(strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFit_sim_data__',...
                'dim_150_150_30_4_04_SNR_',num2str(SNR(snrIdx)),'_noiseType_',noiseType{noiseTypeIdx},'_mEstimate_',num2str(mEsti),'_datasetRNG_1.mat'));
            allCCD(mEsti-1,snrIdx,noiseTypeIdx,1) = myModel{mEsti-1}.CCDParafac2;
            allCCD(mEsti-1,snrIdx,noiseTypeIdx,2) = myModel{mEsti-1}.Parafac2Fit;
        end
    end
end
%%
savepath = '/media/data/Dropbox/Studie/Kandidat/ThesisPhilipJHJ/code/output/paper/final/';

plotname = {'CCD','FIT'}
for i = 1:2
for noiseTypeIdx=1:2
    clf
% for snrIdx =1:4
figure(noiseTypeIdx)
set(gcf,'units','normalized','position',[1 0.5 .16 .25])
bar(2:10,allCCD(:,:,noiseTypeIdx,i))
xlim([1.5 10.5])
xlabel('$R$','Interpreter','latex')
filename = strcat(savepath,plotname{i},'_sim_data_noiseType_',noiseType{noiseTypeIdx},'.png');
            tightfig(gcf)
            saveas(gcf,filename)
% end
end
end
%%
colors = 'rgbm';
Ms=2:10;
for noiseIdx = 1:2
    figure(noiseIdx)
    hold on
    
 ha_CCD=plot(Ms,allCCD(:,:,noiseIdx,1),'d--','linewidth',2,...
                                'markersize',15);
                            addaxislabel(1,'CCD');
                            addaxis([],[],[min(min(allCCD(:,:,noiseIdx,2))) max(max(allCCD(:,:,noiseIdx,2)))])
                            addaxislabel(2,'Fit');
                            ha_fit=addaxisplot(Ms,allCCD(:,:,noiseIdx,2),2,'s:','linewidth',2,...
                                'markersize',15);
%                         yyaxis right

                        
                        

for SNRIdx = 1:4
    set(ha_CCD(SNRIdx),'color',colors(SNRIdx));
    set(ha_fit(SNRIdx),'color',colors(SNRIdx));
end

xlim([2 10])
grid on

end