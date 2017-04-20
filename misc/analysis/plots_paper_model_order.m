close all
warning('off','all')
savepath = '/media/data/Dropbox/Studie/Kandidat/ThesisPhilipJHJ/code/output/paper/final/';

plotType='ELBO'; %{'ELBO','Hinton'};

testName = {'ARD','CV'};
dataType = 'real';
% realData = 'Apple';
realData = 'AminoAcid_amino';

pMethod = {'parafac2svd','vonmises'};
colors = 'rgbm';
markers = 'ox';
linetypes={'-','--'};
noiseType = {'homo','hetero'};
SNR = 0:-4:-12;

Ms = 2:10;
M = numel(Ms);
T = 5;

intervalIdxs = [1:6 8 13 20 21 23 25 26];

if strcmp(dataType,'sim')
    S = numel(SNR);
    N = numel(noiseType);
    I = 1;
else
    S = 1;
    N = 1;
    if strcmp(realData,'Apple')
    I = numel(intervalIdxs);
    else
        I=1;
    end
    if strcmp(realData,'Apple')
        load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/nCompTrue.mat')
    end
end

for testIdx = 2:2
    for intervalIdx = 1:I
        clf
        allTestValues = zeros(M,numel(SNR),numel(noiseType),numel(pMethod));
        allCCDnFit = zeros(M,2,numel(SNR),numel(noiseType));
        alphaMeans = zeros(10,M,numel(SNR),numel(noiseType),numel(pMethod));
        for pMethodIdx = 1:2
            
            for noiseIdx = 1:N
                for SNRIdx = 1:S
                    ELBOs = zeros(M,T);
                    for m = 1:M
                        if strcmp(dataType,'sim')
                            filename = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/',testName{testIdx},'_',dataType,'_data__dim_150_150_30_4_04__pMethod_',pMethod{pMethodIdx},'_SNR_',...
                                num2str(SNR(SNRIdx)),'_noiseType_',noiseType{noiseIdx},'_mEstimate_',num2str(Ms(m)),'_datasetRNG_1.mat');
                            
                        elseif strcmp(dataType,'real')
                            if strcmp(realData,'Apple')
                            filename = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/',testName{testIdx},'_',dataType,'_data__Apple_Int',num2str(intervalIdxs(intervalIdx)),'__pMethod_',...
                                pMethod{pMethodIdx},'_mEsti_',num2str(Ms(m)),'.mat');
                            filenameDF = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFitting_real_data__Apple_Int',num2str(intervalIdxs(intervalIdx)),'_mEsti_',num2str(Ms(m)),'.mat');
                            
                            groundTruthNcomponents=nCompTrue.data(intervalIdxs(intervalIdx),4);
                            else
                                filename = strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/',testName{testIdx},'_',dataType,'_data__',realData,'__pMethod_',...
                                pMethod{pMethodIdx},'_mEsti_',num2str(Ms(m)),'.mat');
                            
                            filenameDF= strcat('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFit_real_data__AminoAcid_amino__mEsti_',num2str(Ms(m)),'.mat');
                            end
                            
                        end
                        if exist(filename, 'file') == 2
                            load(filename);
                        else
                            continue
                        end
                        
                        if strcmp('ARD',testName{testIdx})
                            for t = 1:T
                                ELBOs(m,t) = myModel{m}{t}.dataTrain.ELBO;
                                [v,i]=max(ELBOs,[],2);
                            end
                            
                            alphaMeans(1:(m+1),m,SNRIdx,noiseIdx,pMethodIdx) = myModel{m}{i(m)}.qDistTrain.qAlpha.mean';
                            
                            
                            
                        elseif strcmp('CV',testName{testIdx})
                            allTestValues(m,SNRIdx,noiseIdx,pMethodIdx) = mean(max(myModel{m}.CV_ELBOS_test,[],2));
                        end
                        
                        if strcmp(dataType,'real')
                            load(filenameDF);
                            allCCDnFit(m,1,SNRIdx,noiseIdx) = myModel{m}.CCDParafac2;
                            allCCDnFit(m,2,SNRIdx,noiseIdx) = myModel{m}.Parafac2Fit;
                        end
                    end
                    if strcmp('ARD',testName{testIdx})
                        [v,i]=max(ELBOs,[],2);
                    elseif strcmp('CV',testName{testIdx})
                        v = allTestValues(:,SNRIdx,noiseIdx,pMethodIdx);
                    end
                    
                    if strcmp(plotType,'ELBO')
                        figure(noiseIdx)
%                         set(gcf,'units','normalized','position',[1 0.5 .16 .25])
                        set(gcf,'units','normalized','position',[1 0.5 .60 .60])
%                         set(gca,'fontsize',32)
                        hold on
                        
                        
                        plot(Ms(v<0),v(v<0),linetypes{1},'Color',colors(SNRIdx),'Marker',markers(pMethodIdx),'lineWidth',2,...
                            'markersize',15);
                        if pMethodIdx==2 % Remove when Direct exist for AMINO
%                         yyaxis left
%                         bar(Ms,allCCDnFit(:,:,SNRIdx,noiseIdx))
                            set(gca,{'ycolor'},{'r'})
                            
                            addaxis(Ms,allCCDnFit(:,1,SNRIdx,noiseIdx),...
                                [min(0,min(allCCDnFit(:,1,SNRIdx,noiseIdx))-5) 100],'d:','linewidth',2,...
                                'markersize',15);
                            addaxis(Ms,allCCDnFit(:,2,SNRIdx,noiseIdx),...
                                [min(allCCDnFit(:,2,SNRIdx,noiseIdx))-0.1 100],'s:','linewidth',2,...
                                'markersize',15);
%                         yyaxis right
addaxislabel(1,'ELBO');
                        addaxislabel(2,'CCD');
                        addaxislabel(3,'Fit');
                        end
                        
                        
                        
                        if strcmp(dataType,'real') && strcmp(realData,'Apple')
                            
                            SP=groundTruthNcomponents;
                            line([SP SP],ylim,'Color','b','LineWidth',4,'linestyle','--')
                        end
                        %                 plt = Plot(); % create a Plot object and grab the current figure
                        
                        
                        grid on
                        
                        %                         if strcmp(noiseType{noiseIdx},'homo')
                        %                             ylabel('ELBO')
                        
                        %                     plt.YLabel = 'ELBO';
                        %                         end
                        %                         if strcmp(testName{testIdx},'CV')
                        xlabel('$R$','Interpreter','latex')
                        %                     plt.XLabel = 'R'; % xlabel
                        %                         end
                        xlim([2 10])
                    end
                end
                %title(strcat(noiseType{noiseIdx},'scedastic'))
                if strcmp(plotType,'ELBO')
                    %             legend(strtrim(cellstr(num2str(SNR','SNR=%d'))))
                end
            end
        end
        
        %subplot(1,2,1)
        if strcmp(plotType,'ELBO') && strcmp(dataType,'sim')
            figure(1)
            ylimits(1,1:2)=ylim;
            %subplot(1,2,2)
            figure(2)
            ylimits(2,1:2)=ylim;
            
            ybot = min(ylimits(:,1));
            ytop = max(ylimits(:,2));
            
            %     subplot(1,2,1)
            figure(1)
            ylim([ybot ytop]);
            filename = strcat(savepath,testName{testIdx},'_',dataType,'_data_noiseType_',noiseType{1},'.png');
            saveas(gcf,filename)
            figure(2)
            %     subplot(1,2,2)
            ylim([ybot ytop]);
            filename = strcat(savepath,testName{testIdx},'_',dataType,'_data_noiseType_',noiseType{2},'.png');
            tightfig(gcf)
            saveas(gcf,filename)
        elseif strcmp(plotType,'ELBO') && strcmp(dataType,'real')
            if strcmp(realData,'Apple')
                filename = strcat(savepath,'realdata/',testName{testIdx},'_',dataType,'_data_Apple_Int',num2str(intervalIdxs(intervalIdx)),'.png');
            else
                filename = strcat(savepath,'realdata/',testName{testIdx},'_',dataType,'_data__',realData,'.png');
            end
            tightfig(gcf)
            saveas(gcf,filename)
        end
        
        if strcmp(plotType,'Hinton') && strcmp(testName{testIdx},'ARD')
            %% Hinton
            varianceValues = 1./alphaMeans;
            varianceValues(varianceValues==Inf) = 0;
            varianceValues = sort(varianceValues,1,'ascend');
            
            subpSize = [4 1];
            close all
            [ha, pos] = tight_subplot(subpSize(1),subpSize(2),[.01 .03],[.1 .01],[.01 .01]);
            for i = 1:subpSize(2)
                set(ha(i),'Visible','off')
            end
            
            
            for n = 1:N
                for p = 1:2
                    figure
%                     figure(sub2ind([2 2],p,n))
                    set(gcf,'units','normalized','position',[1 0.5 .16 .25],'InvertHardCopy','off')
                    set(gcf,'Color', [0.5 0.5 0.5],'Colormap', [0 0 0; 1 1 1])
                    if strcmp(dataType,'sim')
                        [ha, pos] = tight_subplot(subpSize(1),subpSize(2),[.05 .01],[.01 .01],[.01 .01]);
                        
                        for i = 1:subpSize(2)
                            set(ha(i),'Visible','off')
                        end
                    end
                    for s = 1:S
                        
                        
                        title=strcat(testName{testIdx},'_',dataType,'_data_pMethod_',pMethod{p},'_noiseType_',noiseType{n});
%                         set(ha(s),'Visible','off')
%                         hinton(varianceValues(sum(varianceValues(:,:,s,n,p),2)>1e-3,:,s,n,p),ha(5-s),[],title);
                        hinton(varianceValues(1:end,2:end,s,n,p),s,subpSize,title);
                        %                         hinton(varianceValues(:,:,s,n,p));
                    end
                    tightfig(gcf)
                    if strcmp(dataType,'sim')
                        filename = strcat(savepath,'hinton_',testName{testIdx},'_',dataType,'_data_pMethod_',pMethod{p},'_noiseType_',noiseType{n},'.png');
                    elseif strcmp(dataType,'real')
                        if strcmp(realData,'Apple')
                        filename = strcat(savepath,'realdata/hinton_',testName{testIdx},'_',dataType,'_data_Apple_Int_',num2str(intervalIdxs(intervalIdx)),'_pMethod_',pMethod{p},'.png');
                    else
                        filename = strcat(savepath,'realdata/hinton_',testName{testIdx},'_',dataType,'_data_',realData,'_',num2str(intervalIdxs(intervalIdx)),'_pMethod_',pMethod{p},'.png');
                        end
                    end
                    saveas(gcf,filename)
                end
            end
            
            
        end
    end
end
warning('on','all')