


% [val_fit,idx_fit] = max(myAnalysis.testResults(:,:,:,:,:,:,:),[],ndims(myAnalysis.testResults));

[val_ELBO,idx_ELBO] = max(myAnalysis.ELBOArray,[],5);

%%

% [i1,i2,i3,i4,i5]=ind2sub(size(myAnalysis.ELBOArray),find(ismember(myAnalysis.ELBOArray(:),nonzeros(val_ELBO(:)))));

idx_max_ELBO = find(ismember(myAnalysis.ELBOArray(:),nonzeros(val_ELBO(:))));

%%

fit1 = myAnalysis.testResults(1,:,:,:,:,:);
fit1 = fit1(idx_max_ELBO);

fit1_final = zeros(size(myAnalysis.ELBOArray));
fit1_final(idx_max_ELBO) = fit1;
fit1_final = sum(fit1_final,5);
fit1_final = sum(fit1_final,4)./sum(fit1_final~=0,4);
%%
fit2 = myAnalysis.testResults(1,:,:,:,:,:);
fit2 = fit2(idx_max_ELBO);

fit2_final = zeros(size(myAnalysis.ELBOArray));
fit2_final(idx_max_ELBO) = fit2;
fit2_final = sum(fit2_final,5);
fit2_final = sum(fit2_final,4)./sum(fit2_final~=0,4);


%%
% avg_fit = sum(val_fit,5)./sum(val_fit~=0,5);
avg_fit_ELBO = sum(val_fit,5)./sum(val_fit~=0,5);

%%

% avg_fit = myAnalysis.testResults(:,:,:,:,7,3);
avg_fit = fit1_final;
% avg_fit = fit2_final;
figure(2)
clf
for pMethod = 1
    for ARD = [1 3]
%         for fitvalue = 1
        dataMethod = avg_fit(pMethod,:,ARD);
        dataMethod = reshape(dataMethod,size(dataMethod,1),size(dataMethod,2)*size(dataMethod,3));
        
        plot(2:7,dataMethod')
        hold on
        axis tight
        
        %     plt.Legend=ARDMethods;
%         end
    end
end

plt = Plot();
        plt.BoxDim = [8, 8];
        plt.LineStyle = {'-.','-.'};
        plt.Colors = {'r','b','b','b','g','g'};
        plt.YLim = [70 100];


%%




idx_max_ELBO = find(ismember(myAnalysis.ELBOArray(:),nonzeros(val_ELBO(:))));

% [i1,i2,i3,i4,i5]=ind2sub(size(myAnalysis.ELBOArray),find(ismember(myAnalysis.ELBOArray(:),nonzeros(val_ELBO(:)))));



%%

[~,best_run] = max(myAnalysis.ELBOArray(:));

[i1,i2,i3,i4,i5]=ind2sub(size(myAnalysis.ELBOArray),best_run);

%%
disp(myAnalysis.testOpts.pMethod(i1))
disp(myAnalysis.testOpts.mEsti(i2))
disp(myAnalysis.testOpts.ARD(i3))
disp(myAnalysis.testOpts.datasetRNG(i4))
disp(myAnalysis.testOpts.initRNG(i5))


%     'parafac2svd'
% 
%     '5'
% 
%     'max'
% 
%     '8'
% 
%     '6'
        