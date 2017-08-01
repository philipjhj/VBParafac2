


load('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFit_real_data__AminoAcid_amino__mEsti_3.mat')
normalModel=myModel{2}

load('/media/data/DataAndResults/VBParafac2paper/amino.mat')
data = dataClass;
data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);
Mesti=3;

I_start=100;
I=50;
CC_all = cell(1,I);

for i = 1:I
    disp(i)
myModel=varBayesModelParafac2(data,Mesti);
% size(myModel.data.X)
%
% myModel=varBayesModelParafac2(Y,100);


% with maxNoARD, avgShared, parafac2svd, 2000 iterations, rng 201, ELBO=-5.69e+06
myModel.opts.verbose = 1;
myModel.opts.debugFlag = 0;
myModel.opts.estimationP= 'parafac2svd';
% myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'maxNoARD';
myModel.opts.estimationNoise = 'avgShared';
myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 100;


% ELBO with rng=16: -4.16e+08, ELBO with rng=31: -4.17e+08

% rng=31 identifies the correct components, but it has a lower ELBO than
% many of the bad solutions. Also, its iterations contains errors with a
% relative magnitude around 1e-12

% bad ones: 8, good ones: 15, 16, 19, 21, 56 CORRECT: 31, 986461368, 986471051
% myModel.opts.rngInput = 986471051;
myModel.opts.rngInput = 200+i;
myModel.opts.maxIter = 2000;
% myModel.opts.maxTime = 4;

myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};
%
rng('default')
rng(1)
myModel.partitionData(myModel.fullData.X)
% tic
myModel.fitTrainingData;
% toc



bestVBmodel=myModel;
% bestVBmodel=myModel{m}{idx};

CC=corrcoef([bestVBmodel.qDistTrain.qC.mean y]);
CC=abs(CC(1:3,4:end));

[~,IDXs]=max(CC);
CCsorted=CC(IDXs,:);

CC_all{i}.CC = CCsorted;
CC_all{i}.ELBO = bestVBmodel.qDistTrain.ELBO;
disp(CC_all{i}.ELBO)
end
%%

all_elbo=cellfun(@(x) x.ELBO, CC_all);

[values,idx]=sort(all_elbo,'descend');


%%
figure(1);
clf;

i_n=25;
j_n=2;

for i = 1:I
    subplot(i_n,j_n,i);
    %disp(idx(i))
    %disp(CC_all{idx(i)}.ELBO)
    CC=CC_all{idx(i)}.CC;

%     CCC=corrcoef([normalModel.C y]);
%     CCC=abs(CCC(1:3,4:end));
% 
%     [~,IDXs]=max(CCC);
%     CCCsorted=CCC(IDXs,:);





imagesc(CC);
title(sprintf('%.4e',CC_all{idx(i)}.ELBO))
set(get(gca,'title'),'Position',[-0.1 2.5 1.00011])
set(gca,'FontSize',10)
colorbar
% figure(2);
%subplot(2,1,2);
%imagesc(CCC);
%colorbar

yticks([])
xticks([])

%pause(5)
end



