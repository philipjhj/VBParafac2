%load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
% set_wd(2)
% set(0,'DefaultFigureWindowStyle','docked')
int_No=19;
load(['/run/user/1001/gvfs/sftp:host=transfer.gbar.dtu.dk,user=phav/work1/phav/VBParafac2paper/data/dataBro/Models and data/Apple data/Int',num2str(int_No),'.mat'])
myInt=I19;
figure; for k=1:36; plot(squeeze(myInt(:,:,k))); hold on; end; axis tight
title(['Int ',num2str(int_No)])
%%

warning on MATLAB:nearlySingularMatrix

%%

load('/media/data/DataAndResults/VBParafac2paper/data/aminoAcid/amino.mat')
data = dataClass;
data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);
Mesti=3;
%%
% myModel=varBayesModelParafac2;

I=450;
J=9550;
K=100;
M=30;
Mesti = M;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
% options.initMethod = 'generative';
options.congruence = 0.4;
% 1e4 1e-3 i ARD tests
options.precision = [1e2 1e-6];
options.SNR = 2;
options.noiseType = 'heteroscedastic';
% [1e4 1e-8] creates problems for qC

% sumSNR = 0;
% for k = 1:100
rng('default')
rng('shuffle')
% rng(1)
data = varBayesModelParafac2.generateDataFromModel(options);

%
% sumSNR = sumSNR+10*log10(norm(data.Xtrue(:),'fro')^2/norm(data.Etrue(:),'fro')^2);
% end
% sumSNR/100

%%

data=dataClass;

data.Xunfolded=rand(450,9550,8);


%%
myInt=I21;
myInt = myInt/(norm(myInt(:),'fro')/numel(myInt));
data=dataClass;
data.Xunfolded = permute(myInt,[2 1 3]);
Mesti=4
%
% %
%%
% data=dataClass;
% data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);

normalModel = normalParafac2(data.X);
% normalModel = normalParafac2(permute(I2,[2 1 3]));
normalModel=normalModel.fitParafac2(Mesti);
[f1_DF,f2_DF]=normalModel.Parafac2Fit(data.Xtrue)

%% Plot CFtPt for comparison of models
rootFolder='/home/philipjhj/gbar_work1/VBParafac2paper/output/b722493/';
dataSetPath='rasmusBroData2/WineSingle'; Mtrue=5;
% dataSetPath='aminoAcid/amino'; Mtrue=3;
clf; 
titles={'Direct Fit','vMF VB Homoscedastic','cMN VB Homoscedastic','vMF VB Heteroscedastic','cMN VB Heteroscedastic'};
legendNames={'Comp. 1','Comp. 2','Comp. 3','Comp. 4','Comp. 5', 'Comp. 6', 'Comp. 7', 'Comp. 8'};

cols=distinguishable_colors(8);


% Get ground truth reconstruction
load([rootFolder,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath,'/trainedModels/DirectFit_',num2str(Mtrue),'_1.mat']);
normalModel=myModel;
[A,C,~,~,FPk] = normalizeParafac2(normalModel.A,normalModel.C,normalModel.F,normalModel.P);
CFPk = bsxfun(@times,permute(C,[2 3 1]),FPk);
T=size(A,2);
X_recon=zeros(prod(size(data.X)),T);
for t=1:T
    Xtemp=mtimesx(A(:,t),CFPk(t,:,:));
    X_recon(:,t) = Xtemp(:);
end
X_reconTrue=X_recon;


Ms=2:8;
% close all
% margin between above/below, left/right, and margin to edge
[ha, pos] = tight_subplot(numel(Ms),5,[.005 0.05],[.05 .05],[.1 .01]);
noiseConfig={'avgShared','avg'};

for M=Ms
%     figure(M-1)
load([rootFolder,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath,'/trainedModels/DirectFit_',num2str(M),'_1.mat']);
normalModel=myModel;
maxELBO=-realmax;
k=0;

for noiseIDX=1:2
% Test which trained model has highest ELBO
VBModelPath=['trainVBParafac2/finalResultsELBOModelOrder/',dataSetPath,'/trainedModels/mle_vonmises_max_0_',noiseConfig{noiseIDX},'_50_'];
for i = 1:5
load([rootFolder,VBModelPath,num2str(M),'_',num2str(i),'.mat'])
if myModel.qDistTrain.ELBO > maxELBO
    maxELBO = myModel.qDistTrain.ELBO;
    k=i;
end
end
k
load([rootFolder,VBModelPath,num2str(M),'_',num2str(k),'.mat'])
VBModel_vMF{noiseIDX}=myModel;


VBModelPath=['trainVBParafac2/finalResultsELBOModelOrder/',dataSetPath,'/trainedModels/mle_parafac2svd_max_0_',noiseConfig{noiseIDX},'_50_'];
for i = 1:5
load([rootFolder,VBModelPath,num2str(M),'_',num2str(i),'.mat'])
if myModel.qDistTrain.ELBO > maxELBO
    maxELBO = myModel.qDistTrain.ELBO;
    k=i;
end
end
k
load([rootFolder,VBModelPath,num2str(M),'_',num2str(k),'.mat'])
VBModel_cMN{noiseIDX}=myModel;

end

clear plts;
Acell={normalModel.A,VBModel_vMF{1}.qDistTrain.qA.mean,VBModel_cMN{1}.qDistTrain.qA.mean,VBModel_vMF{2}.qDistTrain.qA.mean,VBModel_cMN{2}.qDistTrain.qA.mean};
Ccell={normalModel.C,VBModel_vMF{1}.qDistTrain.qC.mean,VBModel_cMN{1}.qDistTrain.qC.mean,VBModel_vMF{2}.qDistTrain.qC.mean,VBModel_cMN{2}.qDistTrain.qC.mean};
Fcell={normalModel.F,VBModel_vMF{1}.qDistTrain.qF.mean,VBModel_cMN{1}.qDistTrain.qF.mean,VBModel_vMF{2}.qDistTrain.qF.mean,VBModel_cMN{2}.qDistTrain.qF.mean};
Pcell={normalModel.P,VBModel_vMF{1}.qDistTrain.qP.mean,VBModel_cMN{1}.qDistTrain.qP.mean,VBModel_vMF{2}.qDistTrain.qP.mean,VBModel_cMN{2}.qDistTrain.qP.mean};

X_reconCell=cell(1,5);
CFPkcell=cell(1,5);

for h=1:5
[A,C,~,~,FPk] = normalizeParafac2(Acell{h},Ccell{h},Fcell{h},Pcell{h});
CFPk = bsxfun(@times,permute(C,[2 3 1]),FPk);
CFPkcell{h}=CFPk;
T=size(A,2);
X_recon=zeros(prod(size(data.X)),T);
for t=1:T
    Xtemp=mtimesx(A(:,t),CFPk(t,:,:));
    X_recon(:,t) = Xtemp(:);
end
X_reconCell{h}=X_recon;

end
disp('1')

idx_cols=cell(1,5);

for g = 1:numel(idx_cols)
CC1=corrcoef([X_reconTrue X_reconCell{g}]);
CC1=abs(CC1(1:Mtrue,Mtrue+1:end));
measure=CC1;
[~, idx1, ~,~,idx1_rows] = greedy_component_match(measure);
idx1=idx1(~isnan(idx1));
idx_rows{g}=idx1_rows;
idx_cols{g}=idx1;
end
% 
% CC2=corrcoef([X_reconTrue X_reconCell{2}]);
% CC2=abs(CC2(1:Mtrue,Mtrue+1:end));
% measure=CC2;
% [~, idx2, ~,~,idx2_rows] = greedy_component_match(measure);
% idx2=idx2(~isnan(idx2));
% idx_rows{2}=idx2_rows;
% idx_cols{2}=idx2;
% 
% CC3=corrcoef([X_reconTrue X_reconCell{3}]);
% CC3=abs(CC3(1:Mtrue,Mtrue+1:end));
% measure=CC3;
% [~, idx3, ~,~,idx3_rows] = greedy_component_match(measure);
% idx3=idx3(~isnan(idx3));
% idx_rows{3}=idx3_rows;
% idx_cols{3}=idx3;

% Correlation plots
% subplot(2,3,sub2ind([3 2],1,2))
% imagesc(CC1(:,idx1))
% 
% colorbar; xticks(1:T); yticks(1:Mtrue)
% 
% subplot(2,3,sub2ind([3 2],2,2))
% imagesc(CC2(:,idx2))
% 
% colorbar; xticks(1:T); yticks(1:Mtrue)
% 
% subplot(2,3,sub2ind([3 2],3,2))
% imagesc(CC3(:,idx3))

% colorbar; xticks(1:T); yticks(1:Mtrue)



% set(ha(1:4),'XTickLabel',''); 

for h = 1:5
%     subplot(numel(Ms),3,sub2ind([3 numel(Ms)],h,M-(min(Ms)-1)))
    axes(ha(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))))
    
    CFPk=CFPkcell{h};
    
%squeeze(CFPk(m,:,:)
colorOrder=[idx_rows{h}(1:min(M,Mtrue)) numel(idx_rows{h})+1:M];
for m=1:M
%     keyboard
    CFPkm=squeeze(CFPk(idx_cols{h}(m),:,:));
    
    % Flip
    [v1,idx_flip1]=max(abs(CFPkm)); 
    for Nidx = 1:numel(idx_flip1)
       CFPkm(:,Nidx)=CFPkm(:,Nidx)*sign(CFPkm(idx_flip1(Nidx),Nidx));
    end
    
    plts(:,idx_cols{h}(m))=plot(CFPkm,'Color',cols(colorOrder(m),:));
    hold on
end


plotScale=0.85;
set(gcf,'Position',[2600 720 plotScale*1000 plotScale*1400]);
set(gca,'tickdir','out');
set(gca,'fontsize',12);
axis tight

if M==Ms(1)
title(titles{h})
elseif M~=Ms(end)
    set(ha,'XTickLabel','')
end

if h==1
    ylabel(['R=',num2str(M)])
end

% legend(plts(1,:),legendNames(1:M))
end
end

%end

%%
% for m = 2:5
    rng('default')
%normalModel.fitParafac2(3)
%
% normalModel.CCDParafac2
%  end
%
%
% parafac2(data.Xtrue,Mesti,[0 0],[0 0 0 0 0]);
[f1,f2]=normalModel.Parafac2Fit(data.Xtrue)
%
%%
Mesti=100;
% rng('default')
myModel=varBayesModelParafac2(data,Mesti);
%
% size(myModel.data.X)
%
% myModel=varBayesModelParafac2(Y,100);

myModel.opts.verbose = 1;
myModel.opts.debugFlag = 0;
myModel.opts.estimationP= 'parafac2svd';
% myModel.opts.estimationP = 'vonmises';
myModel.opts.estimationARD = 'max';
myModel.opts.estimationNoise = 'avg';
myModel.opts.initMethod = 'random';
myModel.opts.noiseLearningDelay=50;
myModel.opts.scaleLearningDelay=0;

myModel.opts.matrixProductPrSlab = 'mtimesx';
myModel.opts.nActiveComponents = 'threshold';
myModel.opts.showIter = 1;
% myModel.opts.rngInput = 15; % bad ones: 8, good ones: 15, 16, 19, 21, CORRECT: 31
myModel.opts.maxIter =10;
myModel.opts.maxTime = 4;

myModel.opts.activeParams = {'qC','qP','qA','qF','qAlpha','qSigma'};
%%
rng('default')
% rng(1)
myModel.partitionData(myModel.fullData.X)
% tic
myModel.fitTrainingData;
% toc 
%
[f1,f2]=myModel.Parafac2Fit(myModel.qDistTrain,data.Xtrue)

% myModel.opts.activeParams = {'qA','qC','qP','qSigma','qF'};
%%
Ms = 5:6
rng('default')
myModel.crossValidateM(Ms)
%%

cellfun( @(S) S.Data.stopReason, myModel.cvRunsTrain, 'uni', false )

ELBOS = cell2mat(cellfun( @(S) S.Data.ELBO, myModel.cvRunsTrain, 'uni', false ));
ELBOS-max(ELBOS,[],2)

%%
best_ELBO_true = squeeze(mean(max(myModel.CV_ELBOS(:,:,3),[],2)));
squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true
(squeeze(mean(max(myModel.CV_ELBOS,[],2)))-best_ELBO_true)/best_ELBO_true

%%
best_ELBO_all = squeeze(mean(max(myModel.CV_ELBOS_test,[],2)));
plot(Ms,best_ELBO_all,'o-')
xlabel('# of components')
ylabel('ELBO')
grid on
% figname=input('Write name of figure: ','s');
% saveas(gcf,strcat('output/paper/temps/',figname,'.jpg'))

%%
tic
myModel.fitTrainingData;
toc
%myModel.qDist.SNR
% myModel.Parafac2Fit
%%
% return
for k=1:myModel{3}{1}.fullData.K
    
clf

myModel{2}{1}.plotSolutionSynthK(k,0)
% myModel{3}{1}.Parafac2Fit(myModel{3}{1}.qDistTrain)
pause
end
%%
clf
myModel.plotELBO([100 myModel.data.iter-1])


%%
for k = 1:4
    
    clf
    myModel.plotSolutionRealMatrixK(k)
   keyboard
   
end
%%
figure
subplot(2,1,1)
imagesc(myModel.qDist.qC.mean)
colorbar
subplot(2,1,2)
imagesc(1./myModel.qDist.qAlpha.mean) % Lav som procent varians i stedet
colorbar
%%
m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
mask = m.mask;
set(0,'DefaultFigureWindowStyle','docked')
%%
close all
k = 1;
m = 1:25;

% close all
myModel.plotSolutionReal3D(k,m,mask)


%%

figure(2)
plot(nonzeros(myModel.ELBO_chain));

%%
myModel.restartqDist
myModel.computeVarDistribution


%%
k=2;
U = myModel.qDist.qA.mean*diag(myModel.qDist.qC.mean(k,:));
V = myModel.qDist.qP.mean(:,:,k)*myModel.qDist.qF.mean;

m=5;
plotComponent(U(:,1:m),V(:,1:m)', mask, [ 53    63    46])



%%

[A,F,C,P,fit]=parafac2(myModel.data.X,myModel.data.M,[0 0]);
%%
clf
plotParafac2SolutionK(1,myModel.data.X,A,C,F,cat(3,P{:}),myModel.data.Atrue,...
    myModel.data.Ctrue,myModel.data.Ftrue,myModel.data.Ptrue)
