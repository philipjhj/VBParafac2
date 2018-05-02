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

clear
% data=dataClass;
% data.Xunfolded = permute(reshape(X,DimX),[2 3 1]);

normalModel = normalParafac2(data.X);
% normalModel = normalParafac2(permute(I2,[2 1 3]));
normalModel=normalModel.fitParafac2(Mesti);
[f1_DF,f2_DF]=normalModel.Parafac2Fit(data.Xtrue)

%% Plot CFtPt for comparison of models
clear; close all; colormap(flipud(colormap('gray')))

%localPath='/media/data/Dropbox/Studie/PhD/publications/VBParafac2/VBParafac2paper/';
localPath='';
savePath=[localPath,'paper/final_results/Reconstructions/'];

% localPath='/home/philipjhj/gbar_';
localPath='/';
rootFolder=[localPath,'work1/VBParafac2paper/output/b722493/'];
rootFolder_directfit=[localPath,'work1/VBParafac2paper/output/fb14caa/'];

% Uncomment to analyze wine data
dataSetName='Wine';
dataSetPath='rasmusBroData2/WineSingle';
dataSetPath_directfit='rasmusBroDataWine/WineFull'; 
Mtrue=5;
max_alpha=0.35;
alpha_min_values=linspace(0.30,0.15,7);

% Uncomment to analyze amino data
% dataSetName='amino';
% dataSetPath=['aminoAcid/',dataSetName]; Mtrue=3;
% dataSetPath_directfit=dataSetPath;
% max_alpha=0.7;
% alpha_min_values=linspace(max_alpha-0.05,0.5,7);

clf; 
titles={'Direct Fit','vMF VB Homoscedastic','cMN VB Homoscedastic','vMF VB Heteroscedastic','cMN VB Heteroscedastic'};
legendNames={'Comp. 1','Comp. 2','Comp. 3','Comp. 4','Comp. 5', 'Comp. 6', 'Comp. 7', 'Comp. 8'};

%cols=distinguishable_colors(8);
cols=linspecer(8,'qualitative');



Ms=2:8;
% close all
% margin between above/below, left/right, and margin to edge

figure(1); clf
[ha, pos] = tight_subplot(numel(Ms),5,[.005 0.05],[.05 .05],[.1 .01]);
%[haa, poss] = tight_subplot(numel(Ms),5,[.005 0.05],[.05 .05],[.1 .01]);
haa=gobjects(numel(Ms)*5,1);


figure(2); clf
[ha2, pos2] = tight_subplot(numel(Ms)-1,5,[.005 0.05],[.05 .05],[.1 .01]);

noiseConfig={'avgShared','avg'};


ylabh = cell(1,length(Ms));
ylab_pos = zeros(length(Ms),3);


% Get ground truth reconstruction

model_fits = [];
for modelID = 1:30
    load([rootFolder_directfit,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath_directfit,'/trainedModels/DirectFit_',num2str(Mtrue),'_',num2str(modelID),'.mat']);
    model_fits = [model_fits myModel.fit];
end
    [~,idx_mf]=max(model_fits);
%     dbstop if isempty(model_fits)
    load([rootFolder_directfit,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath_directfit,'/trainedModels/DirectFit_',num2str(Mtrue),'_',num2str(idx_mf),'.mat']);
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


for M=Ms
% %     figure(M-1)
model_fits = [];
for modelID = 1:30
    load([rootFolder_directfit,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath_directfit,'/trainedModels/DirectFit_',num2str(M),'_',num2str(modelID),'.mat']);
    model_fits = [model_fits myModel.fit];
end
    [~,idx_mf]=max(model_fits);
%     dbstop if isempty(model_fits)
    load([rootFolder_directfit,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath_directfit,'/trainedModels/DirectFit_',num2str(M),'_',num2str(idx_mf),'.mat']);
%     load([rootFolder_directfit,'trainDirectFitParafac2/finalResultsELBOModelOrderDirectFit/',dataSetPath,'/trainedModels/DirectFit_',num2str(M),'_',num2str(1),'.mat']);
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

load([rootFolder,VBModelPath,num2str(M),'_',num2str(k),'.mat'])
VBModel_cMN{noiseIDX}=myModel;

end

clear plts;
Acell={normalModel.A,VBModel_vMF{1}.qDistTrain.qA.mean,VBModel_cMN{1}.qDistTrain.qA.mean,VBModel_vMF{2}.qDistTrain.qA.mean,VBModel_cMN{2}.qDistTrain.qA.mean};
Ccell={normalModel.C,VBModel_vMF{1}.qDistTrain.qC.mean,VBModel_cMN{1}.qDistTrain.qC.mean,VBModel_vMF{2}.qDistTrain.qC.mean,VBModel_cMN{2}.qDistTrain.qC.mean};
Fcell={normalModel.F,VBModel_vMF{1}.qDistTrain.qF.mean,VBModel_cMN{1}.qDistTrain.qF.mean,VBModel_vMF{2}.qDistTrain.qF.mean,VBModel_cMN{2}.qDistTrain.qF.mean};
Pcell={normalModel.P,VBModel_vMF{1}.qDistTrain.qP.mean,VBModel_cMN{1}.qDistTrain.qP.mean,VBModel_vMF{2}.qDistTrain.qP.mean,VBModel_cMN{2}.qDistTrain.qP.mean};

if exist('X_reconCell','var')
    X_reconCell_old=X_reconCell;
end

n_Models=5;
X_reconCell=cell(1,n_Models);
CFPkcell=cell(1,n_Models);



for h=1:n_Models
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

idx_cols=cell(1,n_Models);
idx_rows=cell(1,n_Models);
allCorrelations=cell(1,n_Models);
for g = 1:n_Models
CC1=corrcoef([X_reconTrue X_reconCell{g}]);
CC1=abs(CC1(1:Mtrue,Mtrue+1:end));
measure=CC1;
allCorrelations{g}=measure;
[~, idx1, ~,~,idx1_rows] = greedy_component_match(measure);
idx1=idx1(~isnan(idx1));
idx_rows{g}=idx1_rows;
idx_cols{g}=idx1;
end


% Compute correlation between components of model m and m-1
if exist('X_reconCell_old','var')
    idx_cols_old=cell(1,n_Models);
    idx_rows_old=cell(1,n_Models);
    allCorrelations_old=cell(1,n_Models);
    Mold=M-1;
    for g = 1:n_Models
        CC1=corrcoef([X_reconCell_old{g} X_reconCell{g}]);
        CC1=abs(CC1(1:Mold,Mold+1:end));
        measure=CC1;
        allCorrelations_old{g}=measure;
        [~, idx1_old, ~,~,idx1_rows_old] = greedy_component_match(measure);
        idx1_old=idx1_old(~isnan(idx1_old));
        idx_rows_old{g}=idx1_rows_old;
        idx_cols_old{g}=idx1_old;
    end
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
%figure(1)


for h = 1:5
%     subplot(numel(Ms),3,sub2ind([3 numel(Ms)],h,M-(min(Ms)-1)))
    axes(ha(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))));
    
    CFPk=CFPkcell{h};
    
%squeeze(CFPk(m,:,:)
colorOrder=[idx_rows{h}(1:min(M,Mtrue)) numel(idx_rows{h})+1:M];
% 
% for m=1:M
% CFPkm=squeeze(CFPk(m,:,:));
% [v1,idx_flip1]=max(abs(CFPkm)); 
%     for Nidx = 1:numel(idx_flip1)
%        CFPk(m,:,Nidx)=CFPkm(:,Nidx)*sign(CFPkm(idx_flip1(Nidx),Nidx));
%     end
% end



% Plot correlations as background image scaled with fillPropX/Y
fillPropX=0.45;
fillPropY=0.55;

margin=0.0025;

if exist('allCorrelations_old','var')
    axes(ha2(sub2ind([5 numel(Ms)-1],h,M-(min(Ms)-1)-1)));
    imgCC_old=allCorrelations_old{h}(idx_rows_old{h},idx_cols_old{h})';
    imagesc(imgCC_old);
    hold on
    
    
    
    if M==Ms(1)
    title(titles{h});
    elseif M~=Ms(end)
        set(ha,'XTickLabel','');
    end

    if h==1
        ylabel(['M=',num2str(M)]);
    end
    
end

currentAxes=sub2ind([5 numel(Ms)],h,M-(min(Ms)-1));
axes(ha(currentAxes));
imgCC=allCorrelations{h}(idx_rows{h},idx_cols{h})';
W=size(CFPk,2);
m=size(imgCC,2);
maxX = W*(1-fillPropX/(2*m));
x2=maxX;
x1=x2-(W*fillPropX-2*(W-x2));

n=size(imgCC,1);

yMax=max(max(max(CFPk)));
yMin=abs(min(min(min(CFPk))));
H=(yMax+yMin);

y1=H*(1-fillPropY/(2*n))-yMin;
y_margin=yMax-y1;
if yMin > yMax
    yAdjust=H*(1-fillPropY);
else
    yAdjust=0;
end

y2=y1-H*fillPropY+2*(yMax-y1);

y1=y1-yAdjust;
y2=y2-yAdjust;


imagesc('CData',imgCC,'XData',[x1*(1+margin) x2*(1-margin)],'YData',[y1*(1-margin) y2*(1+margin)],[0 1])

hold on
crnX=x1-W+x2;
lenX=W*fillPropX;
crnY=y2-y_margin;
lenY=H*fillPropY;
r=rectangle('Position',[crnX crnY lenX lenY],'EdgeColor','k','LineWidth',2);
alpha(r,.5);

% Plot invisible line to move above image to right position
%plot(1:W,'LineStyle','None')

% Compute yticks to add if same axes
D=max(1e-2,10^(1-ceil(log10(yMax))));
if min(min(min(CFPk)))>0
    yTickMin=ceil(yMin*D)/D;
else
    yTickMin=floor(min(min(min(CFPk)))*D)/D;
    %disp(yTickMin)
end
yTicks=linspace(yTickMin,round(yMax*D)/D,5);
%disp(yTicks)
% haa(currentAxes)=axes('position',get(ha(currentAxes),'position'));
%haa(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))).YLim=ha(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))).YLim;
for m=1:M
% %     keyboard
%     CFPkm=squeeze(CFPk(idx_cols{h}(m),:,:));
%     
%     % Flip
%     [v1,idx_flip1]=max(abs(CFPkm)); 
%     for Nidx = 1:numel(idx_flip1)
%        CFPkm(:,Nidx)=CFPkm(:,Nidx)*sign(CFPkm(idx_flip1(Nidx),Nidx));
%     end
%     

    %set(haa(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))),'position',get(ha(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))),'position'))
%     plts(:,idx_cols{h}(m))=plot(linspace(1,x2,W),squeeze(CFPk(idx_cols{h}(m),:,:)),'Color',cols(colorOrder(m),:),'LineWidth',1.7);
%     set(plts(:,idx_cols{h}(m)),'visible','off')
     
     %plot(1:(2*W),'LineStyle','None')
     %plts(:,idx_cols{h}(m))=plot(squeeze(CFPk(idx_cols{h}(m),:,:)),'Color',cols(colorOrder(m),:),'LineWidth',1.7);
      %hold on
     %a=plts(:,idx_cols{h}(m));
    
%      alpha_values=linspace(max_alpha,alpha_min_values(M-1),M);
%      for i = 1:length(a)
%          a(i).Color(4)=alpha_values(m);
%      end
% %     
%     
%     yticks(yTicks)
%     yticklabels(yTicks)
%     
%     width=maxX;
%     height=yMax;
    
%     ph=findall(myH,'type','patch');
%     set(ph,'FaceColor','none')
end
% set(haa(currentAxes),'XLim',[0 W]);
% set(haa(currentAxes),'YLim',[-yMin yMax]);
% set(haa(currentAxes),'color','none');
% set(ha(currentAxes),'XLim',[0 W]);
% set(ha(currentAxes),'YLim',[-yMin yMax]);

plotScale=0.85;
%set(gcf,'Position',[2600 720 plotScale*1000 plotScale*1400]);
%set(gca,'tickdir','out');
set(gca,'fontsize',12);
grid on
axis tight


if M==Ms(1)
title(titles{h})
elseif M~=Ms(end)
% set(haa(currentAxes),'XTickLabel','')
end

if h==1
    ylabh{M-1}=ylabel(['M=',num2str(M)]);
    ylab_pos(M-1,:)=get(ylabh{M-1},'position');
end

% legend(plts(1,:),legendNames(1:M))
end

%plot correlations

% for h = 1:5
%     %     subplot(numel(Ms),3,sub2ind([3 numel(Ms)],h,M-(min(Ms)-1)))
%     axes(ha(sub2ind([5 numel(Ms)],h,M-(min(Ms)-1))))
%     hold on
%     %     CFPk=CFPkcell{h};
%     
%     %squeeze(CFPk(m,:,:)
%     % colorOrder=[idx_rows{h}(1:min(M,Mtrue)) numel(idx_rows{h})+1:M];
%     
%     imagesc(allCorrelations{h}(idx_rows{h},idx_cols{h})')
%     
%     plotScale=0.85;
%     set(gcf,'Position',[2600 720 plotScale*1000 plotScale*1400]);
%     set(gca,'tickdir','out');
%     set(gca,'fontsize',12);
%     axis tight
%     
%     if M==Ms(1)
%         title(titles{h})
%     elseif M~=Ms(end)
%         set(ha,'XTickLabel','')
%     end
%     
%     if h==1
%         ylabel(['R=',num2str(M)])
%     end
%     
%     % legend(plts(1,:),legendNames(1:M))
% end
end



%% Save the figures

n_figures=2;
figure_names = {['third_mode_loadings_',dataSetName]
                ['model_order_match_correlations_',dataSetName]};
for i = 1:n_figures
    figure(i);
    fig = gcf;
    
        pb=pbaspect;
    % Set output size
    %fig.PaperUnits = 'inches';
    figsize=40;
    fig.PaperPosition = [0 0 figsize*pb(1) figsize*pb(2)];
    
    % Remove any margins on figure (thirdparty)
%     if i == 2
      tightfig;
%     end
    
    % Correct ylabel position
    for ii = 1:length(Ms)
        ylab_pos=get(ylabh{ii},'position');
    end
    [~,min_id]=min(ylab_pos(:,1));

    set(ylabh{min_id},'Units','normalized');
    max_pos = get(ylabh{min_id},'position');
    for ii = 1:length(Ms)
        set(ylabh{ii},'Units','normalized');
        local_pos = get(ylabh{ii},'position');
        shift=local_pos(1)-max_pos(1);
        set(ylabh{ii},'position',local_pos - [shift 0 0]);
    end

    

   
    print(gcf,[savePath,figure_names{i}],'-depsc')
    %saveas(gcf,[savePath,figure_names{i}],'depsc');
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
