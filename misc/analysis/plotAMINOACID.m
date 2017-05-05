

XreconComp=normalModel.plotComponents;

%% Direct fitting
close all
A =normalModel.A;
C = normalModel.C;
P=normalModel.P;
F=normalModel.F;


%% find best init


%%
A =myModel{2}{idx}.qDistTrain.qA.mean;
C = myModel{2}{idx}.qDistTrain.qC.mean;
P=myModel{2}{idx}.qDistTrain.qP.mean;
F=myModel{2}{idx}.qDistTrain.qF.mean;


%%
xEmission=250:450;
xExcitation=240:300;
figure
plot(xEmission,A)

    figure
for k =1:5
   subplot(5,1,k)
    plot(xExcitation,mtimesx(P(:,:,k),F)')
end
figure
plot(C)
%%
%load('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_real_data__AminoAcid_amino__pMethod_vonmises_mEsti_3.mat')
% 2 3 1
xEmission=250:450;
myModel{2}{idx}.plotComponents(xEmission);
%%

figure
normalModel.plotComponents(xEmission);


%%

% xticklabels(xEmission)
% xticks(0:20:201)
% ylabel('Excitation')

load('/media/data/DataAndResults/VBParafac2paper/results_paper/DirectFit_real_data__AminoAcid_amino__mEsti_3.mat')
normalModel=myModel{2}
%%
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_real_data__AminoAcid_amino__pMethod_parafac2svd_mEsti_3.mat')
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/ARD_real_data__AminoAcid_amino__pMethod_vonmises_mEsti_3.mat')
% close all
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/long_ARD_real_data__AminoAcid_amino__pMethod_vonmises_mEsti_3.mat')
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/long_ARD_real_data__AminoAcid_amino__pMethod_parafac2svd_mEsti_3.mat')
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/SharedScale_real_data__AminoAcid_amino__pMethod_parafac2svd_mEsti_3.mat')
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/SharedScale_real_data__AminoAcid_amino__pMethod_vonmises_mEsti_3.mat')
% load('/media/data/DataAndResults/VBParafac2paper/results_paper/CV_SharedScale_real_data__AminoAcid_amino__pMethod_vonmises_mEsti_3.mat')
m=3;
T=50
ELBO=[]
for t=1:T
    ELBO(t)=myModel{m}{t}.dataTrain.ELBO;
end
ELBO

[~,idx]=max(ELBO)
% diff(myModel{m}{idx}.dataTrain.ELBO_chain(end-1:end))
%%

bestVBmodel=myModel;
% bestVBmodel=myModel{m}{idx};

CC=corrcoef([bestVBmodel.qDistTrain.qC.mean y]);
CC=abs(CC(1:3,4:end));

[~,IDXs]=max(CC);
CCsorted=CC(IDXs,:);

CCC=corrcoef([normalModel.C y]);
CCC=abs(CCC(1:3,4:end));

[~,IDXs]=max(CCC);
CCCsorted=CCC(IDXs,:);

% figure(1);
subplot(2,1,1);
imagesc(CC); 
colorbar
% figure(2);
subplot(2,1,2);
imagesc(CCC);
colorbar



