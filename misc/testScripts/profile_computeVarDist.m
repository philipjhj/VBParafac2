load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')
%%


 profile on
tic;
myModel=varBayesModelParafac2(Y(:,:,:),10);



myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'};

myModel.computeVarDistribution;


myprofile = profile('info');

outdir=strcat('output/profiles/computeVarDist/',datestr(now),'_DIM',sprintf('_%d',size(myModel.data.X)),num2str(myModel.data.M,'_%d'));

profsave(myprofile,outdir)

% use profview(0,myprofile) to read results
% 
% profview(0,myprofile)