
profile on

myModel = varBayesModelParafac2;



myModel.qDist.debugflag = 0;
myModel.verbose = 1;
% myModel.qDist.method = 'vonmises';
myModel.qDist.method = 'parafac2svd';

myModel.qDist.activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'};

myModel.computeVarDistribution;

myprofile = profile('info');


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 1;
else
    poolsize = p.NumWorkers;
end

outdir=strcat('output/profiles/computeVarDist/',datestr(now),...
    '_DIM',num2str(size(myModel.data.X),'_%d'),num2str(myModel.data.M,'_%d'),...
    '_nWorkers_',num2str(poolsize));
profsave(myprofile,outdir)

% use profview(0,myprofile) to read results
% 
% profview(0,myprofile)