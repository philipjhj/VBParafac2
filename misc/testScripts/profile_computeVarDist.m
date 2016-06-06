
profile on

myModel = varBayesModelParafac2;

myModel.computeVarDistribution



%gibbs = gibbs;

myprofile = profile('info');

outdir='output/profiles/computeVarDist/'
save(strcat(outdir,'test'),'myprofile')

%end

% use profview(0,myprofile) to read results

profview(0,myprofile)