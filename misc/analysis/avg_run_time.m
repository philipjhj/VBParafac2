


testDir = '/media/data/DataAndResults/Thesis/output/results/results_RealData_tests3/';
% testDir = '/media/data/DataAndResults/Thesis/output/results/results_ARD_final/';
files=dir(testDir);

evaltime = [];

not_converged = [];

for file = files'
%     if regexp(file.name, regexptranslate('wildcard','Int*'))
    if regexp(file.name,'.*__.*')

        %             disp(file.name)
        load(strcat(file.folder,'/',file.name))
       
        evaltime = [evaltime myModel.evaltime(end)];
        
        not_converged = [not_converged 1e-7<((myModel.ELBO_chain(end)-myModel.ELBO_chain(end-1))/myModel.ELBO_chain(end))];
        
        
    end
end
%%

mean(evaltime)/60/60
sum(not_converged)

%%


plot(evaltime/60/60)


%%


% formatsfpec = '%C%d%b%d%d'
% mytable = readtable('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/nCompTrue.csv','Delimiter',',')

nCompTrue.data = csvread('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/nCompTrue_numeric.csv',1,0)

nCompTrue.columns = {'dataset','interval','badfit','nCompTrue1','nCompTrue2'}
nCompTrue.DataSets = {'Apple','Cheese','Wine','Aroma'}



