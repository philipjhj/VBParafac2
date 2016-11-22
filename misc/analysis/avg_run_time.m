


% testDir = '/media/data/DataAndResults/Thesis/output/results/results_RealData_tests/';
testDir = '/media/data/DataAndResults/Thesis/output/results/ARD_tests/';
files=dir(testDir);

evaltime = [];

for file = files'
    if regexp(file.name, regexptranslate('wildcard','Int*'))
        %             disp(file.name)
        load(strcat(file.folder,'/',file.name))
       
        evaltime = [evaltime myModel.evaltime(end)];
        
    end
end


%%


% formatsfpec = '%C%d%b%d%d'
% mytable = readtable('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/nCompTrue.csv','Delimiter',',')

nCompTrue.data = csvread('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/nCompTrue_numeric.csv',1,0)

nCompTrue.columns = {'dataset','interval','badfit','nCompTrue1','nCompTrue2'}
nCompTrue.DataSets = {'Apple','Cheese','Wine','Aroma'}



