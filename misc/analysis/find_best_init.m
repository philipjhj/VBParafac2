
myDir = '/media/data/DataAndResults/Thesis/output/results/';
files = dir(myDir);

%%

n_files = numel({files.name});
i=1;

max_ELBO = zeros(1,n_files);

cc = hsv(n_files);
for file = files'
    if numel(file.name) > 5
    load(strcat(myDir,file.name));
    end
    max_ELBO(i) = max(myModel.ELBO_chain);
    
    i=i+1
end

%%
