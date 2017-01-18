myDir = '/media/data/DataAndResults/Thesis/output/results/results_ARD_tests2';
files=dir(myDir);

% addpath(genpath('../../VBParafac2'))

%%


for i = 1:numel(files)
   if regexp(files(i).name,'.*__.*\.mat')
%     disp(files(i).name)
    
    load(strcat(myDir,'/',files(i).name))
    
    
    
    fullname = regexp(files(i).name,'(.*)(__.*)','tokens');
    
    FtF = myModel.data.Ftrue'*myModel.data.Ftrue;
    
    cong = 10*FtF(1,2);
    
    cong = ['0' num2str(cong)];
    
    
    newname=strcat(fullname{1}{1},'_',cong,fullname{1}{2});
%     disp(newname)
    
%     disp(['mv ',myDir,'/',files(i).name,' ',myDir,'/',newname])
    system(['mv ',myDir,'/',files(i).name,' ',myDir,'/',newname])

    
   end
    
end