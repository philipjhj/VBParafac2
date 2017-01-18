

myAnalysis = analysisVBParafac2([]);
%%


load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Fmax_it_models.mat')
load('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/Apple data/Int1.mat')

normalModel = normalParafac2(I1,A{1,5},C{1,5},H{1,5},P{1,5});
%%
myAnalysis = analysisVBParafac2([]);
saveFlag = 1;

% myAnalysis.plotReconElutionProfiles(myModel,'test',saveFlag)
% myAnalysis.plotReconElutionProfiles(normalModel,'NormalParafac2_M_10',saveFlag)
% myAnalysis.plotReconElutionProfiles(myModel,'VBParafac2_apple_int2_best',saveFlag)
myAnalysis.plotReconElutionProfiles(myModel,'Parafac2_apple_int2_true',saveFlag)

%%
close all


%%
myAnalysis = analysisVBParafac2([]);

resultsPath = '/media/data/DataAndResults/Thesis/output/results/';

% testTitle = 'results_RealData_tests';
testTitle = 'results_RealData_tests3';
% testTitle= 'results_ARD_tests';
% testTitle='results_ARD_tests2';
% testTitle= 'first_dones';
% testTitle= 'test_folder';

myAnalysis = myAnalysis.computeTableELBOAll(resultsPath,testTitle);

load gong.mat;
sound(y)

%%

myAnalysis.find_max_ELBO_lin_idx;
myAnalysis.find_max_ELBO_sub_idx;

myAnalysis.find_max_ELBO_filenames;
%%
myAnalysis.dispBestRuns
%%
myAnalysis=myAnalysis.computeTableResultsFull;

load gong.mat;
sound(y)
%%
myAnalysis.sortOrder = [1];
myAnalysis.computeTableFoundComponents;

%%
myAnalysis.fontsize = 12;   
myAnalysis.plotTableResults
% close all


%% TABLE OUTPUT

clear input;
fprintf('\n\nExample 3: using an array as data input\n\n');

% numeric values you want to tabulate:
% this field has to be an array or a MATLAB table
% in this example we use an array
input.data = [1.12345 NaN 3.12345; ...
    4.12345 5.12345 6.12345; ...
    7.12345 8.12345 9.12345; ...
    10.12345 11.12345 12.12345];

% Optional fields:

% Set column labels (use empty string for no label):
input.tableColLabels = {'col1','col2','col3'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'row1','row2','','row4'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.3f',2,'%.1f',1}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 1;

% Uses booktabs basic formating rules ('1' = using booktabs, '0' = not using booktabs). 
% Note that this option requires the booktabs package being available in your LaTex. 
% Also, setting the booktabs option to '1' overwrites input.tableBorders if it exists.
% input.booktabs = 0;


% LaTex table caption:
input.tableCaption = 'MyTableCaption';

% LaTex table label:
input.tableLabel = 'MyTableLabel';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latex = latexTable(input);

% myDir=

% save LaTex code as file
fid=fopen('MyLatex.tex','w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');





%% OLD CODE BELOW THIS POINT
imagesc((myAnalysis.nFoundComponents))
colorbar
axis image
yticks(1:36)


testConfig = sortrows(...
    unique(myAnalysis.max_ELBO_sub_idx(:,1:3),'rows'),myAnalysis.sortOrder);
ynames = cell(1,size(testConfig,1));

for i = 1:size(testConfig,1)
    
    myString = '';
    for j = 1:3 
        myString = sprintf('%s %s',myString,...
            myAnalysis.testOpts.(myAnalysis.testOpts_names{j}){...
            testConfig(i,j)});
    end
    ynames{i} = myString;
end

% ynames = unique(cat(1,ynames{:}));

set(gca,'YTickLabel',ynames)
xlabel('Data Sets')
%%





close all



%%
%%
myAnalysis.plotARDtest

%%

subplot(2,2,1)
imagesc(normalModel.A)
colorbar
title('Normal Parafac2 - A')
subplot(2,2,2); 
imagesc(myModel.qDist.qA.mean)
title('VB Parafac2 - A')
colorbar

subplot(2,2,3)
imagesc(normalModel.C)
colorbar
title('Normal Parafac2 - C')
subplot(2,2,4); 
imagesc(myModel.qDist.qC.mean)
colorbar
title('VB Parafac2 - C')






