% Variational Factor combinations

%load('matlab.mat')
import matlab.unittest.TestSuite
testDir = 'tests/';
myTests = TestSuite.fromClass(?qDistributionTest,'Tag','CAVI');
%
nFails=numel([results.Failed]);

nParams = 6;
nTestRepeats = 10;

param_mat = zeros(2^nParams,nParams);
k=1;
for i = 1:nTestRepeats:nTestRepeats*2^nParams
 
    param_mat(k,:) = [myTests(i).Parameterization(1:nParams).Value];
    k=k+1;
end

%
k=1;
counts = zeros(2^nParams,1);

for i = find([results.Failed])
    idx = ismember(param_mat,[myTests(i).Parameterization(1:nParams).Value],'rows');
    counts(idx) = counts(idx)+1;
end


%


paramNames={myTests(1).Parameterization(1:nParams).Property};

%imagesc(bsxfun(@times,counts,param_mat))
% subplot(1,8,2:8)

% imagesc(param_mat)
pcolor((10*[counts./max(counts) param_mat zeros(2^nParams,1); zeros(1,nParams+2)]))
% shading(gca,'faceted')
axis ij
axis square
set(gca,'XTick',1.5:7.5)
set(gca,'XTickLabel',['Fails',paramNames])
% colormap('hsv')
% subplot(1,8,1)
% imagesc(counts)
colorbar
% axis square
% xticklabels(paramNames)

% set(gca,'YTick',1:2^nParams)
% set(gca,'YTickLabel',counts)
% colorbar
%reshape([myTests(1).Parameterization.Value],[1 9])

%results.Name([1:10])


%%