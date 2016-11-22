function [meanCorr mci domfrac domseq fracunmatch]=calcperf(Sest,Strue)
%Strue and Sest both noc x voxels x subjects
mci=zeros(size(Strue,3),size(Strue,1)); %zeros(M,L)
matching=zeros(size(Strue,3),size(Strue,1)); %zeros(M,L)

for i=1:size(Sest,3) % for each subject
    C=corr(Strue(:,:,i)',Sest(:,:,i)'); %corrcoef(S_est S_true(mask))
    nic=size(Strue,1); %number of components
    %r=zeros(size(Si.ic,1),1);
    oc=1:size(Strue,1); % 1:L
    ec=1:size(Strue,1); % 1:L
    for j=1:size(Strue,1) % for each component
        [rj ridx]=max(abs(C)); % find max value in each col and their row index
        [rj cidx]=max(rj); % find max col index and highest correlation
        ridx=ridx(cidx);  % find max row index
        mci(i,oc(ridx))=rj; % mci = correlation of each component for each subject
        C(ridx,:)=[]; % remove row
        C(:,cidx)=[]; % remove col
        matching(i,oc(ridx))=ec(cidx); % components for each subject
        oc(ridx)=[];
        ec(cidx)=[];
    end
end

um=unique(matching,'rows'); % returns all unique rows
for i=1:size(um,1)
    % Estimate how many times each unique row is present
    nm(i)=sum(sum(abs(bsxfun(@minus,matching,um(i,:))),2)==0);
%sum(sum(abs(matching - unique row),col) == 0)
end
[moc midx]=max(nm); % find first max repetiton. moc = value, midx = index.
domfrac=moc/size(Strue,3); % max repetion/subjects = dominating fraction
domseq=um(midx,:); % one of the dominating sequences
disp(sprintf('mean correlation %.3f, dominating sequence %sat %.1f%%',...
    mean(mci(:)),sprintf('%i ',um(midx,:)),domfrac*100.));
% M - sum((matching - domseq) ==0) (how many times is each component the
% same as the component in domseq). 30 comoponents total - how many are
% matched to the domseq. so fraction of unmatched components.
fracunmatch=(size(Strue,3)-sum(bsxfun(@minus,matching,um(midx,:))==0))/size(Strue,3);
disp(sprintf('occurance of non-matching over components:%s',...
    sprintf(' %.1f%%',fracunmatch*100)));
meanCorr = mean(mci(:)); % mean correlation