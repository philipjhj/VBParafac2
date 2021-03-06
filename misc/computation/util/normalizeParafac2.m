function [A,C,F,P,FtPt] = normalizeParafac2(A,C,F,P)

M = size(A,2);
K = size(P,3);

FtPt = mtimesx(F',permute(P,[2 1 3]));
FtPtScale = max(abs(FtPt),[],2);

%FtPtScale(FtPtScale<eps) = 0;
FtPt = bsxfun(@times,1./FtPtScale,FtPt);
FtPt(isnan(FtPt)) = 0;

AScale = max(abs(A),[],1);
A = bsxfun(@times,1./AScale,A);
A(isnan(A))=0;

C = C*diag(AScale).*permute(FtPtScale,[3 1 2]);


%[Csorted,sortOrder] = sort(abs(C),2,'descend');

%TODO; permute to have same sign?

end
