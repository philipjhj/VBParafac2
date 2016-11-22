function plotParafac2SolutionK(k,X,A,C,F,P,Atrue,Ctrue,Ftrue,Ptrue,MLEflag,nActive,sortOrder)


M = size(A,2);
Mtrue = size(Atrue,2);

% K = size(X,3);




xRecon = A*diag(C(k,:))*F'*P(:,:,k)';
FtPt = F'*P(:,:,k)';

FtPtTrue = Ftrue'*Ptrue(:,:,k)';



% CtrueScale = max(abs(Ctrue)');
% Ctrue = bsxfun(@times,1./CtrueScale,Ctrue')';
% [~,sortOrderTrue] = sort(abs(Ctrue(k,:)),'descend');
% 


AtrueScale = max(abs(Atrue));
Atrue = bsxfun(@times,1./AtrueScale,Atrue);

FtPtTrueScale = max(abs(FtPtTrue),[],2);
FtPtTrue = bsxfun(@times,1./FtPtTrueScale,FtPtTrue);


Ctrue(k,:) = Ctrue(k,:).*AtrueScale*diag(FtPtTrueScale);

[~,sortOrderTrue] = sort(abs(Ctrue(k,:)),'descend');

% CScale = max(abs(C)');
% C = bsxfun(@times,1./CScale,C')';
% [~,sortOrder] = sort(abs(C(k,:)),'descend');

AScale = max(abs(A));
A = bsxfun(@times,1./AScale,A);

FtPtScale = max(abs(FtPt),[],2);
FtPtScale(FtPtScale<10e-16) = 0;

FtPt = bsxfun(@times,1./FtPtScale,FtPt);
FtPt(isinf(FtPt)) = 0;

C(k,:) = C(k,:).*AScale*diag(FtPtScale);

[~,sortOrder] = sort(abs(C(k,:)),'descend');

% sort out ratio
K = numel(sortOrderTrue);

% CRatioScale = 1./C(k,sortOrder(1:K)).*abs(Ctrue(k,sortOrderTrue));
% C(k,sortOrder(1:K)) = sign(C(k,sortOrder(1:K))).*abs(Ctrue(k,sortOrderTrue));
% A(:,sortOrder(1:K)) = bsxfun(@times,1./CRatioScale,A(:,sortOrder(1:K)));
% 


signATrue = sign(Atrue(1,:));
signA = sign(A(1,:));

signFtPtTrue = sign(FtPtTrue(:,1));
signFtPt = sign(FtPt(:,1));

signCTrue = sign(Ctrue(k,:));
signC = sign(C(k,:));

% signC(sortOrder)
% signCTrue(sortOrderTrue)
signA(sortOrder)
signATrue(sortOrderTrue)
% for i = 1:numel(signAFtPtTrue)
%     disp(signAFtPtTrue(i) ~= signAFtPt(i))
%     if signAFtPtTrue(i) ~= signAFtPt(i)
%        A(:,sortOrder(i)) = -1*A(:,sortOrder(i));
%        FtPt(sortOrder(i),:) = -1*FtPt(sortOrder(i),:);
%     end
% end


% for i = 1:numel(signCTrue)
% %     disp('sign')
% %     disp(signCTrue(sortOrderTrue(i)) ~= signC(sortOrder(i)))
%     disp('A sign')
%     disp(signATrue(sortOrderTrue(i)) ~= signA(sortOrder(i)))
%     if signCTrue(sortOrderTrue(i)) ~= signC(sortOrder(i))
%         
%         C(k,sortOrder(i)) = -1*C(k,sortOrder(i));
%         
%         if signATrue(sortOrderTrue(i)) ~= signA(sortOrder(i))
%             A(:,sortOrder(i)) = -1*A(:,sortOrder(i));
%         else
%             FtPt(sortOrder(i),:) = -1*FtPt(sortOrder(i),:);
%         end
%     elseif signATrue(sortOrderTrue(i)) ~= signA(sortOrder(i))
%             A(:,sortOrder(i)) = -1*A(:,sortOrder(i));
%             FtPt(sortOrder(i),:) = -1*FtPt(sortOrder(i),:);
%     end
% end

estiCell = {xRecon,A,C,FtPt,F,P};
trueCell = {X,Atrue,Ctrue,FtPtTrue,Ftrue,Ptrue};


titleTrueCell ={'Xk, k=',
    'A True',
    'diag(C_k True)',
    'FtPt True',
    'F True',
    'Pk True'};


titleEstiCell = {'Xk estimate, k =',
    'qA mean',
    'diag(qC_k mean)',
    'FtPt'
    'qF mean',
    'qPk mean'};

if MLEflag
    [A_MLE,F_MLE,C_MLE,P_MLE,~]=parafac2(X,M,[0 0],[0 0 0 0 1]);
    
    P_MLE = cat(3,P_MLE{:});
    
    XreconMLE = A_MLE*diag(C_MLE(k,:))*F_MLE'*P_MLE(:,:,k)';
    
    mleCell = {XreconMLE,A_MLE,C_MLE,F_MLE,P_MLE};
    
    %sortData(A_MLE,C_MLE,F_MLE,P_MLE);
    
    titleMLECell = {'Xk MLE',
        'A MLE',
        'diag(Ck) MLE',
        'F MLE',
        'P MLE'};
    
else
    mleCell = [];
end


if nargin < 11
    MLEflag = 0;
end

equalColorScale = 0;

if equalColorScale
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end

if MLEflag
    nRows = 3;
else
    nRows = 2;
end


nCols = 6;
ha=tight_subplot(nRows,nCols,[0.1 0.02],[0.01 0.05],[0.015 0.015]);


if M == Mtrue;
[componentGraph,signs] = sortData(trueCell,estiCell,mleCell);


    compOrderTrue=1:Mtrue;%sortData(Atrue,Ctrue,Ftrue,Ptrue);
    compOrderEsti=componentGraph(:,1);
    compOrderMLE=componentGraph(:,2);
else
    compOrderTrue=1:Mtrue;
    compOrderEsti=1:M;
    compOrderMLE=1:M;
    signs = ones(2,M);
end

% ### X plot
axes(ha(1))
displayImageValues(trueCell{1}(:,:,k),strcat(titleTrueCell{1},num2str(k)),colorInterval)%
colorbar
axes(ha(nCols+1))
displayImageValues(estiCell{1},strcat(titleEstiCell{1},num2str(k)),colorInterval)
colorbar
if MLEflag
    axes(ha(11))
    displayImageValues(mleCell{1},titleMLECell{1},colorInterval)
    colorbar
end


% ### Others
for i = 2:numel(trueCell)
    
    if all(i ~= [3,4,6])
        myImageTrue = trueCell{i}(:,sortOrderTrue);
        myImageEsti = estiCell{i}(:,sortOrder);
        if MLEflag
            myImageMLE = mleCell{i}(:,sortOrder);
        end
    elseif i == 3 % qC
        myImageTrue = diag(trueCell{i}(k,sortOrderTrue));
        myImageEsti = diag(estiCell{i}(k,sortOrder));
        if MLEflag
            myImageMLE = diag(mleCell{i}(k,:));
        end
    elseif i == 4 % FtPt
        myImageTrue = trueCell{i}(sortOrderTrue,:);
        myImageEsti = estiCell{i}(sortOrder,:);
    elseif i == 6
        myImageTrue = trueCell{i}(:,sortOrderTrue,k);
        myImageEsti = estiCell{i}(:,sortOrder,k);
        if MLEflag
            myImageMLE = mleCell{i}(:,:,k);
        end
    end
    
    axes(ha(i))
    displayImageValues(myImageTrue,titleTrueCell{i},colorInterval)
    
    if i == 3
        addValuesToImage(myImageTrue)
    end
    
    axes(ha(i+nCols))
    displayImageValues(myImageEsti,titleEstiCell{i},colorInterval)
    
    if i == 3
        addValuesToImage(myImageEsti)
    end
    
    if MLEflag
        axes(ha(i+10))
        displayImageValues(myImageMLE,titleMLECell{i},colorInterval)
        if i == 3
            addValuesToImage(myImageMLE)
        end
    end
end

% for i = 2:numel(trueCell)
%     
%     if all(i ~= [3,6])
%         myImageTrue = trueCell{i}(:,compOrderTrue);
%         myImageEsti = estiCell{i}(:,compOrderEsti)*diag(signs(1,:));
%         if MLEflag
%             myImageMLE = mleCell{i}(:,compOrderMLE)*diag(signs(2,:));
%         end
%     elseif i == 3
%         myImageTrue = diag(trueCell{i}(k,compOrderTrue));
%         myImageEsti = diag(estiCell{i}(k,compOrderEsti))*diag(signs(1,:));
%         if MLEflag
%             myImageMLE = diag(mleCell{i}(k,compOrderMLE))*diag(signs(2,:));
%         end
%     elseif i == 6
%         myImageTrue = trueCell{i}(:,compOrderTrue,k);
%         myImageEsti = estiCell{i}(:,compOrderEsti,k)*diag(signs(1,:));
%         if MLEflag
%             myImageMLE = mleCell{i}(:,compOrderMLE,k)*diag(signs(2,:));
%         end
%     end
%     
%     axes(ha(i))
%     displayImageValues(myImageTrue,titleTrueCell{i},colorInterval)
%     
%     if i == 3
%         addValuesToImage(myImageTrue)
%     end
%     
%     axes(ha(i+nCols))
%     displayImageValues(myImageEsti,titleEstiCell{i},colorInterval)
%     
%     if i == 3
%         addValuesToImage(myImageEsti)
%     end
%     
%     if MLEflag
%         axes(ha(i+10))
%         displayImageValues(myImageMLE,titleMLECell{i},colorInterval)
%         if i == 3
%             addValuesToImage(myImageMLE)
%         end
%     end
% end


end


function [componentGraph,signs] = sortData(trueCell,estiCell,mleCell)

Mtrue = size(trueCell{2},2);
M = size(estiCell{2},2);

componentGraph = repmat((1:M)',1,2);
signs = ones(2,M);

for m = 1:Mtrue
    trueSign = sign(trueCell{2}(:,m));%*diag(trueCell{3}(k,m))*trueCell{4}(:,m)');
    for n = 1:M
        estiSign = sign(estiCell{2}(:,n));%*diag(estiCell{3}(k,m))*estiCell{4}(:,m)');        
        if all(trueSign~=estiSign)
            componentGraph(m,1) = n;
            signs(1,m) = -1;
            break
        elseif all(trueSign==estiSign)
            componentGraph(m,1) = n;
            break
        end
    end
    for n = 1:M
        if ~isempty(mleCell)
            mleSign = sign(mleCell{2}(:,n));
            if all(trueSign~=mleSign)
                componentGraph(m,2) = n;
                signs(2,m) = -1;
                break
            elseif all(trueSign==mleSign)
                componentGraph(m,2) = n;
                break
            end
            
        end
    end
end

end
