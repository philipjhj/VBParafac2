function plotParafac2SolutionK(k,X,A,C,F,P,Atrue,Ctrue,Ftrue,Ptrue,MLEflag,Etrue)


M = size(A,2);
Mtrue = size(Atrue,2);

xRecon = A*diag(C(k,:))*F'*P(:,:,k)';

estiCell = {xRecon,A,C,F,P};
trueCell = {X,Atrue,Ctrue,Ftrue,Ptrue};


titleTrueCell ={'Xk, k=',
    'A True',
    'diag(C_k True)',
    'F True',
    'Pk True'};


titleEstiCell = {'Xk estimate, k =',
    'qA mean',
    'diag(qC_k mean)',
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


ha=tight_subplot(nRows,5,[0.1 0.02],[0.01 0.05],[0.015 0.015]);


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
displayImageValues(trueCell{1}(:,:,k),strcat(titleTrueCell{1},num2str(k)),colorInterval)

axes(ha(6))
displayImageValues(estiCell{1},strcat(titleEstiCell{1},num2str(k)),colorInterval)

if MLEflag
    axes(ha(11))
    displayImageValues(mleCell{1},titleMLECell{1},colorInterval)
end


% ### Others
for i = 2:numel(trueCell)
    
    if all(i ~= [3,5])
        myImageTrue = trueCell{i}(:,compOrderTrue);
        myImageEsti = estiCell{i}(:,compOrderEsti)*diag(signs(1,:));
        if MLEflag
            myImageMLE = mleCell{i}(:,compOrderMLE)*diag(signs(2,:));
        end
    elseif i == 3
        myImageTrue = diag(trueCell{i}(k,compOrderTrue));
        myImageEsti = diag(estiCell{i}(k,compOrderEsti))*diag(signs(1,:));
        if MLEflag
            myImageMLE = diag(mleCell{i}(k,compOrderMLE))*diag(signs(2,:));
        end
    elseif i == 5
        myImageTrue = trueCell{i}(:,compOrderTrue,k);
        myImageEsti = estiCell{i}(:,compOrderEsti,k)*diag(signs(1,:));
        if MLEflag
            myImageMLE = mleCell{i}(:,compOrderMLE,k)*diag(signs(2,:));
        end
    end
    
    axes(ha(i))
    displayImageValues(myImageTrue,titleTrueCell{i},colorInterval)
    
    axes(ha(i+5))
    displayImageValues(myImageEsti,titleEstiCell{i},colorInterval)
    
    
    if MLEflag
        axes(ha(i+10))
        displayImageValues(myImageMLE,titleMLECell{i},colorInterval)
    end
end


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
