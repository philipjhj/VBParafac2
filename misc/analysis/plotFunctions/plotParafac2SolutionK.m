function plotParafac2SolutionK(k,X,A,C,F,P,Atrue,Ctrue,Ftrue,Ptrue,MLEflag)



if nargin < 11
    MLEflag = 0;
end

equalColorScale = 0;



if MLEflag
    nRows = 3;
else
    nRows = 2;
end


ha=tight_subplot(nRows,5,[0.1 0.02],[0.01 0.05],[0.015 0.015]);

xRecon = A*diag(C(k,:))*F'*P(:,:,k)';

% ### X plot
myImageTrue = X(:,:,k);
myImageEsti = xRecon;

if 0
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end

figure(1)
% subplot(2,5,1)
axes(ha(1))
titleText = strcat('Xk, k=',num2str(k));
displayImageValues(myImageTrue,titleText,colorInterval)

% subplot(2,5,6)
axes(ha(6))
titleText = strcat('Xk estimate, k =',num2str(k));
displayImageValues(myImageEsti,titleText,colorInterval)

% ### A Plot
myImageTrue = Atrue;
myImageEsti = A;

if equalColorScale
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end

% figure(3)
% subplot(2,5,3)
axes(ha(2))
titleText = 'A True';
displayImageValues(myImageTrue,titleText,colorInterval)

% subplot(2,5,8)
axes(ha(7))
titleText = 'qA mean';
displayImageValues(myImageEsti,titleText,colorInterval)



% ### C Plot
myImageTrue = diag(Ctrue(k,:));
myImageEsti = diag(C(k,:));

if equalColorScale
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end
% figure(4)
% subplot(2,5,4)
axes(ha(3))
titleText = 'diag(C_k True)';
displayImageValues(myImageTrue,titleText,colorInterval)

% subplot(2,5,9)
axes(ha(8))
titleText = 'diag(qC_k mean)';
displayImageValues(myImageEsti,titleText,colorInterval)


% ### F Plot
myImageTrue = Ftrue';
myImageEsti = F';

if equalColorScale
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end

% figure(2)
% subplot(2,5,2)
axes(ha(4))
titleText = 'F True';
displayImageValues(myImageTrue,titleText,colorInterval)

% subplot(2,5,7)
axes(ha(9))
titleText = 'qF mean';
displayImageValues(myImageEsti,titleText,colorInterval)

% ### P Plot
myImageTrue = Ptrue(:,:,k)';
myImageEsti = P(:,:,k)';

if equalColorScale
    colorInterval = findColorInterval(myImageTrue,myImageEsti);
else
    colorInterval = [];
end

% figure(5)
% subplot(2,5,5)
axes(ha(5))
titleText = 'Pk True';
displayImageValues(myImageTrue,titleText,colorInterval)

% subplot(2,5,10)
axes(ha(10))
titleText = 'qPk mean';
displayImageValues(myImageEsti,titleText,colorInterval)


if MLEflag
    
    
    [A,F,C,P,fit]=parafac2(X,size(A,2),[0 0],[0 0 0 0 1]);
    
    P = cat(3,P{:});
    
    XreconMLE = A*diag(C(k,:))*F'*P(:,:,k)';
    
    
    
    % ### X Plot
    myImageEsti = XreconMLE;
    
    axes(ha(11))
    titleText = 'Xk MLE';
    displayImageValues(myImageEsti,titleText,colorInterval)
    
    
    % ### A Plot
    myImageEsti = A;
    
    axes(ha(12))
    titleText = 'A MLE';
    displayImageValues(myImageEsti,titleText,colorInterval)
    
    
    % ### C Plot
    myImageEsti = diag(C(k,:));
    
    axes(ha(13))
    titleText = 'Ck MLE';
    displayImageValues(myImageEsti,titleText,colorInterval)
    
    % ### A Plot
    myImageEsti = F';
    
    axes(ha(14))
    titleText = 'F MLE';
    displayImageValues(myImageEsti,titleText,colorInterval)
    
    % ### A Plot
    myImageEsti = P(:,:,k)';
    
    axes(ha(15))
    titleText = 'Pk MLE';
    displayImageValues(myImageEsti,titleText,colorInterval)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end









end