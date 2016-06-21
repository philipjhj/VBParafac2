
figure(1)
clf

ha=tight_subplot(2,5,[0.1 0.02],[0.01 0.05],[0.015 0.015]);

equalColorScale = 0;

for k = 1:myModel.data.K;   
% subplot(2,3,[1 4])

xRecon = myModel.qDist.qA.mean*diag(myModel.qDist.qC.mean(k,:))*myModel.qDist.qF.mean'*myModel.qDist.qP.mean(:,:,k)';

% ### X plot
myImageTrue = myModel.data.X(:,:,k);
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
myImageTrue = myModel.data.Atrue;
myImageEsti = myModel.qDist.qA.mean;

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
myImageTrue = diag(myModel.data.Ctrue(k,:));
myImageEsti = diag(myModel.qDist.qC.mean(k,:));

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
myImageTrue = myModel.data.Ftrue';
myImageEsti = myModel.qDist.qF.mean';

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
myImageTrue = myModel.data.Ptrue(:,:,k)';
myImageEsti = myModel.qDist.qP.mean(:,:,k)';

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

pause

end


%%





for k = 1:myModel.data.K
    figure(4)
    subplot(1,2,1)
    imagesc(myModel.data.Ptrue(:,:,k))
    colorbar
    subplot(1,2,2)
    imagesc(myModel.qDist.qP.mean(:,:,k))
    colorbar
    pause(1)
end


%%


    figure(5)
    subplot(1,2,1)
    imagesc(myModel.data.Ctrue)
    colorbar
    subplot(1,2,2)
    imagesc(myModel.qDist.qC.mean)
    colorbar
    pause(3)
    
    %%
    
    myModel.qDist.qSigma
    myModel.qDist.qAlpha.mean









