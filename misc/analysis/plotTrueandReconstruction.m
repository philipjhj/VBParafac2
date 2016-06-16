

figure(1)

for k = 1:myModel.data.K;   
% subplot(2,3,[1 4])
figure(1)
subplot(2,1,1)
imagesc(myModel.data.X(:,:,k))
colorbar
title(num2str(k))
% subplot(2,3,[2 5])
subplot(2,1,2)
imagesc(myModel.qDist.qA.mean*myModel.qDist.eD(:,:,k)*myModel.qDist.qF.mean'*myModel.qDist.qP.mean(:,:,k)')
colorbar
title(sprintf('variance %e',sqrt(myModel.qDist.qSigma.mean(k))))

% old 
% subplot(2,3,3)
% imagesc((myModel.data.X(:,:,k)-myModel.qDist.qA.mean*myModel.qDist.eD(:,:,k)*myModel.qDist.qF.mean'*myModel.qDist.qP.mean(:,:,k)')./myModel.data.X(:,:,k))
% colorbar
% title('diff')
% subplot(2,3,6)
% imagesc(myModel.data.Etrue(:,:,k))
% colorbar
% title('true error')

figure(2)
subplot(2,1,1)
imagesc(myModel.data.Ftrue'/max(max(myModel.data.Ftrue)))
title('F True')
colorbar
subplot(2,1,2)

imagesc(myModel.qDist.qF.mean'/max(max(myModel.qDist.qF.mean)))
title('qF mean')
colorbar

figure(3)
subplot(2,1,1)
imagesc(myModel.data.Atrue)
title('A True')
colorbar
subplot(2,1,2)
imagesc(myModel.qDist.qA.mean)
title('qA mean')
colorbar


figure(4)
subplot(2,1,1)
imagesc(myModel.data.Ctrue)
title('C True')
colorbar
subplot(2,1,2)
imagesc(myModel.qDist.qC.mean)
colorbar
title('qC mean')


figure(5)
subplot(2,1,1)
imagesc(myModel.data.Ptrue(:,:,k))
colorbar
title('Pk True')
subplot(2,1,2)
imagesc(myModel.qDist.qP.mean(:,:,k))
colorbar
title('qPk mean')


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









