



for k = 1:1;   
subplot(2,3,[1 4])
imagesc(myModel.data.X(:,:,k))
colorbar
title(num2str(k))
subplot(2,3,[2 5])
imagesc(myModel.qDist.qA.mean*myModel.qDist.eD(:,:,k)*myModel.qDist.qF.mean'*myModel.qDist.qP.mean(:,:,k)')
colorbar
title(sprintf('variance %e',sqrt(myModel.qDist.qSigma.mean(k))))

subplot(2,3,3)
imagesc((myModel.data.X(:,:,k)-myModel.qDist.qA.mean*myModel.qDist.eD(:,:,k)*myModel.qDist.qF.mean'*myModel.qDist.qP.mean(:,:,k)')./myModel.data.X(:,:,k))
colorbar
title('diff')
subplot(2,3,6)
imagesc(myModel.data.Etrue(:,:,k))
colorbar
title('true error')

pause(1)
end