load('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat')

%%
for k = 1:1%size(Y,3)
X = Y(:,:,k);
[~,~,~,LL]=EM_OPCA(X,2);
%%
plot(LL(2:end))
title(k)
pause(0.3)
% end