addpath ~/Documents/MATLAB/LatentVariableViz/

load motor_normalized_all_subs
%%
Z = reshape(Y,[size(Y,1),size(Y,2)*size(Y,3)]); % concatenating in time
tic
[U,S,V] = svd(Z,'econ'); % run svd
toc

%% plot first 10 PC's
plotComponent(U(:,1:10),V(:,1:10)', mask, [ 53    63    46])

