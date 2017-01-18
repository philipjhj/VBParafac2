myModel.compute_reconstruction;

Xdiff=myModel.data.Xtrue-myModel.data.Xrecon;

norm(Xdiff(:),'fro')^2/norm(myModel.data.Xtrue(:),'fro')^2


% 0.3417
% 0.2308
% 0.3181
% 0.3758
% 0.4069

%%
normalModel.SSE_Fit(data.Xtrue)

% rng 4..8
% 3.8571
% 3.3431
% 3.7295
% 2.8050
% 3.2113