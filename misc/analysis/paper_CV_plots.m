


%%
Ms = 2:10;
sum_test_ELBO = zeros(1,9);
for i = 1:9

[~,idx]=max(myModel{i}.CV_ELBOS_train,[],2);

% for k = 1:30
    sum_test_ELBO(i) = mean(max(myModel{i}.CV_ELBOS_test,[],2));
% end


end


plot(Ms,sum_test_ELBO/30,'.-')
grid on