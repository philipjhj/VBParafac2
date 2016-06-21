function [componentOrder] = sortData(A,C,F,P)







I = size(A,1);
M = size(A,2);
K = size(C,2);

B = zeros(I,M,K);
compVar = zeros(M,K);


for k = 1:K
   B(:,:,k)=A*diag(C(k,:))*F';  

   compVar(:,k)=var(B(:,:,k));
end


[v,componentOrder]=sort(mean(compVar,1),'descend');
% [v,componentOrder]=sort(var(A),'descend');


% disp(v)

end