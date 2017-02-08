function C=khatriRaoProduct(obj,A,B)
% Khatri Rao product
% input
% B m x n matrix
% C p x n matrix
% output
% A m*p x n matrix
sb=size(A,1);
sc=size(B,1);
C=zeros(sb*sc,size(A,2));
for k=1:size(A,2)
    C(:,k)=reshape(B(:,k)*A(:,k)',sb*sc,1);
end