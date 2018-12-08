function [P,F,C,A,fit] = Parafac2(X,M)
% Xk = Pk*F*Dk*A';
% D: M x M x K diagonal matrix
% F: M x M (Factor scores)
% A: J x M (loadings)
% P: I x M x K columnwise orthonormal
[I,J,K] = size(X);

% Initializing according to Kiers direct fitting paper
D = eye(M);
D = repmat(D,1,1,K);
F = eye(M);

XX = zeros(size(X,2));
for k = 1:K
    XX = XX + X(:,:,k)'*X(:,:,k);
end
A = pcacov(XX);
A = A(:,1:M);

%Allocating space
P = zeros(I,M,K);
C = zeros(K,M);
XP = zeros(M,J);

SSE = inf;
dSSE = inf;
SSE_old = SSE;
iter = 0;
maxIter = 2000;
while dSSE >= 1e-9*SSE_old && maxIter > iter;
    SSE_old=SSE;
    iter = iter+1;
    
    for k = 1:K
        [U,~,V]=svd(F*D(:,:,k)*A'*X(:,:,k)');
        P(:,:,k) = V(:,1:M)*U'; % Only use M first columns of V?
        XP(:,:,k) = P(:,:,k)'*X(:,:,k);
    end
    % C from the matrix notation of Parafac2 (n'th row is diag(D(:,:,n))
    C = undiagnolize(D);
    
    % Parafac1 Cycle for sum over k for norm2(Pk'*Xk-FDkA')^2 cost function
    [F,A,C] = parafac1(XP,F,A,C);
    
    D = diagnolize(C);
    
    SSE = 0;
    for k = 1:K
        SSE = SSE + sum(sum((X(:,:,k)-P(:,:,k)*F*D(:,:,k)*A').^2));
    end
    
    dSSE=SSE_old-SSE;
    display(SSE)
    %display(dSSE)
    
    SSE_old=SSE;
    
end

fit = SSE;
end

function B = undiagnolize(A)
B = zeros(size(A,3),size(A,1));
for i = 1:size(A,3)
    B(i,:) = diag(A(:,:,i))';
end
end

function A = diagnolize(B)
A = zeros(size(B,2),size(B,2),size(B,1));
for i = 1:size(A,3)
    A(:,:,i) = diag(B(i,:));
end
end

function [A,B,C] = parafac1(X,A,B,C)

A=matricizing(X,1)*(krprod(C,B)*((C'*C).*(B'*B))^(-1));
B=matricizing(X,2)*(krprod(C,A)*((C'*C).*(A'*A))^(-1));
C=matricizing(X,3)*(krprod(B,A)*((B'*B).*(A'*A))^(-1));


% ORDER ACCORDING TO VARIANCE
Tuck     = diag((A'*A).*(B'*B).*(C'*C));
[out,ID] = sort(Tuck);
A        = A(:,ID);
B        = B(:,ID);
C        = C(:,ID);
% NORMALIZE A AND C (variance in B)
Fac = size(A,2);
for f=1:Fac,normC(f) = norm(C(:,f));end
for f=1:Fac,normA(f) = norm(A(:,f));end
B        = B*diag(normC)*diag(normA);
A        = A*diag(normA.^(-1));
C        = C*diag(normC.^(-1));

% APPLY SIGN CONVENTION
SignA = sign(sum(sign(A))+eps);
SignC = sign(sum(sign(C))+eps);
A = A*diag(SignA);
C = C*diag(SignC);
B = B*diag(SignA)*diag(SignC);


end
