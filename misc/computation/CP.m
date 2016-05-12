function [A,B,C]=CP(X,D)

% The CandeComp/PARAFAC model optimized by alternating least squares for a
% 3-way array
%
% X_ijk=\sum_d A_id B_jd C_kd
%
% Usage:
% FACT=CP(X,D)
%
% Input:
% X             n-way array to decompose
% D             number of components
%
% Output:
% A             First mode loadings
% B             Second mode loadings
% C             Third mode loadings
optVerbose = 0;
N=size(X);

% Initialization
A=randn(N(1),D);
B=randn(N(2),D);
C=randn(N(3),D);

SSE=inf;
SST=sum(X(:).^2);
dSSE=inf;
maxiter=250;
tic;
iter=0;

if optVerbose
    disp([' '])
    disp(['3-way CP optimization'])
    disp(['A ' num2str(D) ' component model will be fitted']);
    dheader = sprintf('%12s | %12s | %12s | %12s ','Iteration','Expl. var.','dSSE','Time');
    dline = sprintf('-------------+--------------+--------------+--------------+');
end
while dSSE>=1e-9*SSE && iter<maxiter
    if mod(iter,100)==0 && optVerbose
        disp(dline); disp(dheader); disp(dline);
    end
    iter=iter+1;
    SSE_old=SSE;
    
    % Update factors making use of the functions
    %   matrizicing.m
    %   krprod.m
    A=matricizing(X,1)*(krprod(C,B)*((C'*C).*(B'*B))^(-1));
    B=matricizing(X,2)*(krprod(C,A)*((C'*C).*(A'*A))^(-1));
    C=matricizing(X,3)*(krprod(B,A)*((B'*B).*(A'*A))^(-1));
    
    % Sum of Square Error
    SSE = sum(sum(sum((matricizing(X,1)-A*(krprod(C,B)')).^2)));
    dSSE=SSE_old-SSE;
    
    if mod(iter,5)==0 && optVerbose
        disp(sprintf('%12.0f | %12.4f | %12.4f | %12.4e |',iter, (SST-SSE)/SST,dSSE,toc));
        tic;
    end
end
% Display final iteration
if optVerbose
    disp(sprintf('%12.0f | %12.4f | %12.4f | %12.4e |',iter, (SST-SSE)/SST,dSSE,toc));
end
