function [A,E_Z,sigma_sq,LL]=EM_OPCA(X,D)

% Model:
%  X = A*E_Z+E
%
% Where 
%    Z~UniformMatrixVonMisses
%    E_ij~Normal(X-A*E_Z,sigma_sq)
% 
% Z solved by VB-type of update, i.e. E-step
% A and sigma_sq solved by maximum likelihood, i.e. M-step
% 
% Input:
%  X  n x p data matrix
%  D  number of components
%
% Output:
%  A            n x D matrix
%  E_Z          D x p matrix
%  sigma_sq     noise variance
%  LL           EM-lower bound

[U,S,V]=svd(X,'econ');
A=U(:,1:D)*S(1:D,1:D);
dS=diag(S.^2);
sigma_sq=(sum(dS)-sum(dS(1:3)))/numel(X);
p=size(X,2);
maxiter=100;

const=0;
for j=1:D
    const=const+gammaln(0.5*(p-j+1));
end
logUniformPrior=-D*log(2)-0.5*p*D*log(pi)+1/4*D*(D-1)*log(pi)+const;

for k=1:maxiter
   
   % Estep - this step can be adapted directly to PARAFAC2 VB
   F=A'*X/sigma_sq;
   [UU,SS,VV]=svd(F,'econ');
   %diag(SS)
   [f,V,lF]=hyperg(p,diag(SS),5);
   E_Z=UU*diag(f)*VV';                 % Expectation of Z
   H_Z= lF-sum(sum(F.*E_Z));           % Entropy

   % Mstep
   A=(X*E_Z');                         % note since E_ZZ=I this term disappears
   sigma_sq=(sum(sum((X.^2)))-2*sum(sum(X.*(A*E_Z)))+sum(sum(A.*A)))/numel(X);
   
   % EM-Lower Bound
   LL(k)=-0.5*numel(X)*log(2*pi*sigma_sq)-0.5*numel(X)+logUniformPrior+H_Z; 
   %              logLikelihood                 + logPrior      + Entropy           
    
end