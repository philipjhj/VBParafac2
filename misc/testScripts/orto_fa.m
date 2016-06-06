function [A,S,costs] = orto_fa(X,d)
% Orthogonal factor analysis optimized with Manopt toolbox
% Factor A is constrained to be orthogonal and thus lies on the Stiefel
% manifold
% 
%   Problem: 
%       min L = || X - AS ||_F ^2
%           st: A'*A = I
%
%   Euclidean derivative of A:
%           dL/dA = dL/dA( tr(X'X - 2 X'AS + S'A'AS)
%                 = -2XS' + ASS' + ASS'
%                 = 2(AS -X)S'  
%   
% Written by: Søren Føns Vind Nielsen
% Copyright....none whatsoever - go to town
addpath(genpath('~/Documents/MATLAB/manopt'))

maxiter = 100;
cost_tol = 1e-6;

[p,n] = size(X);

A = randn(p,d);
S = randn(d,n);

% Define cost function, gradient and helper functions
% argument to be optimized should be first argument in both cost and
% gradient evaluation
    function [cost,store] = evalCost(A,store)
        R = X - A*S;
        cost = sum(sum(R.*R));
        store.cost = cost;
    end
    
    % Calculates gradient wrt A
    function [G,store] = evalGrad_A(A,store)
        if ~isfield(store,'G')
            store.G =  2*(A*S - X)*S';
        end
        G = store.G;
    end


% define struct for manopt
St = stiefelfactory(p, d);
problem.M = St;
problem.egrad = @evalGrad_A; 
problem.cost = @evalCost;
store = struct();

[cost,store] = evalCost(A,store);
[~,store] = evalGrad_A(A,store);

costs = nan(1,maxiter);
cost_diff = Inf;
i = 0;
while i<maxiter && abs(cost_diff)/abs(cost)>cost_tol
    % solve for A
    A = trustregions(problem);
    
    % compute S
    S = A'*X;
    
    % eval cost function
    i = i+1;
    costs(i) = evalCost(A,store);
    fprintf('L: %d , iter %d',costs(i),i)
    cost_diff = costs(i)-cost;
    cost = costs(i);
end



%eof
end





