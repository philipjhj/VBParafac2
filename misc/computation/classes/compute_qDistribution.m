function obj = compute_qDistribution(obj)
% Implementation of CAVI to compute the variational distribution
% of the probalistic Parafac2 model
% Input:
%       obj;  class of type VBParafac2
%       Used properties:
%           X; Data matrix
%
% Output:
%       obj;  class of type VBParafac2
%       Used properties:
%           qDistribution; Computed variatonal distribution


qDistribution = obj.qDistribution;

% Initialize variational factors
% ...


% Compute ELBO
ELBO = compute_ELBO(qDistribution);
ELBO_prev = 0;

% Update Variational Factors until ELBO has converged
while ELBO-ELBO_prev > 1e-6
    
    % Update all variational factors except qP
    % ...
    
    % Approximate qP's parameters
    % ...
    
    % Compute ELBO
    ELBO_prev = ELBO;
    ELBO = compute_ELBO(qDistribution);
    
    % Output progress
    % ...
    
end

obj.qDistribution = qDistribution;

end

function ELBO = compute_ELBO(qDistribution)

% Compute the expected of log p(x,z) w.r.t q(z)
    E_pxz = compute_ElogpX(qDistribution)+compute_ElogpA(qDistribution)+...
        compute_ElogpC(qDistribution)+compute_ElogpF(qDistribution)+...
        compute_ElogpP(qDistribution)+compute_ElogpSigma(qDistribution)+...
        compute_ElogpAlpha(qDistribution);

% Compute the expected of log q(z) w.r.t q(z)
    E_qz = ...;

    ELBO = E_pxz-E_qz;

end

function E = compute_ElogpX(qDistribution)
    entropyIG =  @(a,b) a+log(b*gamma(a))-(1+a)*psi(a);
    meanIG = @(a,b) b/(a-1);
    meanXXt = @(Mu,Sigma){Sigma+Mu*Mu'};
    
%     entropySigma = qDistribution.J/2*cellfun(@(a,b){entropyIG(a,b)},...
%         qDistribution.qSigma(1,1,:),qDistribution.qSigma(1,2,:));
%     entropySigma = sum(cat(3,entropySigma{:}),3);
%     meanSigma = -1/2*cellfun(@(a,b){meanIG(a,b)},...
%         qDistribution.qSigma(1,1,:),qDistribution.qSigma(1,2,:));
%     meanPall = qDistribution.J*cellfun(@(Mu,Sigma){Sigma+Mu*Mu'},...
%         qDistribution.qP(1,1,:,:),qDistribution.qP(1,2,:,:));
%     meanPk = cell(1,qDistribution.K);
%     meanFk = cell(1,qDistribution.K);
%     for i = 1:qDistribution.K;
%         meanPk{i} = sum(cat(3,meanPall{i,:}),3);
%         meanFk = cellfun(@(Mu,Sigma,meanP){trace(meanP*Sigma)+Mu'*meanP*Mu},...
%             qDistribution.qF(1,1,:),qDistribution.qF(1,2,:));
%     end

    entropySigma=qDistribution.J/2*entropyIG(qDistribution.qSigma{1,1,:},...
        qDistribution.qSigma{1,2,:});
    meanSigma = -1/2*meanIG(qDistribution.qSigma{1,1,:},...
        qDistribution.qSigma{1,2,:});
    
    meanPall = qDistribution.J*meanXXt(qDistribution.qP{1,1,:,:},...
        qDistribution.qP{1,2,:,:});
    meanC =
    
    meanA=TR(
    

end

function E = compute_ElogpA(qDistribution)
    %moments=cellfun(@(Mu,Sigma){Sigma+Mu*Mu'},...
    %    qDistribution.qA(1,1,:),qDistribution.qA(1,2,:));
    %E = -1/2*sum(cat(3,moments{:}),3);
    
    E = -1/2*qDistribution.I*qDistribution.K; % chi-squared mean
end


function E = compute_ElogpC(qDistribution)
    entropyIG =  @(a,b) a+log(b*gamma(a))-(1+a)*psi(a);
    
    % Cell computations
%     E = cellfun(@(a,b){entropyIG(a,b)+qDistribution.M*(b^2*gamma(a+1/2)/gamma(a))},...
%         qDistribution.qAlpha(1,1,:),qDistribution.qAlpha(1,2,:));
%     E = qDistribution.K*sum(cat(3,E{:}),3);

    % Array computations
    ElogpC = @(a,b) qDistribution.K*sum(entropyIG(a,b)+qDistribution.M*...
        (b^2*gamma(a+1/2)/gamma(a)));
    
    E = ElogpC(qDistribution.qAlpha{1,1,:},qDistribution.qAlpha{1,2,:});
        
        
    
end


function E = compute_ElogpF(qDistribution)
    E = -1/2*qDistribution.M*qDistribution.M;
end

function E = compute_ElogpP(qDistribution)

end

