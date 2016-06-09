classdef varDistributionC < handle
    properties (Access = private)
        % Data reference property
        data
    end
    
    properties
        ELBO
        
        % Probability Distributions
        pX
        pSigma
        pAlpha
        pP
        
        % Variational Factors
        qA
        qF
        qC
        qP
        qSigma
        qAlpha
        
        % Computed values
        XInnerProduct
        
        qPvonmisesEntropy
        
        
        % Settings
        method='manopt' % Method used to approximate E(qP)
        debugflag = 1
        activeParams = {'qA','qC','qF','qP','qSigma','qAlpha'}
    end
    
    properties
        % Terms in ELBO
        ePxz
        eQz
    end
    
    properties %(Access = private)
        % For updating moments, per k'th slab
        X
        eP
        eD
        eF
        qPmean
        
        eAsquared
        eFsquared
        eCsquared
        ePtP
        ePtPcond
        eFtPtPF
        eDFtPtPFD
        eAiDFtPtPFDAi
        eCtC
        eSigmaInv
        eAlphaInv
        
        % # ELBO terms
        % ## Means
        XqMeanLog
        AqMeanLog
        CqMeanLog
        FqMeanLog
        PqMeanLog
        SigmaqMeanLog
        AlphaqMeanLog
        % ## Entropies
        AqEntropy
        CqEntropy
        FqEntropy
        PqEntropy
        SigmaqEntropy
        AlphaqEntropy
    end
    
    methods
        function obj = varDistributionC(modelobj)
            % Class to controls all the variational factors and their
            % operations.
            % Constructed by an 'varBayesModelParafac2' object
            
            % Initialize Data
            obj.data = modelobj.data;
            
            obj.qA = multiNormalDist('qA',[obj.data.I obj.data.M]);
            obj.qC = multiNormalDist('qC',[obj.data.K obj.data.M]);
            obj.qF = multiNormalDist('qF',[obj.data.M obj.data.M]);
            obj.qP = multiNormalDist('qP',[obj.data.J obj.data.M obj.data.K]);
            obj.qSigma = inverseGammaDist('qSigma',[1 obj.data.K]);
            obj.qAlpha = inverseGammaDist('qAlpha',[1 obj.data.M]);
            
            
            obj.pSigma = inverseGammaDist('pSigma',[1 obj.data.K]);
            obj.pAlpha = inverseGammaDist('pAlpha',[1 obj.data.M]);
            obj.pP = multiNormalDist('pP',[obj.data.K obj.data.M]);
            
            
            if sum(sum(sum(obj.data.X)))==0
                obj.data.SigmaAtrue = 1;
                obj.data.SigmaBtrue = 1;
                obj.data.AlphaAtrue = 1;
                obj.data.AlphaBtrue = 1;
                
%                 obj.data.Sigmatrue = repmat(1e9,1,obj.data.K);
                               obj.data.Sigmatrue = 1./gamrnd(obj.data.SigmaAtrue,1/obj.data.SigmaBtrue,1,obj.data.K);
                obj.data.Alphatrue = 1./gamrnd(obj.data.AlphaAtrue,1/obj.data.AlphaBtrue,1,obj.data.Mtrue);
                
                obj.data.Atrue = 5*randn(obj.data.I,obj.data.Mtrue);
                obj.data.Ftrue = 5*randn(obj.data.Mtrue,obj.data.Mtrue);
                obj.data.Ptrue = 5*randn(obj.data.J,obj.data.Mtrue,obj.data.K);
                
                obj.data.Ctrue = 5*randn(obj.data.K,obj.data.Mtrue).*repmat(sqrt(obj.data.Alphatrue),obj.data.K,1);
                
                obj.data.Etrue = zeros(obj.data.I,obj.data.J,obj.data.K);
                obj.data.X = zeros(obj.data.I,obj.data.J,obj.data.K);
                
                for k = 1:obj.data.K
                    obj.data.Etrue(:,:,k) = randn(obj.data.I,obj.data.J)*sqrt(obj.data.Sigmatrue(k));
                    obj.data.X(:,:,k) = obj.data.Atrue*diag(obj.data.Ctrue(k,:))*obj.data.Ftrue'*obj.data.Ptrue(:,:,k)'+obj.data.Etrue(:,:,k);
                end
                
                
            end
            
            
            
            obj.X = obj.data.X;
            
            
            obj.XInnerProduct = obj.computeXInnerProduct;
            
        end
        
        function SNR(obj)
            
            disp(norm(obj.data.X(:))^2/norm(obj.data.Etrue(:))^2)
        end
        
        % #################################################################
        % # ELBO computations
        
        function value = get.ELBO(obj)
            % Compute the expected of 'log p(x,z)-log q(z)' w.r.t q(z)
            value = obj.ePxz+obj.eQz;
        end
        
        % ## Mean terms
        function value = get.ePxz(obj)
            % Recompute updated terms
            if any(ismember({'qA','qC','qF','qP','qSigma'},obj.activeParams)) || isempty(obj.XqMeanLog)
                obj.XqMeanLog = obj.computeXqMeanLog;
            end
            
            if ismember('qA',obj.activeParams) || isempty(obj.AqMeanLog)
                obj.AqMeanLog = obj.computeAqMeanLog;
            end
            
            if any(ismember({'qC','qAlpha'},obj.activeParams)) || isempty(obj.CqMeanLog)
                obj.CqMeanLog = obj.computeCqMeanLog;
            end
            
            if ismember('qF',obj.activeParams) || isempty(obj.FqMeanLog)
                obj.FqMeanLog = obj.computeFqMeanLog;
            end
            
            if ismember('qP',obj.activeParams) || isempty(obj.PqMeanLog)
                obj.PqMeanLog = obj.computePqMeanLog;
            end
            
            if ismember('qSigma',obj.activeParams) || isempty(obj.SigmaqMeanLog)
                obj.SigmaqMeanLog = obj.computeSigmaqMeanLog;
            end
            
            if ismember('qAlpha',obj.activeParams) || isempty(obj.AlphaqMeanLog)
                obj.AlphaqMeanLog = obj.computeAlphaqMeanLog;
            end
            
            % Compute expected log of P dists (first term ELBO)
            value = obj.XqMeanLog+obj.AqMeanLog+obj.CqMeanLog+...
                obj.FqMeanLog+obj.PqMeanLog+obj.SigmaqMeanLog+...
                obj.AlphaqMeanLog;
        end
        
        % ## Entropy terms
        function value = get.eQz(obj)
            % Recompute updated terms
            if ismember('qA',obj.activeParams) || isempty(obj.AqEntropy)
                obj.AqEntropy = obj.qA.entropy;
            end
            
            if ismember('qC',obj.activeParams) || isempty(obj.CqEntropy)
                obj.CqEntropy = obj.qC.entropy;
            end
            
            if ismember('qF',obj.activeParams) || isempty(obj.FqEntropy)
                obj.FqEntropy = obj.qF.entropy;
            end
            
            if ismember('qP',obj.activeParams) || isempty(obj.PqEntropy)
                if strcmp(obj.method,'vonmises')
                    if isempty(obj.qPvonmisesEntropy)
                        obj = obj.qPmean;
                    end
                    obj.PqEntropy = obj.qPvonmisesEntropy;
                else
                    obj.PqEntropy = obj.qP.entropy;
                end
            end
            
            if ismember('qSigma',obj.activeParams) || isempty(obj.SigmaqEntropy)
                obj.SigmaqEntropy = obj.qSigma.entropy;
            end
            
            if ismember('qAlpha',obj.activeParams) || isempty(obj.AlphaqEntropy)
                obj.AlphaqEntropy = obj.qAlpha.entropy;
            end
            
            
            % Compute sum of entropies
            value = obj.AqEntropy+obj.CqEntropy+obj.FqEntropy+...
                obj.PqEntropy+obj.SigmaqEntropy+obj.AlphaqEntropy;
        end
        
        % #################################################################
        % Mean values for ELBO
        function value = computeXqMeanLog(obj)
            t3 = zeros(1,obj.data.K);
            for k = 1:obj.data.K
                t3(k) = sum(sum(obj.data.X(:,:,k)*obj.eP(:,:,k)*obj.eF*obj.eD(:,:,k).*obj.qA.mean));
            end
            
            value = obj.data.J/2*sum(obj.qSigma.entropy)-1/2*sum(obj.qSigma.meanGamma.*(...
               sum(obj.eAiDFtPtPFDAi)+sum(obj.XInnerProduct)-2*t3));
        end
        
        function value = computeAqMeanLog(obj)
            value = -1/2*obj.qA.meanInnerProductSumComponent;
        end
        
        function value = computeCqMeanLog(obj)
            
            value = 1/2*obj.qAlpha.entropy-1/2*trace(obj.qC.mean*diag(obj.qAlpha.meanGamma)*obj.qC.mean');
        end
        
        function value = computeFqMeanLog(obj)
            value = -1/2*obj.qF.meanInnerProductSumComponent;
        end
        
        function value = computePqMeanLog(obj)
            value = -1/2*obj.qP.meanInnerProductSumComponent;
        end
        
        function value = computeSigmaqMeanLog(obj)
            value = sum((-obj.pSigma.alpha-1).*obj.qSigma.mean-obj.qSigma.meanGamma.*obj.pSigma.beta);
        end
        
        function value = computeAlphaqMeanLog(obj)
            value = sum((-obj.pAlpha.alpha-1).*obj.qAlpha.mean-obj.qAlpha.meanGamma.*obj.pAlpha.beta);
        end
        
        % #################################################################
        % # Moment Updates
        
        % ## Update control function
        function updateMoments(obj)
            if ismember('qA',obj.activeParams)
                obj.updateA;
%                 check_variance_matrix(obj.qA.variance);
            end
            if ismember('qC',obj.activeParams)
                obj.updateC;
%                 check_variance_matrix(obj.qC.variance);
            end
            if ismember('qF',obj.activeParams)
                obj.updateF;
%                 check_variance_matrix(obj.qF.variance);
            end
            if ismember('qP',obj.activeParams)
                obj.updateP;
%                 check_variance_matrix(obj.qP.variance);
            end
            if ismember('qSigma',obj.activeParams)
                obj.updateSigma;
            end
            if ismember('qAlpha',obj.activeParams)
                obj.updateAlpha;
            end
        end
        
        % ## Normal distributions
        % ### Variational Factor A
        function updateA(obj)
            
            ELBO_prev = obj.ELBO;
            for i = 1:obj.data.I
                obj.qA.variance(:,:,i) = inv(sum(multiplyTensor(obj.eDFtPtPFD,obj.qSigma.meanGamma),3)...
                    +eye(obj.data.M));
                
                sum_k=0;
                for k = 1:obj.data.K
                    sum_k = sum_k + obj.qSigma.meanGamma(k)*obj.eD(:,:,k)*obj.eF'*obj.eP(:,:,k)'*obj.data.X(i,:,k)';
                end
                
                check_ELBO(obj,ELBO_prev,obj.qA.varname,'variance',obj.debugflag)
                ELBO_prev = obj.ELBO;
                obj.qA.mean(i,:) = obj.qA.variance(:,:,i)*sum_k;
                check_ELBO(obj,ELBO_prev,obj.qA.varname,'mean',obj.debugflag)
            end
            
        end
        
        % ### Variational Factor C
        function updateC(obj)
            for k = 1:obj.data.K
                
                obj.qC.variance(:,:,k) = inv(obj.qSigma.meanGamma(k)*obj.qA.meanOuterProduct.*obj.eFtPtPF(:,:,k) + diag(obj.qAlpha.meanGamma));
                
                
                obj.qC.mean(k,:) = obj.qSigma.meanGamma(k)*diag(obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(:,:,k)'*obj.qA.mean)'*obj.qC.variance(:,:,k);
            end
        end
        
        % ### Variational Factor F
        function updateF(obj)
            for t = 1:1
            for m = 1:obj.data.M
                t1 = 0;
                for k = 1:obj.data.K
                    t1 = t1+obj.qSigma.meanGamma(k)*obj.eCtC(:,:,k).*obj.qA.meanOuterProduct*obj.ePtP(m,m,k);
                end
                obj.qF.variance(:,:,m) = inv(t1+eye(obj.data.M));
                
                
            end
            
            for m = 1:obj.data.M
                ELBO_prev = obj.ELBO;
                allButM = 1:obj.data.M~=m;
                tempMean = obj.qF.mean;
                t1 = zeros(1,obj.data.M);
                t2 = zeros(1,obj.data.M);
                for k = 1:obj.data.K
%                     if ~issymmetric(obj.eCtC(:,:,k).*obj.qA.meanOuterProduct)
%                         disp('Error')
%                         disp('Not symmetric') 
%                     end
                    t1 = t1+obj.qSigma.meanGamma(k)*(obj.eCtC(:,:,k).*obj.qA.meanOuterProduct*...
                        sum(repmat(obj.ePtP(m,allButM,k)',1,obj.data.M).*tempMean(allButM,:))')';
                    
                    t2 = t2+obj.qSigma.meanGamma(k)*obj.qP.mean(:,m,k)'*obj.data.X(:,:,k)'*obj.qA.mean*obj.eD(:,:,k);
                end
               check_ELBO(obj,ELBO_prev,obj.qF.varname,'variance',obj.debugflag) 
               ELBO_prev = obj.ELBO;
               obj.qF.mean(m,:) = (t2-t1)*obj.qF.variance(:,:,m);
               check_ELBO(obj,ELBO_prev,obj.qF.varname,'mean',obj.debugflag)
            end
            
            
%             obj.qF.mean=matricizing(X,1)*(krprod(C,B)*((C'*C).*(B'*B))^(-1))
            
            
            end
        end
        % ### Variational Factor P
        function updateP(obj)
            if ~strcmp(obj.method,'vonmises')
                for k = 1:obj.data.K
                    for j = 1:obj.data.J
                      obj.qP.variance(:,:,j,k) = inv(obj.qSigma.meanGamma(k)*(obj.qF.computeMeanInnerProductScaledSlabs(obj.eCtC(:,:,k).*obj.qA.meanOuterProduct))+eye(obj.data.M));
                    end
                end
            end
            obj.qPmean;
        end
        
        function obj = get.qPmean(obj)
            if strcmp(obj.method,'manopt')
               
                manifold = stiefelfactory(obj.data.J,obj.data.M);
                problem.M = manifold;
                problem.cost = @costFunc;
                problem.egrad = @gradFunc;
%                 k=1;
%                 checkgradient(problem)
                warning('off', 'manopt:getHessian:approx')
                options.verbosity=0;
                
                for k = 1:obj.data.K
                    
                    
                    % Constant terms in cost and grad function
                    gradconstant=(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k))';
                    costconstant=obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k);%*obj.qP.mean(:,:,k);
                    
                    cost = costFunc(obj.qP.mean(:,:,k));
                    
                    
                    [x,xcost] = trustregions(problem,[],options);
%                     fprintf('\n Slab %d with cost diff: %.2f \n',k,cost-xcost)
                    obj.qP.mean(:,:,k) = x;
                    
                end
                
            elseif strcmp(obj.method,'vonmises')
                obj.qPvonmisesEntropy = 0;
                for k=1:obj.data.K
                    D = obj.data.M;
                    X = obj.data.X(:,:,k);%obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)';
                    %[U,S,V]=svd(X,'econ');
                    %A = U(:,1:D)*S(1:D,1:D);
                    A = obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean'; %U(:,1:D)*S(1:D,1:D);
                    %dS=diag(S.^2);
                    %sigma_sq=(sum(dS)-sum(dS(1:3)))/numel(X);
                    p=obj.data.J;
                    
                    
                    % Estep - this step can be adapted directly to PARAFAC2 VB
                    F=A'*X/obj.qSigma.mean(k); %/sigma_sq;
                    [UU,SS,VV]=svd(F,'econ');
                    [f,~,lF]=hyperg(p,diag(SS),3);
                    %E_Z=UU*VV';                 % Expectation of Z
                    E_Z=UU*diag(f)*VV';                 % Expectation of Z
                    H_Z= lF-sum(sum(F.*E_Z));           % Entropy
                    %H_Z =0;
                    obj.qP.mean(:,:,k) = E_Z';
                    obj.qPvonmisesEntropy = obj.qPvonmisesEntropy+H_Z;
                end
                %
            elseif strcmp(obj.method,'parafac2svd')
                for k = 1:obj.data.K
                    [U,~,V] = svd(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k),'econ');
                    
                    obj.qP.mean(:,:,k) = V(:,1:obj.data.M)*U';
                end
            
            end
            
            
            function [cost] = costFunc(x)
                % if ~isfield(store,'cost')
                %                     disp(size(x))
                %                     size(obj.qP.mean(:,:,k))
%                 obj.qP.mean(:,:,k) = x;
%                 cost = -(obj.XqMeanLog+obj.PqMeanLog);
%                 cost = -(-1/2*obj.qSigma.mean(k)*(sum(obj.eAiDFtPtPFDAi(:,k))-...
%                     2*sum(sum(obj.data.X(:,:,k)*obj.qP.mean(:,:,k)*obj.qF.mean*obj.eD(:,:,k).*obj.qA.mean)))-1/2*...
%                     trace(obj.qP.meanInnerProductMatrix(:,:,k)));
             
%                 cost = -sum(sum(obj.data.X(:,:,k)*x*obj.qF.mean*obj.eD(:,:,k).*obj.qA.mean));
                cost = -trace(costconstant*x);
                %   cost = -(t3cost+t4cost*(sum(obj.eAiDFtPtPFDAi(:,k))+t5cost-2*sum(sum(t1cost*obj.eP(:,:,k)*t2cost)))...
            %        +obj.PqMeanLog);
                %end
                %cost = store.cost;
            end
            
            function [grad] = gradFunc(x)
                %if ~isfield(store,'cost')
                % [~,store] = costFunc(x);
                %end
                
                
                %if ~isfield(store,'grad')
                %obj.qP.mean(:,:,k) = x;
                
                %B = store.B;
                
%                 t1grad = 1/2*obj.qSigma.meanGamma(k);
%                 t2grad = 2*(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k))';
                %B = obj.qF.computeMeanInnerProductScaledSlabs(obj.eCtC(:,:,k).*obj.qA.meanOuterProduct);
                
%                 grad = t1grad*(x*B+x*B'+t2grad)+x;
%                 grad = -(-1/2*obj.qSigma.meanGamma(k)*(2*x*B...
%                     -2*(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k))')-x);
                grad = -gradconstant;
                %grad = -2*(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k))';
                %end
                %   grad = store.grad;
            end
        end
        
        % ## (Inverse-)Gamma distributions
        % ### Variational Factor Sigma
        function updateSigma(obj)
            for k = 1:obj.data.K
                ELBO_prev = obj.ELBO;
                obj.qSigma.alpha(k) = obj.pSigma.alpha(k)+obj.data.J*obj.data.I/2;
%                 check_ELBO(obj,ELBO_prev,obj.qSigma.varname,'alpha') 
                ELBO_prev = obj.ELBO;
                obj.qSigma.beta(k) = obj.pSigma.beta(k)+1/(1/2*sum(obj.eAiDFtPtPFDAi(:,k))...
                    +1/2*trace(obj.data.X(:,:,k)*obj.data.X(:,:,k)')...
                    -sum(sum(obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)'.*obj.data.X(:,:,k))));
%                 check_ELBO(obj,ELBO_prev,obj.qSigma.varname,'beta') 
            end
        end
        
        % ### Variational Factor Alpha
        function updateAlpha(obj)
            for m = 1:obj.data.M
                obj.qAlpha.alpha(m) = obj.pAlpha.alpha(m)+1/2*obj.data.K;
                obj.qAlpha.beta(m) = 1/(obj.pAlpha.beta(m)+1/2*sum(obj.eCsquared(:,m)));
            end
        end
        
        % #################################################################
        % # Terms for moment updates and means in the ELBO
        
        
        % ## First order
        function value = get.eP(obj)
            value = obj.qP.mean;
        end
        
        function value = get.eF(obj)
            value = obj.qF.mean;
        end
        
        function value = get.eD(obj)
            %value = cell(1,obj.data.K);
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            for k = 1:obj.data.K
                value(:,:,k) = diag(obj.qC.mean(k,:));
            end
        end
        
        % ## Second or Higher Order
        function value = get.eAsquared(obj)
            value = obj.qA.distSquared;
        end
        
        function value = get.eFsquared(obj)
            value = obj.qF.distSquared;
        end
        
        function value = get.eCsquared(obj)
            %             value = obj.qC.variance+obj.qC.mean.^2;
            value = obj.qC.distSquared;
        end
        
        function value = get.ePtP(obj)
            if strcmp(obj.method,'vonmises')
                value = repmat(eye(obj.data.M),1,1,obj.data.K);
            else
                value = squeeze(sum(obj.qP.variance,3))+repmat(obj.data.J*eye(obj.data.M),1,1,obj.data.K);%repmat(eye(obj.data.M),1,1,obj.data.K);
            end
        end
        
        function value = get.ePtPcond(obj)
            if strcmp(obj.method,'vonmises')
                value = repmat(eye(obj.data.M),1,1,obj.data.K);
            else
                
                varCond = zeros(obj.data.K,obj.data.M);
                
                for k=1:obj.data.K
                    for m = 1:obj.data.M
                        for i=1:obj.data.I
                            % Take the covariance matrix for i'th variable in k'th slab of P_i^k
                            covSigma = obj.qP.variance(:,:,i,k);
                            
                            % calculate conditional covariance for m'th element of P_i^k given the rest
                            varCond(k,m) = varCond(m) + covSigma(m,m) - covSigma(m,[1:(m-1) (m+1):end])*...
                                inv(covSigma([1:(m-1) (m+1):end],[1:(m-1) (m+1):end]))*covSigma([1:(m-1) (m+1):end],m);
                        end
                    end
                end
                value = zeros(obj.data.M,obj.data.M,obj.data.K);
                
                for k = 1:obj.data.K
                    value(:,:,k) = obj.qP.mean(:,:,k)'*obj.qP.mean(:,:,k)+diag(varCond(k,:));
                end
            end
        end
                
            function value = get.eFtPtPF(obj)
            % For all components
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            for k = 1:obj.data.K
                %cVar = sum(obj.qP.variance(:,:,:,k),3);
                cVar = obj.ePtP(:,:,k);
                for m = 1:obj.data.M
                    for n = 1:obj.data.M
                        if m == n
                            value(:,:,k) = value(:,:,k)+cVar(m,n)*(obj.qF.mean(m,:)'...
                                *obj.qF.mean(m,:)+obj.qF.variance(:,:,m));
                        else
                            value(:,:,k) = value(:,:,k)+cVar(m,n)*(obj.qF.mean(m,:)'...
                                *obj.qF.mean(n,:));
                        end
                    end
                end
            end
            % value = obj.qF.computeMeanInnerProductScaledSlabs(obj.ePtP);
        end
        
        function value = get.eCtC(obj)
            value = obj.qC.meanOuterProductSingle;
            %value = obj.qC.covMatrix;
        end
        
        function value = get.eDFtPtPFD(obj)
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            for k = 1:obj.data.K;
                %value = cellfun(@(a,b){a.*b},obj.eCtC,obj.eFtPtPF);
                value(:,:,k) = obj.eCtC(:,:,k).*obj.eFtPtPF(:,:,k);
            end
        end
        
        function value = get.eAiDFtPtPFDAi(obj)
            expected = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD);
            %expected = sum(expected,3);
            value = zeros(obj.data.I,obj.data.K);
            for k = 1:obj.data.K
                value(:,k) = diag(expected(:,:,k));
            end
        end
        
        function value = computeXInnerProduct(obj)
            value = zeros(obj.data.I,obj.data.K);
            for i = 1:obj.data.I
                for k = 1:obj.data.K
                    value(i,k) = obj.data.X(i,:,k)*obj.data.X(i,:,k)';
                end
            end
        end
        
        % #################################################################
        % # Display methods
        function displayPrivateDependent(obj)
            privateDependentProps = {'ePtP','eFtPtPF','eCtC','eSigmaInv'};
            
            for i = 1:numel(privateDependentProps)
                disp(obj.(privateDependentProps{i}))
            end
        end
        
    end
end

function check_variance_matrix(variance,debugflag)
if debugflag
I = size(variance,3);

if ndims(variance)>3
    K = size(variance,4);
else
    K = 1;
end

for k = 1:K
    for i = 1:I
        [~,p] = chol(variance(:,:,i,k));
        if p
            disp('error')
            keyboard
        end
    end
end
end
end

function check_ELBO(obj,ELBO_prev,var,updatedparam,debugflag)

if debugflag
diff = obj.ELBO-ELBO_prev;

if diff < -1e-6
    disp(diff)
    disp(var)
    disp(updatedparam)
%     keyboard
end
end
end





