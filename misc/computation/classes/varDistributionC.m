classdef varDistributionC < handle
    properties (Access = private)
        % Data reference property
        data
    end
    
    properties
        
        
        % Probability Distributions
        pX
        pSigma
        pAlpha
        
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
        activeParams_opt %= {'qP','qF','qC','qA','qSigma','qAlpha'}
    end
    
    properties (Dependent)
        % Terms in ELBO
        
        ELBO
        ePxz
        eQz
    end
    
    properties %\(Access = private)
        % For updating moments, per k'th slab
        X
        eD
        qPmean
        
        eCsquared
        ePtP
        eFtPtPF
        eDFtPtPFD
        eAiDFtPtPFDAi
        eCtC
        
        
        % # ELBO terms
        % ## Means
        qXMeanLog
        qAMeanLog
        qCMeanLog
        qFMeanLog
        qPMeanLog
        qSigmaMeanLog
        qAlphaMeanLog
        % ## Entropies
        qAEntropy
        qCEntropy
        qFEntropy
        qPEntropy
        qSigmaEntropy
        qAlphaEntropy
        
        activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'}
    end
    
    methods
        function obj = varDistributionC(modelobj)
            % Class to controls all the variational factors and their
            % operations.
            % Constructed by an 'varBayesModelParafac2' object
            
            % Initialize Data
            obj.data = modelobj.data;
            
            obj.qA = multiNormalDist('qA',[obj.data.I obj.data.M],true);
            obj.qC = multiNormalDist('qC',[obj.data.K obj.data.M]);
            obj.qF = multiNormalDist('qF',[obj.data.M obj.data.M]);
            obj.qP = multiNormalDist('qP',[obj.data.J obj.data.M obj.data.K],true);
            obj.qSigma = GammaDist('qSigma',[1 obj.data.K]);
            obj.qAlpha = GammaDist('qAlpha',[1 obj.data.M]);
            
            
            obj.pSigma = GammaDist('pSigma',[1 obj.data.K]);
            obj.pAlpha = GammaDist('pAlpha',[1 obj.data.M]);
            
            obj.pAlpha.alpha = obj.pAlpha.alpha;
            
            if sum(sum(sum(obj.data.X)))==0
                obj.data.SigmaAtrue = 1;
                obj.data.SigmaBtrue = 1;
                obj.data.AlphaAtrue = 1;
                obj.data.AlphaBtrue = 1;
                
                %                 obj.data.Sigmatrue = repmat(1e9,1,obj.data.K);
                obj.data.Sigmatrue = gamrnd(obj.data.SigmaAtrue,obj.data.SigmaBtrue,1,obj.data.K);
                obj.data.Alphatrue = gamrnd(obj.data.AlphaAtrue,obj.data.AlphaBtrue,1,obj.data.Mtrue);
                
                
                obj.data.Atrue = mvnrnd(zeros(obj.data.I,obj.data.Mtrue),eye(obj.data.Mtrue));
                obj.data.Ftrue = mvnrnd(zeros(obj.data.M,obj.data.Mtrue),eye(obj.data.Mtrue));
                
                obj.data.Ctrue = mvnrnd(zeros(obj.data.K,obj.data.Mtrue),1./obj.data.Alphatrue*eye(obj.data.Mtrue));
                
                obj.data.Etrue = zeros(obj.data.I,obj.data.J,obj.data.K);
                obj.data.X = zeros(obj.data.I,obj.data.J,obj.data.K);
                obj.data.Ptrue = zeros(obj.data.J,obj.data.M,obj.data.K);
                
                for k = 1:obj.data.K
                    
                    obj.data.Ptrue(:,:,k) = orth(mvnrnd(zeros(obj.data.J,obj.data.Mtrue),eye(obj.data.Mtrue)));
                
                    obj.data.Etrue(:,:,k) = mvnrnd(zeros(obj.data.I,obj.data.J)...
                        ,eye(obj.data.J)*1./obj.data.Sigmatrue(k));
                    obj.data.X(:,:,k) = obj.data.Atrue*diag(obj.data.Ctrue(k,:))*...
                        obj.data.Ftrue'*obj.data.Ptrue(:,:,k)'+obj.data.Etrue(:,:,k);
                end
                
                
            end
            obj.X = obj.data.X;
            obj.XInnerProduct = obj.computeXInnerProduct;
%             
            
            
            
            
            % Initialize Sufficient Stats
            obj.compute_eD;
            obj.compute_eCtC;
            obj.compute_eCsquared;
            obj.compute_ePtP;
            obj.compute_eFtPtPF;
            obj.compute_eDFtPtPFD;
            obj.compute_eAiDFtPtPFDAi;
            
            % Initialize ELBO terms
            all_params = {'qA','qC','qF','qP','qSigma','qAlpha'};
            
            obj.computeqXMeanLog
            for i  = 1:numel(all_params)
                obj.(all_params{i}).updateStatistics;
                
                methodStr = strcat('compute',all_params{i},'MeanLog');
%                 disp(methodStr)
                obj.(methodStr);
                
                methodStr = strcat(all_params{i},'Entropy');
                obj.(methodStr) = obj.(all_params{i}).entropy;
                
            end
            
%             obj.qA.mean = obj.data.Atrue;
%             obj.qC.mean = obj.data.Ctrue;
% %             obj.qF.mean = obj.data.Ftrue;
            obj.qP.mean = obj.data.Ptrue;
%             obj.qAlpha.mean = obj.data.Alphatrue;
%             obj.qSigma.mean = obj.data.Sigmatrue;
            
            
            
        end
        
        
        function set.activeParams_opt(obj,value)
            all_params = {'qP','qF','qC','qA','qSigma','qAlpha'};
            obj.activeParams_opt = value;
            obj.activeParams = all_params(ismember(all_params,obj.activeParams_opt));
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
            % Recompute for active parameters
            for i  = 1:numel(obj.activeParams)
                methodStr = strcat('compute',obj.activeParams{i},'MeanLog');
                obj.(methodStr);
            end
            
            if any(~ismember(obj.activeParams,'qAlpha'))
                obj.computeqXMeanLog % Depends on everything but qAlpha
            end
            
            % C depends on Alpha, so update even if C is not active
            if ismember('qAlpha',obj.activeParams) && ~ismember('qC',obj.activeParams)
                obj.computeqCMeanLog;
            end
            
            % Compute expected log of P dists (first term ELBO)
            value = obj.qXMeanLog+obj.qAMeanLog+obj.qCMeanLog+...
                obj.qFMeanLog+obj.qPMeanLog+obj.qSigmaMeanLog+...
                obj.qAlphaMeanLog;
        end
        
        % ## Entropy terms
        function value = get.eQz(obj)
            % Recompute updated terms
            
            for i = 1:numel(obj.activeParams)
                methodStr = strcat(obj.activeParams{i},'Entropy');
                obj.(methodStr) = obj.(obj.activeParams{i}).entropy;
            end
            
            % Compute sum of entropies
            value = obj.qAEntropy+obj.qCEntropy+obj.qFEntropy+...
                obj.qPEntropy+obj.qSigmaEntropy+obj.qAlphaEntropy;
        end
        
        % #################################################################
        % Mean values for ELBO
        function computeqXMeanLog(obj)
            t3 = zeros(obj.data.I,obj.data.M,obj.data.K);
            
            dataX = obj.data.X;
            qPMean = obj.qP.mean;
            qFMean = obj.qF.mean;
            qDMean = obj.eD;
            
            %             for k = 1:obj.data.K
            %                 t3(:,:,k) = dataX(:,:,k)*qPMean(:,:,k)*qFMean*qDMean(:,:,k);
            %             end
            
            t3 = mtimesx(mtimesx(mtimesx(dataX,qPMean),qFMean),qDMean);
            
            t3sum = sum(sum(sum(multiplyTensor(t3,obj.qSigma.mean),3).*obj.qA.mean));
            
            obj.qXMeanLog = obj.data.J/2*sum(obj.qSigma.entropy)-1/2*sum(obj.qSigma.mean.*(...
                sum(obj.eAiDFtPtPFDAi)+sum(obj.XInnerProduct)))+t3sum;
        end
        
        function computeqAMeanLog(obj)
            obj.qAMeanLog = -1/2*obj.qA.meanInnerProductSumComponent;
        end
        
        function computeqCMeanLog(obj)
            if isempty(obj.qAlpha.entropy)
                 obj.qAlpha.updateStatistics;
            end
            obj.qCMeanLog = 1/2*obj.qAlpha.entropy-1/2*trace(obj.qC.mean*diag(obj.qAlpha.mean)*obj.qC.mean');
        end
        
        function computeqFMeanLog(obj)
            obj.qFMeanLog = -1/2*obj.qF.meanInnerProductSumComponent;
        end
        
        function computeqPMeanLog(obj)
            obj.qPMeanLog = -1/2*obj.qP.meanInnerProductSumComponent;
        end
        
        function computeqSigmaMeanLog(obj)
            obj.qSigmaMeanLog = sum((-obj.pSigma.alpha-1).*obj.qSigma.mean-obj.qSigma.mean.*obj.pSigma.beta);
        end
        
        function computeqAlphaMeanLog(obj)
            obj.qAlphaMeanLog = sum((-obj.pAlpha.alpha-1).*obj.qAlpha.mean-obj.qAlpha.mean.*obj.pAlpha.beta);
        end
        
        % #################################################################
        % # Moment Updates
        
        % ## Update control function
        function updateMoments(obj)
            for i = 1:numel(obj.activeParams)
                methodStr = strcat('update',obj.activeParams{i});
                obj.(methodStr);
                obj.(obj.activeParams{i}).updateStatistics;
                
                if strcmp(obj.activeParams{i},'qP')
                    obj.compute_ePtP;
                elseif strcmp(obj.activeParams{i},'qF')
                    obj.compute_eFtPtPF;
                elseif strcmp(obj.activeParams{i},'qC')
                    obj.compute_eD;
                    obj.compute_eCtC;
                    obj.compute_eCsquared;
                    obj.compute_eDFtPtPFD;
                elseif strcmp(obj.activeParams{i},'qA')
                    obj.compute_eAiDFtPtPFDAi;
                end
            end
        end
        
        % ## Normal distributions
        % ### Variational Factor A
        function updateqA(obj)
            %
            %             ELBO_prev = obj.ELBO;
            obj.qA.variance = inv(sum(multiplyTensor(obj.eDFtPtPFD,obj.qSigma.mean),3)...
                +eye(obj.data.M));
            
            qSigmaMean = obj.qSigma.mean;
            qD = obj.eD;
            qFMeanT = obj.qF.mean';
            qPMeanT = permute(obj.qP.mean,[2 1 3]);
            dataXT = permute(obj.data.X,[2 1 3]);
            
            %             qAmean = obj.qA.mean;
            qAvariance = obj.qA.variance;
            
            K = obj.data.K;
            
            %             for i = 1:obj.data.I
            %
            %                 sum_k=0;
            %
            %                 for k = 1:K
            %                     sum_k = sum_k + qSigmaMean(k)*qD(:,:,k)*qFMeanT*qPMeanT(:,:,k)*dataXT(:,i,k);
            %                 end
            %
            %
            %
            %                 qAmean(i,:) = (qAvariance*sum_k)';
            %             end
            
            sum_k = sum(mtimesx(reshape(qSigmaMean,1,1,obj.data.K),mtimesx(mtimesx(mtimesx(qD,qFMeanT),qPMeanT),dataXT)),3);
            obj.qA.mean = (qAvariance*sum_k)';
            %             obj.qA.mean = qAmean;
            
        end
        
        % ### Variational Factor C
        function updateqC(obj)
            %             sum_k = 0;
            for k = 1:obj.data.K
                %                 ELBO_prev = obj.ELBO;
                obj.qC.variance(:,:,k) = inv(obj.qSigma.mean(k)*obj.qA.meanOuterProduct.*obj.eFtPtPF(:,:,k) + diag(obj.qAlpha.mean));
                %                 sum_k = sum_k+check_ELBO(obj,ELBO_prev,obj.qC.varname,'variance',obj.debugflag);
                %                 ELBO_prev = obj.ELBO;
                obj.qC.mean(k,:) = obj.qSigma.mean(k)*diag(obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(:,:,k)'*obj.qA.mean)'*obj.qC.variance(:,:,k);
                %                 sum_k = sum_k+check_ELBO(obj,ELBO_prev,obj.qC.varname,'mean',obj.debugflag);
            end
            %             disp(sum_k)
        end
        
        % ### Variational Factor F
        function updateqF(obj)
            %                 for m = 1:obj.data.M
            %                     t1 = 0;
            %                     for k = 1:obj.data.K
            %                         t1 = t1+obj.qSigma.mean(k)*obj.eCtC(:,:,k).*obj.qA.meanOuterProduct*obj.ePtP(m,m,k);
            %                     end
            %
            %
            %
            % %                     obj.qF.variance(:,:,m) = inv(t1+eye(obj.data.M));
            %                 end
            
            p = obj.data.M;
            n = obj.data.K;
            
            % Below code replaces a loop over m.
            % the bsxfun as index for ePtP gets the diagonals
            obj.qF.variance = multinv(squeeze(sum(bsxfun(@times,mtimesx(reshape(obj.qSigma.mean,1,1,obj.data.K),...
                bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct)),...
                reshape(obj.ePtP(bsxfun(@plus,[1:p+1:p*p]',[0:n-1]*p*p))',1,1,n,p)),3)));
            
            
            
            t2=squeeze(sum(mtimesx(reshape(obj.qSigma.mean,1,1,obj.data.K),...
                mtimesx(permute(obj.qP.mean,[4 1 3 2]),...
                mtimesx(permute(obj.data.X,[2 1 3]),...
                mtimesx(obj.qA.mean,obj.eD)))),3));
            
            
            for m = 1:obj.data.M
                allButM = 1:obj.data.M~=m;
                tempMean = obj.qF.mean;
                
                t1=sum(mtimesx(reshape(obj.qSigma.mean,1,1,obj.data.K),...
                    mtimesx(bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct),...
                    sum(bsxfun(@times,obj.ePtP(m,allButM,:),tempMean(allButM,:)'),2))),3)';
                %                             sum(repmat(,1,obj.data.M).*,1)')';
                
                %                                     t2=sum(mtimesx(reshape(obj.qSigma.mean,1,1,2),...
                %                                         mtimesx(permute(obj.qP.mean(:,m,:),[2 1 3]),...
                %                                         mtimesx(permute(obj.data.X,[2 1 3]),...
                %                                         mtimesx(obj.qA.mean,obj.eD)))),3);
                
                
                obj.qF.mean(m,:) = (t2(:,m)'-t1)*obj.qF.variance(:,:,m);
            end
            
            
            %                 for m = 1:obj.data.M
            %                     allButM = 1:obj.data.M~=m;
            %                     t1=squeeze(sum(mtimesx(reshape(obj.qSigma.mean,1,1,2),...
            %                         mtimesx(bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct),...
            %                         sum(bsxfun(@times,permute(obj.ePtP(:,allButM,:),[4 2 3 1]),tempMean(allButM,:)'),2))),3));
            %
            %                     t2=squeeze(sum(mtimesx(reshape(obj.qSigma.mean,1,1,2),...
            %                         mtimesx(permute(obj.qP.mean,[4 1 3 2]),...
            %                         mtimesx(permute(obj.data.X,[2 1 3]),...
            %                         mtimesx(obj.qA.mean,obj.eD)))),3));
            %
            %                     obj.qF.mean = squeeze(mtimesx(reshape(t2-t1,1,obj.data.M,obj.data.M),obj.qF.variance))';
            %                 end
            
        end
        % ### Variational Factor P
        function updateqP(obj)
            if ~strcmp(obj.method,'vonmises')
                for k = 1:obj.data.K
                    obj.qP.variance(:,:,1,k) = inv(obj.qSigma.mean(k)*...
                        (obj.qF.computeMeanInnerProductScaledSlabs(obj.eCtC(:,:,k).*...
                        obj.qA.meanOuterProduct))+eye(obj.data.M));
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
                    
                    %                     cost = costFunc(obj.qP.mean(:,:,k));
                    
                    
                    [x,~] = trustregions(problem,[],options);
                    %                     fprintf('\n Slab %d with cost diff: %.2f \n',k,cost-xcost)
                    obj.qP.mean(:,:,k) = x;
                    
                end
                
            elseif strcmp(obj.method,'vonmises')
                obj.qPvonmisesEntropy = 0;
                for k=1:obj.data.K
                    A = obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean'; %U(:,1:D)*S(1:D,1:D);
                    
                    % Estep - this step can be adapted directly to PARAFAC2 VB
                    F=A'*obj.data.X(:,:,k)/obj.qSigma.mean(k); %/sigma_sq;
                    [UU,SS,VV]=svd(F,'econ');
                    [f,~,lF]=hyperg(obj.data.J,diag(SS),3);
                    E_Z=UU*diag(f)*VV';                 % Expectation of Z
                    H_Z= lF-sum(sum(F.*E_Z));           % Entropy
                    obj.qP.mean(:,:,k) = E_Z';
                    obj.qPvonmisesEntropy = obj.qPvonmisesEntropy+H_Z;
                end
            elseif strcmp(obj.method,'parafac2svd')
                
                for k = 1:obj.data.K
                    [U,~,V] = svd(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k),'econ');
                    
                    obj.qP.mean(:,:,k) = V(:,1:obj.data.M)*U';
                end
                
            end
            
            function [cost] = costFunc(x)
                cost = -trace(costconstant*x);
            end
            
            function [grad] = gradFunc(x)
                grad = -gradconstant;
            end
        end
        
        % ## (Inverse-)Gamma distributions
        % ### Variational Factor Sigma
        function updateqSigma(obj)
            for k = 1:obj.data.K
                %ELBO_prev = obj.ELBO;
                obj.qSigma.alpha(k) = obj.pSigma.alpha(k)+obj.data.J*obj.data.I/2;
                %                 check_ELBO(obj,ELBO_prev,obj.qSigma.varname,'alpha')
                %ELBO_prev = obj.ELBO;
%                 obj.qSigma.beta(k) = obj.pSigma.beta(k)+1/(1/2*sum(obj.eAiDFtPtPFDAi(:,k))...
%                     +1/2*sum(obj.XInnerProduct)...
%                     -sum(sum(obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)'.*obj.data.X(:,:,k))));
                %                 check_ELBO(obj,ELBO_prev,obj.qSigma.varname,'beta')
            end
            obj.qSigma.beta = obj.pSigma.beta+1./(1/2*sum(obj.eAiDFtPtPFDAi,1)...
                    +1/2*sum(obj.XInnerProduct)...
                    -squeeze(sum(sum(...
                    bsxfun(@times,...
                    mtimesx(obj.qA.mean,mtimesx(obj.eD,mtimesx(obj.qF.mean',permute(obj.qP.mean,[2 1 3]))))...
                    ,obj.data.X),1),2))');
%                     -mtimesx(obj.qA.mean,mtimesx(obj.eD,mtimesx(obj.qF.mean',mtimesx(obj.qP.mean',obj.data.X')))));
        end
        
        % ### Variational Factor Alpha
        function updateqAlpha(obj)
            for m = 1:obj.data.M
                obj.qAlpha.alpha(m) = obj.pAlpha.alpha(m)+1/2*obj.data.K;
                obj.qAlpha.beta(m) = 1/(obj.pAlpha.beta(m)+1/2*sum(obj.eCsquared(:,m)));
            end
        end
        
        % #################################################################
        % # Terms for moment updates and means in the ELBO
        
        
        % ## First order
        
        function compute_eD(obj)
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            for k = 1:obj.data.K
                value(:,:,k) = diag(obj.qC.mean(k,:));
            end
            obj.eD = value;
        end
        
        % ## Second or Higher Order
        function compute_eCsquared(obj)
            obj.eCsquared = obj.qC.distSquared;
        end
        
        function compute_ePtP(obj)
            if strcmp(obj.method,'vonmises')
                obj.ePtP = repmat(eye(obj.data.M),1,1,obj.data.K);
            else
                obj.ePtP = obj.data.J*squeeze(obj.qP.variance)+repmat(eye(obj.data.M),1,1,obj.data.K);%repmat(eye(obj.data.M),1,1,obj.data.K);
            end
        end
        
        function compute_eFtPtPF(obj)
            % For all components
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            qFmeanT=obj.qF.mean';
            qFmean = obj.qF.mean;
            qFvariance = obj.qF.variance;
            
            ePtP = obj.ePtP; %#ok<PROP>
            
            M = obj.data.M;
            
            for k = 1:obj.data.K
                cVar = ePtP(:,:,k); %#ok<PROP>
                for m = 1:M
                    for n = 1:M
                        if cVar(m,n) ~= 0
                            if m == n
                                value(:,:,k) = value(:,:,k)+cVar(m,n)*(qFmeanT(:,m)...
                                    *qFmean(m,:)+qFvariance(:,:,m));
                            else
                                value(:,:,k) = value(:,:,k)+cVar(m,n)*(qFmeanT(:,m)...
                                    *qFmean(n,:));
                            end
                        end
                    end
                end
            end
            obj.eFtPtPF = value;
        end
        
        function compute_eCtC(obj)
            obj.eCtC = obj.qC.meanOuterProductSingle;
        end
        
        function compute_eDFtPtPFD(obj)
            value = zeros(obj.data.M,obj.data.M,obj.data.K);
            for k = 1:obj.data.K;
                value(:,:,k) = obj.eCtC(:,:,k).*obj.eFtPtPF(:,:,k);
            end
            obj.eDFtPtPFD = value;
        end
        
        function compute_eAiDFtPtPFDAi(obj)
            %             expected = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD);
            %             value = zeros(obj.data.I,obj.data.K);
            %             for k = 1:obj.data.K
            %                 value(:,k) = diag(expected(:,:,k));
            %             end
            obj.eAiDFtPtPFDAi = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD,1);
        end
        
        function value = computeXInnerProduct(obj)
            value = zeros(1,obj.data.K);
            %             for i = 1:obj.data.I
            for k = 1:obj.data.K
                value(k) = sum(sum(obj.data.X(:,:,k).^2));
            end
            %             end
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

% function check_variance_matrix(variance,debugflag)
% if debugflag
%     I = size(variance,3);
%
%     if ndims(variance)>3
%         K = size(variance,4);
%     else
%         K = 1;
%     end
%
%     for k = 1:K
%         for i = 1:I
%             [~,p] = chol(variance(:,:,i,k));
%             if p
%                 disp('error')
%                 keyboard
%             end
%         end
%     end
% end
% end

function val = check_ELBO(obj,ELBO_prev,var,updatedparam,debugflag)
val = 0;
if debugflag
    
    if strcmp(updatedparam,'qP')
        obj.compute_ePtP;
    elseif strcmp(updatedparam,'qF')
        obj.compute_eFtPtPF;
    elseif strcmp(updatedparam,'qC')
        obj.compute_eD;
        obj.compute_eCtC;
        obj.compute_eCsquared;
        obj.compute_eDFtPtPFD;
    elseif strcmp(updatedparam,'qA')
        obj.compute_eAiDFtPtPFDAi;
    end
    
    diff = obj.ELBO-ELBO_prev;
    
    if diff < -1e-6
        disp(diff)
        disp(var)
        disp(updatedparam)
        %     keyboard
        val = 1;
    end
    
end
end





