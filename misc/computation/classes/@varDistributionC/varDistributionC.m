classdef varDistributionC < handle
    properties (Access = private)
        % Data reference property
        data
        
    end
    
    properties
        
        
        % Probability Distributions
        pSigma
        pAlpha
        
        % Variational Factors
        qA
        qF
        qC
        qP
        qSigma
        qAlpha
        
        % Settings
        methodPesti='parafac2svd' % Method used to approximate E(qP)
        methodMatrixProduct = 'mtimesx'
        debugflag = 1
        activeParams_opt = {'qA','qC','qF','qP','qSigma','qAlpha'}
        
    end
    
    properties (Dependent)
        % ELBO Terms
        ELBO
        ePxz
        eQz
    end
    
    properties %(Access = private)
        
        % Computed values
        XInnerProduct
        qPvonmisesEntropy
        
        % Shared terms between moments
        eD
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
    end
    properties %(Access = protected)
        activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'}
    end
    
    methods
        function obj = varDistributionC(modelobj)
            % Class to controls all the variational factors and their
            % operations.
            % Constructed by an 'varBayesModelParafac2' object
            
            % Initialize Data
            obj.data = modelobj.data;
            
            % Initialize distribution classes
            obj.qA = multiNormalDist('qA',[obj.data.I obj.data.M],true);
            obj.qC = multiNormalDist('qC',[obj.data.K obj.data.M]);
            obj.qF = multiNormalDist('qF',[obj.data.M obj.data.M]);
            obj.qP = multiNormalDist('qP',[obj.data.J obj.data.M obj.data.K],true);
            obj.qSigma = GammaDist('qSigma',[1 obj.data.K]);
            obj.qAlpha = GammaDist('qAlpha',[1 obj.data.M]);
            
            obj.pSigma = GammaDist('pSigma',[1 1]);
            obj.pAlpha = GammaDist('pAlpha',[1 1]);
            
            
            
            if strcmp(obj.methodMatrixProduct,'gpu')
                all_params = {'qA','qC','qF','qP'};
                for i  = 1:numel(all_params)
                    obj.(all_params{i}).mean = gpuArray(obj.(all_params{i}).mean);
                    obj.(all_params{i}).variance = gpuArray(obj.(all_params{i}).variance);
                end
            end
            
            
            obj.XInnerProduct = obj.computeXInnerProduct;
            
            % Use same initialization as the original parafac2 code
            %             [A,F,C,P]=obj.parafac2([0 0],[0, -1, 0,0,1]);
            %
            %             obj.qA.mean = A;
            %             obj.qC.mean = C;
            %             obj.qF.mean = F;
            %             obj.qP.mean = cat(3,P{:});
            
            % Initialize Shared Terms
            obj.compute_eD;
            obj.compute_eCtC;
            obj.compute_eCsquared;
            obj.compute_ePtP;
            obj.compute_eFtPtPF;
            obj.compute_eDFtPtPFD;
            obj.compute_eAiDFtPtPFDAi;
            
            % Initialize ELBO terms
            
            
            all_params = {'qA','qC','qF','qP','qSigma','qAlpha'};
            for i  = 1:numel(all_params)
                obj.(all_params{i}).updateStatistics;
                
                methodStr = strcat('compute',all_params{i},'MeanLog');
                obj.(methodStr);
                
                methodStr = strcat(all_params{i},'Entropy');
                obj.(methodStr) = obj.(all_params{i}).entropy;
                
            end
            
            obj.computeqXMeanLog
            
            % Uncomment to set values to the ground truth
            %                         obj.qA.mean = [obj.data.Atrue zeros(obj.data.I,obj.data.M-obj.data.Mtrue)];
            %                         obj.qC.mean = [obj.data.Ctrue zeros(obj.data.K,obj.data.M-obj.data.Mtrue)];
%                                     obj.qF.mean = obj.data.Ftrue;
            %                         obj.qP.mean = obj.data.Ptrue;
%                                     obj.qAlpha.mean = [obj.data.Alphatrue 1e9*ones(1,obj.data.M-obj.data.Mtrue)];
%                                     obj.qSigma.mean = obj.data.Sigmatrue;
            
            
        end
        
        
        function set.activeParams_opt(obj,value)
            all_params = {'qP','qF','qC','qA','qSigma','qAlpha'};
            obj.activeParams_opt = value;
            obj.activeParams = all_params(ismember(all_params,obj.activeParams_opt));
        end
        
        
        function SNR(obj)
            disp(norm(obj.data.X(:))^2/norm(obj.data.Etrue(:))^2)
        end
        
        % Function by Rasmus Bro to get initial values
        [A,H,C,P,fit,AddiOutput] = parafac2(obj,Constraints,Options,A,H,C,P);
        
        
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
            
            t3 = matrixProductPrSlab(obj,matrixProductPrSlab(obj,matrixProductPrSlab(obj,dataX,qPMean),qFMean),qDMean);
            
            t3sum = sum(sum(sum(multiplyTensor(t3,obj.qSigma.mean),3).*obj.qA.mean));
            
            obj.qXMeanLog = obj.data.J*obj.data.I/2*sum(obj.qSigma.MeanLog)-1/2*sum(obj.qSigma.mean.*(...
                sum(obj.eAiDFtPtPFDAi)+obj.XInnerProduct))+t3sum;
        end
        
        function computeqAMeanLog(obj)
            obj.qAMeanLog = -1/2*obj.qA.meanInnerProductSumComponent;
        end
        
        function computeqCMeanLog(obj)
            if isempty(obj.qAlpha.entropy)
                obj.qAlpha.updateStatistics;
            end
            obj.qCMeanLog = 1/2*sum(obj.qAlpha.MeanLog)-1/2*trace(obj.qC.mean*diag(obj.qAlpha.mean)*obj.qC.mean');
        end
        
        function computeqFMeanLog(obj)
            obj.qFMeanLog = -1/2*obj.qF.meanInnerProductSumComponent;
        end
        
        function computeqPMeanLog(obj)
            obj.qPMeanLog = -1/2*obj.qP.meanInnerProductSumComponent;
        end
        
        function computeqSigmaMeanLog(obj)
            obj.qSigmaMeanLog = sum((obj.pSigma.alpha-1).*obj.qSigma.MeanLog-obj.qSigma.mean.*1/obj.pSigma.beta);
        end
        
        function computeqAlphaMeanLog(obj)
            obj.qAlphaMeanLog = sum((obj.pAlpha.alpha-1)*obj.qAlpha.MeanLog-obj.qAlpha.mean.*1/obj.pAlpha.beta);
        end
        
        % #################################################################
        % # Moment Updates
        
        % ## Update control function
        function updateMoments(obj)
            any_error = 0;
            first_error = 0;
            
            for i = 1:numel(obj.activeParams)
                %                 ELBO_prev = obj.ELBO;
                
                if any(~ismember(obj.activeParams,'qAlpha'))
                    XMeanLog_prev = obj.qXMeanLog; % Depends on everything but qAlpha
                else
                    XMeanLog_prev= 0;
                end
                
                if strcmp(obj.activeParams{i},'qAlpha')
                    CMeanLog_prev = obj.qCMeanLog; % Depends on everything but qAlpha
                else
                    CMeanLog_prev= 0;
                end
                
                methodStr = strcat(obj.activeParams{i},'MeanLog');
                
                ELBO_prev = XMeanLog_prev+CMeanLog_prev+...
                    obj.(methodStr)+obj.(obj.activeParams{i}).entropy;
                
                ELBO_prev2 = obj.ELBO;
                
                if ~ismember(obj.activeParams{i},{'qAlpha','qSigma'}) || obj.data.iter>=25% || obj.data.iter==0
                    
                    methodStr = strcat('update',obj.activeParams{i});
                    obj.(methodStr);
                    obj.(obj.activeParams{i}).updateStatistics;
                end
                
                
                
                % Update shared terms
                if strcmp(obj.activeParams{i},'qP')
                    obj.compute_ePtP;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qF',obj.activeParams)
                        obj.compute_eFtPtPF;
                    elseif all(~ismember({'qC','qF'},obj.activeParams))
                        obj.compute_eDFtPtPFD;
                    elseif all(~ismember({'qC','qF','qA'},obj.activeParams))
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.activeParams{i},'qF')
                    obj.compute_eFtPtPF;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qC',obj.activeParams)
                        obj.compute_eDFtPtPFD;
                    elseif all(~ismember({'qC','qA'},obj.activeParams))
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.activeParams{i},'qC')
                    obj.compute_eD;
                    obj.compute_eCtC;
                    obj.compute_eCsquared;
                    obj.compute_eDFtPtPFD;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qA',obj.activeParams)
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.activeParams{i},'qA')
                    obj.compute_eAiDFtPtPFDAi;
                end
                
                
                if any(~ismember(obj.activeParams,'qAlpha'))
                    obj.computeqXMeanLog;
                    XMeanLog_new = obj.qXMeanLog; % Depends on everything but qAlpha
                else
                    XMeanLog_new=0;
                end
                
                if strcmp(obj.activeParams{i},'qAlpha')
                    obj.computeqCMeanLog;
                    CMeanLog_new = obj.qCMeanLog; % Depends on everything but qAlpha
                else
                    CMeanLog_new= 0;
                end
                
                
                methodStr = strcat('compute',obj.activeParams{i},'MeanLog');
                obj.(methodStr);
                
                methodStr = strcat(obj.activeParams{i},'MeanLog');
                
                ELBO_new = XMeanLog_new+CMeanLog_new+...
                    obj.(methodStr)+obj.(obj.activeParams{i}).entropy;
                
                ELBO_new2 = obj.ELBO;
                
                if obj.debugflag && (ELBO_new-ELBO_prev)/abs(ELBO_new) < -1e-8 && obj.data.iter > 0
                    any_error = 1;
                    if ~first_error
                        first_error = 1;
                        fprintf('\n');
                    end
                    warning('off','backtrace')
                    warning(sprintf('varDist:Update:%s',obj.activeParams{i}),...
                        'Update problem; ELBO absolute/relative change for %s update is %f / %f \t %f\n',obj.activeParams{i},...
                        ELBO_new-ELBO_prev,(ELBO_new-ELBO_prev)/abs(ELBO_new),(ELBO_new2-ELBO_prev2)/abs(ELBO_new2));
                    warning('on','backtrace')
                    
                end
                
                obj.data.ELBOall = [obj.data.ELBOall obj.ELBO];
            end
            if any_error
                fprintf('\n')
            end
        end
        
        % ## Normal distributions
        % ### Variational Factor A
        function updateqA(obj)
            
            
%             obj.qA.variance = inv(sum(multiplyTensor(obj.eDFtPtPFD,obj.qSigma.mean),3)...
%                     +eye(obj.data.M));
                
            
%             for i = 1:obj.data.I
                
%                 sum_k=0;
%                 for k = 1:obj.data.K
%                     sum_k = sum_k + obj.qSigma.mean(k)*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(i,:,k)';
%                 end
%                 
% %                 check_ELBO(obj,ELBO_prev,obj.qA.varname,'variance',obj.debugflag)
% %                 ELBO_prev = obj.ELBO;
%                 obj.qA.mean(i,:) = obj.qA.variance*sum_k;
% %                 check_ELBO(obj,ELBO_prev,obj.qA.varname,'mean',obj.debugflag)
%             end

            
            %
            %             ELBO_prev = obj.ELBO;
            obj.qA.variance = inv(sum(multiplyTensor(obj.eDFtPtPFD,obj.qSigma.mean),3)...
                +eye(obj.data.M));
            
            qSigmaMean = obj.qSigma.mean;
            qD = obj.eD;
            qFMeanT = obj.qF.mean';
            qPMeanT = permute(obj.qP.mean,[2 1 3]);
            dataXT = permute(obj.data.X,[2 1 3]);
            
            qAvariance = obj.qA.variance;
            
            K = obj.data.K;
            sum_k = sum(bsxfun(@times,reshape(qSigmaMean,1,1,obj.data.K),...
                matrixProductPrSlab(obj,matrixProductPrSlab(obj,...
                matrixProductPrSlab(obj,qD,qFMeanT),qPMeanT),dataXT)),3);
            obj.qA.mean = (qAvariance*sum_k)';
%             
        end
        
        % ### Variational Factor C
        function updateqC(obj)
            for k = 1:obj.data.K
                
                obj.qC.variance(:,:,k) = inv(obj.qSigma.mean(k)*obj.qA.meanOuterProduct.*obj.eFtPtPF(:,:,k) + diag(obj.qAlpha.mean));
                
                obj.qC.mean(k,:) = obj.qSigma.mean(k)*diag(obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(:,:,k)'*obj.qA.mean)'*obj.qC.variance(:,:,k);
            end
        end
        
        % ### Variational Factor F
        function updateqF(obj)
            
            p = obj.data.M;
            n = obj.data.K;
            
            % Below code replaces a loop over m.
            % the bsxfun as index for ePtP gets the diagonals
            
            obj.qF.variance = multinv(bsxfun(@plus,squeeze(sum(bsxfun(@times,bsxfun(@times,...
                reshape(obj.qSigma.mean,1,1,obj.data.K),...
                bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct)),...
                reshape(obj.ePtP(bsxfun(@plus,[1:p+1:p*p]',[0:n-1]*p*p))',1,1,n,p)),3)),eye(obj.data.M)));
            
            
            
            t2=squeeze(sum(bsxfun(@times,reshape(obj.qSigma.mean,1,1,obj.data.K),...
                matrixProductPrSlab(obj,permute(obj.qP.mean,[2 1 3]),...
                matrixProductPrSlab(obj,permute(obj.data.X,[2 1 3]),...
                matrixProductPrSlab(obj,obj.qA.mean,obj.eD)))),3));
            
            
            for m = 1:obj.data.M
                allButM = 1:obj.data.M~=m;
                tempMean = obj.qF.mean;
                
                t1=sum(bsxfun(@times,reshape(obj.qSigma.mean,1,1,obj.data.K),...
                    matrixProductPrSlab(obj,bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct),...
                    sum(bsxfun(@times,obj.ePtP(m,allButM,:),tempMean(allButM,:)'),2))),3)';
                
                obj.qF.mean(m,:) = (t2(m,:)-t1)*obj.qF.variance(:,:,m);
            end
            
            
            
        end
        % ### Variational Factor P
        function updateqP(obj)
            if ~strcmp(obj.methodPesti,'vonmises')
                for k = 1:obj.data.K
                    obj.qP.variance(:,:,1,k) = inv(obj.qSigma.mean(k)*...
                        (obj.qF.computeMeanInnerProductScaledSlabs(obj.eCtC(:,:,k).*...
                        obj.qA.meanOuterProduct))+eye(obj.data.M));
                end
            end
            obj.computeqPmean;
        end
        
        function computeqPmean(obj)
            if strcmp(obj.methodPesti,'manopt')
                
                manifold = stiefelfactory(obj.data.J,obj.data.M);
                problem.M = manifold;
                problem.cost = @costFunc;
                problem.egrad = @gradFunc;
                %                 checkgradient(problem)
                warning('off', 'manopt:getHessian:approx')
                options.verbosity=0;
                
                for k = 1:obj.data.K
                    
                    
                    % Constant terms in cost and grad function
                    gradconstant=(obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k))';
                    costconstant=obj.qF.mean*obj.eD(:,:,k)*obj.qA.mean'*obj.data.X(:,:,k);%*obj.qP.mean(:,:,k);
                    
                    % cost = costFunc(obj.qP.mean(:,:,k));
                    
                    
                    [x,~] = trustregions(problem,[],options);
                    % fprintf('\n Slab %d with cost diff: %.2f \n',k,cost-xcost)
                    obj.qP.mean(:,:,k) = x;
                    
                end
                
            elseif strcmp(obj.methodPesti,'vonmises')
                obj.qPvonmisesEntropy = 0;
                for k=1:obj.data.K
                    A = obj.qA.mean*obj.eD(:,:,k)*obj.qF.mean';
                    
                    % Estep
                    F=A'*obj.data.X(:,:,k)/obj.qSigma.mean(k); %/sigma_sq;
                    [UU,SS,VV]=svd(F,'econ');
                    [f,~,lF]=hyperg(obj.data.J,diag(SS),3);
                    E_Z=UU*diag(f)*VV';                 % Expectation of Z
                    H_Z= lF-sum(sum(F.*E_Z));           % Entropy
                    obj.qP.mean(:,:,k) = E_Z';
                    obj.qPvonmisesEntropy = obj.qPvonmisesEntropy+H_Z;
                end
            elseif strcmp(obj.methodPesti,'parafac2svd')
                
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
                obj.qSigma.alpha(k) = obj.pSigma.alpha+obj.data.I*obj.data.J/2;
            end
            obj.qSigma.beta = 1./(1./obj.pSigma.beta+1/2*sum(obj.eAiDFtPtPFDAi,1)...
                +1/2*obj.XInnerProduct...
                -squeeze(sum(sum(...
                bsxfun(@times,...
                matrixProductPrSlab(obj,obj.qA.mean,...
                matrixProductPrSlab(obj,obj.eD,...
                matrixProductPrSlab(obj,obj.qF.mean',...
                permute(obj.qP.mean,[2 1 3]))))...
                ,obj.data.X),1),2))');
        end
        
        
        
        % ### Variational Factor Alpha
        function updateqAlpha(obj)
            for m = 1:obj.data.M
                obj.qAlpha.alpha(m) = obj.pAlpha.alpha+1/2*obj.data.K;
                obj.qAlpha.beta(m) = 1/(1/obj.pAlpha.beta+1/2*sum(obj.eCsquared(:,m)));
            end
        end
        
        % #################################################################
        % # Terms for moment updates and means in the ELBO
        
        
        % ## First order
        
        function compute_eD(obj)
%             value = zeros(obj.data.M,obj.data.M,obj.data.K);
%             for k = 1:obj.data.K
%                 value(:,:,k) = diag(obj.qC.mean(k,:));
%             end
%             obj.eD = value;
            obj.eD = bsxfun(@mtimes,reshape(obj.qC.mean',1,...
                obj.data.M,obj.data.K),...
                repmat(eye(obj.data.M),1,1,obj.data.K));
        end
        
        % ## Second or Higher Order
        function compute_eCsquared(obj)
            obj.eCsquared = obj.qC.distSquared;
        end
        
        function compute_ePtP(obj)
            if strcmp(obj.methodPesti,'vonmises')
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
            obj.eAiDFtPtPFDAi = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD,1);
        end
        
        function value = computeXInnerProduct(obj)
            value = zeros(1,obj.data.K);
            for k = 1:obj.data.K
                value(k) = sum(sum(obj.data.X(:,:,k).^2));
            end
        end
        
        % #################################################################
        % # Utility Functions
        
        function C = matrixProductPrSlab(obj,A,B)
            % Computes the matrix product C = A*B with the defined method
           
            switch obj.methodMatrixProduct
                case 'gpu'
                    C = pagefun(@mtimes,A,B);
                    % C = gather(C);
                    
                case 'naive'
                    K = max(size(A,3),size(B,3));
                    C = zeros(size(A,1),size(B,2),K);
                    for k = 1:K
                        if size(A,3) > 1 && size(B,3) > 1
                            C(:,:,k) = A(:,:,k)*B(:,:,k);
                        elseif size(A,3) > 1
                            C(:,:,k) = A(:,:,k)*B;
                        else
                            C(:,:,k) = A*B(:,:,k);
                        end
                            
                    end
                    
                case 'mtimesx'
                    C = mtimesx(A,B);
                    
                case 'mmx'
                    C = mmx('mult',A,B);
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





