classdef varDistributionC < handle
    properties (Access = private)
        % Data reference property
        data
        util
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
        
        % Options
        opts    % Struct with options
        % see constructor for possible options
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
        orderParams = {'qP','qF','qC','qA','qSigma','qAlpha'}
    end
    
    methods
        function obj = varDistributionC(modelobj)
            % Class to controls all the variational factors and their
            % operations.
            % Constructed by a 'varBayesModelParafac2' object
            
            % Initialize Data
            obj.data = modelobj.data;
            obj.opts = modelobj.opts;
            obj.util = modelobj.util;
            
            % Initialize distribution classes
            obj.qA = multiNormalDist('qA',[obj.data.I obj.data.M],true,obj.util);
            obj.qC = multiNormalDist('qC',[obj.data.K obj.data.M],false,obj.util);
            obj.qF = multiNormalDist('qF',[obj.data.M obj.data.M],false,obj.util);
            obj.qP = multiNormalDist('qP',[obj.data.J obj.data.M obj.data.K],true,obj.util);
            obj.qSigma = GammaDist('qSigma',[1 obj.data.K]);
            obj.qAlpha = GammaDist('qAlpha',[1 obj.data.M]);
            
            obj.pSigma = GammaDist('pSigma',[1 1]);
            obj.pAlpha = GammaDist('pAlpha',[1 1]);
            
        end
        
        
        function obj = initDist(obj)
            
            if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                all_params = {'qA','qC','qF','qP'};
                for i  = 1:numel(all_params)
                    obj.(all_params{i}).mean = gpuArray(obj.(all_params{i}).mean);
                    obj.(all_params{i}).variance = gpuArray(obj.(all_params{i}).variance);
                end
                obj.data.X = gpuArray(obj.data.X);
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
            obj.qPvonmisesEntropy=0;
            % Uncomment to set values to the ground truth
            %                         obj.qA.mean = [obj.data.Atrue zeros(obj.data.I,obj.data.M-obj.data.Mtrue)];
            %                         obj.qC.mean = [obj.data.Ctrue zeros(obj.data.K,obj.data.M-obj.data.Mtrue)];
            %                                     obj.qF.mean = obj.data.Ftrue;
            %                         obj.qP.mean = obj.data.Ptrue;
            %                                     obj.qAlpha.mean = [obj.data.Alphatrue 1e9*ones(1,obj.data.M-obj.data.Mtrue)];
            %                                     obj.qSigma.mean = obj.data.Sigmatrue;
            
            
            
        end
        
        
        %
        %         function set.opts.activeParams(obj,value)
        %             all_params = {'qP','qF','qC','qA','qSigma','qAlpha'};
        %             obj.opts.activeParams = value;
        %             obj.activeParams = all_params(ismember(all_params,obj.opts.activeParams));
        %         end
        %
        %         function set.opts.MatrixProduct(obj,value)
        %            obj.methodMatrixProduct = value;
        %            obj.util.opt_matrixProductPrSlab = obj.methodMatrixProduct;
        %         end
        %
        %         function set.opts(obj,value)
        %             obj = obj;
        %             disp(value);
        %         end
        %
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
            if isa(value,'gpuArray')
                value = gather(value);
            end
        end
        
        % ## Mean terms
        function value = get.ePxz(obj)
            % Recompute for active parameters
            for i  = 1:numel(obj.opts.activeParams)
                methodStr = strcat('compute',obj.opts.activeParams{i},'MeanLog');
                obj.(methodStr);
                
            end
            
            if any(~ismember(obj.opts.activeParams,'qAlpha'))
                obj.computeqXMeanLog % Depends on everything but qAlpha
            end
            
            % C depends on Alpha, so update even if C is not active
            if ismember('qAlpha',obj.opts.activeParams) && ~ismember('qC',obj.opts.activeParams)
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
            
            for i = 1:numel(obj.opts.activeParams)
                methodStr = strcat(obj.opts.activeParams{i},'Entropy');
                
                obj.(methodStr) = obj.(obj.opts.activeParams{i}).entropy;
            end
            
            if strcmp(obj.opts.estimationP,'vonmises')
                obj.qPEntropy = obj.qPvonmisesEntropy;
            end
            
            % Compute sum of entropies
            value = obj.qAEntropy+obj.qCEntropy+obj.qFEntropy+...
                obj.qPEntropy+obj.qSigmaEntropy+obj.qAlphaEntropy;
        end
        
        % #################################################################
        % Mean values for ELBO
        function computeqXMeanLog(obj)
            
            %t3 = zeros(obj.data.I,obj.data.M,obj.data.K);
            
            %  if ~strcmp(obj.opts.matrixProductPrSlab,'gpu')
            dataX = obj.data.X;
            qPMean = obj.qP.mean;
            qFMean = obj.qF.mean;
            qDMean = obj.eD;
            %   else
            %	dataX = gpuArray(obj.data.X);
            %     qPMean = gpuArray(obj.qP.mean);
            %     qFMean = gpuArray(obj.qF.mean);
            %     qDMean = gpuArray(obj.eD);
            % end
            
            
            t3 = obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(dataX,qPMean),qFMean),qDMean);
            
            t3sum = sum(sum(sum(bsxfun(@times,t3,reshape(obj.qSigma.mean,1,1,obj.data.K)),3).*obj.qA.mean));
            
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
            
            % Check order of updates
            if ~all(cellfun(@strcmp,obj.orderParams(...
                    ismember(obj.orderParams,obj.opts.activeParams)),...
                    obj.opts.activeParams))
                obj.opts.activeParams = obj.orderParams(ismember(obj.orderParams,obj.opts.activeParams));
            end
            
            any_error = 0;
            first_error = 0;
            
            for i = 1:numel(obj.opts.activeParams)
                %                 ELBO_prev = obj.ELBO;
                
                %                 if any(~ismember(obj.opts.activeParams,'qAlpha'))
                %                     XMeanLog_prev = obj.qXMeanLog; % Depends on everything but qAlpha
                %                 else
                %                     XMeanLog_prev= 0;
                %                 end
                
                %                 if strcmp(obj.opts.activeParams{i},'qAlpha')
                %                     CMeanLog_prev = obj.qCMeanLog; % Depends on everything but qAlpha
                %                 else
                %                     CMeanLog_prev= 0;
                %                 end
                
                %                 methodStr = strcat(obj.opts.activeParams{i},'MeanLog');
                
                %                 ELBO_prev = XMeanLog_prev+CMeanLog_prev+...
                %                     obj.(methodStr)+obj.(obj.opts.activeParams{i}).entropy;
                
                %                 ELBO_prev2 = obj.ELBO;
                
                if ~ismember(obj.opts.activeParams{i},{'qAlpha','qSigma'}) || obj.data.iter>=25% || obj.data.iter==0
                    
                    methodStr = strcat('update',obj.opts.activeParams{i});
                    obj.(methodStr);
                    obj.(obj.opts.activeParams{i}).updateStatistics;
                end
                
                
                
                % Update shared terms
                if strcmp(obj.opts.activeParams{i},'qP')
                    obj.compute_ePtP;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qF',obj.opts.activeParams)
                        obj.compute_eFtPtPF;
                    elseif all(~ismember({'qC','qF'},obj.opts.activeParams))
                        obj.compute_eDFtPtPFD;
                    elseif all(~ismember({'qC','qF','qA'},obj.opts.activeParams))
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.opts.activeParams{i},'qF')
                    obj.compute_eFtPtPF;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qC',obj.opts.activeParams)
                        obj.compute_eDFtPtPFD;
                    elseif all(~ismember({'qC','qA'},obj.opts.activeParams))
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.opts.activeParams{i},'qC')
                    obj.compute_eD;
                    obj.compute_eCtC;
                    obj.compute_eCsquared;
                    obj.compute_eDFtPtPFD;
                    
                    % Update dependent terms if some not active
                    if ~ismember('qA',obj.opts.activeParams)
                        obj.compute_eAiDFtPtPFDAi;
                    end
                    
                elseif strcmp(obj.opts.activeParams{i},'qA')
                    obj.compute_eAiDFtPtPFDAi;
                end
                
                
                %                 if any(~ismember(obj.opts.activeParams,'qAlpha'))
                %                     obj.computeqXMeanLog;
                %                     XMeanLog_new = obj.qXMeanLog; % Depends on everything but qAlpha
                %                 else
                %                     XMeanLog_new=0;
                %                 end
                %
                %                 if strcmp(obj.opts.activeParams{i},'qAlpha')
                %                     obj.computeqCMeanLog;
                %                     CMeanLog_new = obj.qCMeanLog; % Depends on everything but qAlpha
                %                 else
                %                     CMeanLog_new= 0;
                %                 end
                %
                %
                %                 methodStr = strcat('compute',obj.opts.activeParams{i},'MeanLog');
                %                 obj.(methodStr);
                %
                %                 methodStr = strcat(obj.opts.activeParams{i},'MeanLog');
                %
                %                 ELBO_new = XMeanLog_new+CMeanLog_new+...
                %                     obj.(methodStr)+obj.(obj.opts.activeParams{i}).entropy;
                %
                %                 ELBO_new2 = obj.ELBO;
                
                %                 if obj.opts.debugFlag > 1 && (ELBO_new-ELBO_prev)/abs(ELBO_new) < -1e-8 && obj.data.iter > 0
                %                     any_error = 1;
                %                     if ~first_error
                %                         first_error = 1;
                %                         fprintf('\n');
                %                     end
                %                     warning('off','backtrace')
                %                     warning(sprintf('varDist:Update:%s',obj.opts.activeParams{i}),...
                %                         'Update problem; ELBO absolute/relative change for %s update is %f / %f \t %f\n',obj.opts.activeParams{i},...
                %                         ELBO_new-ELBO_prev,(ELBO_new-ELBO_prev)/abs(ELBO_new),(ELBO_new2-ELBO_prev2)/abs(ELBO_new2));
                %                     warning('on','backtrace')
                %                 end
                
                obj.data.ELBOall = [obj.data.ELBOall obj.ELBO];
            end
            if any_error
                fprintf('\n')
            end
        end
        
        % ## Normal distributions
        % ### Variational Factor A
        function updateqA(obj)
            obj.qA.variance = inv(sum(bsxfun(@times,obj.eDFtPtPFD,reshape(obj.qSigma.mean,1,1,obj.data.K)),3)...
                +eye(obj.data.M));
            
            qSigmaMean = obj.qSigma.mean;
            qD = obj.eD;
            qFMeanT = obj.qF.mean';
            qPMeanT = permute(obj.qP.mean,[2 1 3]);
            dataXT = permute(obj.data.X,[2 1 3]);
            
            qAvariance = obj.qA.variance;
            
            K = obj.data.K;
            sum_k = sum(bsxfun(@times,reshape(qSigmaMean,1,1,obj.data.K),...
                obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(...
                obj.util.matrixProductPrSlab(qD,qFMeanT),qPMeanT),dataXT)),3);
            obj.qA.mean = (qAvariance*sum_k)';
            %
        end
        
        % ### Variational Factor C
        function updateqC(obj)
            
            obj.qC.variance=obj.util.matrixInversePrSlab(...
                bsxfun(@plus,bsxfun(@times,...
                reshape(obj.qSigma.mean,1,1,obj.data.K),...
                bsxfun(@times,obj.qA.meanOuterProduct,obj.eFtPtPF))...
                ,diag(obj.qAlpha.mean)));
            
            obj.qC.mean=squeeze(obj.util.matrixProductPrSlab(...
                bsxfun(@times,reshape(obj.qSigma.mean,1,1,obj.data.K),...
                sum(bsxfun(@times,eye(obj.data.M),...
                obj.util.matrixProductPrSlab(obj.qF.mean',...
                obj.util.matrixProductPrSlab(permute(obj.qP.mean,[2 1 3]),...
                obj.util.matrixProductPrSlab(permute(obj.data.X,[2 1 3]),...
                obj.qA.mean)))),1)),...
                obj.qC.variance))';
            
            
        end
        
        % ### Variational Factor F
        function updateqF(obj)
            
            p = obj.data.M;
            n = obj.data.K;
            
            % Below code replaces a loop over m.
            % the bsxfun as index for ePtP gets the diagonals
            
            
            if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                obj.ePtP = gather(obj.ePtP);
            end
            
            ePtPdiag = reshape(obj.ePtP(bsxfun(@plus,[1:p+1:p*p]',[0:n-1]*p*p))',1,1,n,p);
            
            varInv = gather(bsxfun(@plus,squeeze(sum(bsxfun(@times,bsxfun(@times,...
                reshape(obj.qSigma.mean,1,1,obj.data.K),...
                bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct)),...
                ePtPdiag),3)),eye(obj.data.M)));
            
            obj.qF.variance = multinv(varInv);
            
            
            
            t2=squeeze(sum(bsxfun(@times,reshape(obj.qSigma.mean,1,1,obj.data.K),...
                obj.util.matrixProductPrSlab(permute(obj.qP.mean,[2 1 3]),...
                obj.util.matrixProductPrSlab(permute(obj.data.X,[2 1 3]),...
                obj.util.matrixProductPrSlab(obj.qA.mean,obj.eD)))),3));
            
            t2 = gather(t2);
            
            %
            ePtProws = reshape(permute(obj.ePtP,[2 1 3]),1,obj.data.M,obj.data.M,obj.data.K);
            %allButM = bsxfun(@(A,B) A~=B,repmat(1:obj.data.M,obj.data.M,1),reshape(1:obj.data.M,1,1,obj.data.M));
            
            if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                obj.qF.mean = gather(obj.qF.mean);
                obj.qF.variance = gather(obj.qF.variance);
            end
            
            for m = 1:obj.data.M
                allButM = repmat(1:obj.data.M,obj.data.M,1)~=m;
                
                if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                    a = bsxfun(@times,ePtProws,gpuArray(obj.qF.mean'));
                else
                    a = bsxfun(@times,ePtProws,obj.qF.mean');
                end
                b = bsxfun(@times,a,allButM);
                
                if isa(b,'gpuArray')
                    c = gather(squeeze(sum(b,2)));
                else
                    c = squeeze(sum(b,2));
                end
                
                t1=sum(bsxfun(@times,reshape(obj.qSigma.mean,1,1,obj.data.K),...
                    mtimesx(gather(bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct)),...
                    c(:,m,:))),3)';
                
                if isa(t1,'gpuArray')
                    t1=gather(t1);
                end
                
                obj.qF.mean(m,:) = (t2(m,:)-t1)*obj.qF.variance(:,:,m);
            end
            
            if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                obj.qF.mean = gpuArray(obj.qF.mean);
                obj.qF.variance = gpuArray(obj.qF.variance);
            end
            
        end
        % ### Variational Factor P
        function updateqP(obj)
            if ~strcmp(obj.opts.estimationP,'vonmises')
                %
                obj.qP.variance = permute( obj.util.matrixInversePrSlab(bsxfun(@plus,bsxfun(@times,...
                    reshape(obj.qSigma.mean,1,1,obj.data.K),...
                    (obj.qF.computeMeanInnerProductScaledSlabs(...
                    permute(bsxfun(@times,obj.eCtC,obj.qA.meanOuterProduct),[1 2 4 3])...
                    )))...
                    ,eye(obj.data.M))),[1 2 4 3]);
            end
            obj.computeqPmean;
        end
        
        function computeqPmean(obj)
            if strcmp(obj.opts.estimationP,'manopt')
                
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
                
            elseif strcmp(obj.opts.estimationP,'vonmises')
                obj.qPvonmisesEntropy = 0;
                
                A = obj.util.matrixProductPrSlab(obj.qA.mean,...,
                    obj.util.matrixProductPrSlab(obj.eD,obj.qF.mean'));
                
                % Estep
                F=bsxfun(@times,obj.util.matrixProductPrSlab(permute(A,[2 1 3])...
                    ,obj.data.X),...
                    reshape(obj.qSigma.mean,1,1,obj.data.K)); %/sigma_sq;
                
                if isa(F,'gpuArray')
                    F = gather(F);
                end
                
                if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                    obj.qP.mean = gather(obj.qP.mean);
                end
                
                
                for k=1:obj.data.K
                    [UU,SS,VV]=svd(F(:,:,k),'econ');
                    [f,~,lF]=hyperg(obj.data.J,diag(SS),3);
                    E_Z=UU*diag(f)*VV';                 % Expectation of Z
                    H_Z= lF-sum(sum(F(:,:,k).*E_Z));           % Entropy
                    obj.qP.mean(:,:,k) = E_Z';
                    obj.qPvonmisesEntropy = obj.qPvonmisesEntropy+H_Z;
                end
                
                
                if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                    obj.qP.mean = gpuArray(obj.qP.mean);
                end
                
                
            elseif strcmp(obj.opts.estimationP,'parafac2svd')
                
                
                qFeDqAX=obj.util.matrixProductPrSlab(obj.qF.mean,...
                    obj.util.matrixProductPrSlab(obj.eD,...
                    obj.util.matrixProductPrSlab(obj.qA.mean',obj.data.X)));
                
                if isa(qFeDqAX,'gpuArray')
                    qFeDqAX = gather(qFeDqAX);
                end
                
                if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                    obj.qP.mean = gather(obj.qP.mean);
                end
                
                for k = 1:obj.data.K
                    [U,~,V] = svd(qFeDqAX(:,:,k),'econ');
                    
                    obj.qP.mean(:,:,k) = V(:,1:obj.data.M)*U';
                end
                
                if strcmp(obj.opts.matrixProductPrSlab,'gpu')
                    obj.qP.mean = gpuArray(obj.qP.mean);
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
            obj.qSigma.alpha(:) = obj.pSigma.alpha+obj.data.I*obj.data.J/2;
            
            obj.qSigma.beta = 1./(1./obj.pSigma.beta+1/2*sum(obj.eAiDFtPtPFDAi,1)...
                +1/2*obj.XInnerProduct...
                -squeeze(sum(sum(...
                bsxfun(@times,...
                obj.util.matrixProductPrSlab(obj.qA.mean,...
                obj.util.matrixProductPrSlab(obj.eD,...
                obj.util.matrixProductPrSlab(obj.qF.mean',...
                permute(obj.qP.mean,[2 1 3]))))...
                ,obj.data.X),1),2))');
            
            if isa(obj.qSigma.beta,'gpuArray')
                obj.qSigma.beta = gather(obj.qSigma.beta);
            end
            
        end
        
        
        
        % ### Variational Factor Alpha
        function updateqAlpha(obj)
            
            obj.qAlpha.alpha(:) = obj.pAlpha.alpha+1/2*obj.data.K;
            obj.qAlpha.beta = 1./(1/obj.pAlpha.beta+1/2*sum(obj.eCsquared,1));
            
            if isa(obj.qAlpha.beta,'gpuArray')
                obj.qAlpha.beta = gather(obj.qAlpha.beta);
            end
        end
        
        % #################################################################
        % # Terms for moment updates and means in the ELBO
        
        
        % ## First order
        
        function compute_eD(obj)
            
            obj.eD = bsxfun(@mtimes,reshape(obj.qC.mean',1,...
                obj.data.M,obj.data.K),...
                repmat(eye(obj.data.M),1,1,obj.data.K));
        end
        
        % ## Second or Higher Order
        function compute_eCsquared(obj)
            obj.eCsquared = obj.qC.distSquared;
        end
        
        function compute_ePtP(obj)
            if strcmp(obj.opts.estimationP,'vonmises')
                obj.ePtP = repmat(eye(obj.data.M),1,1,obj.data.K);
            else
                obj.ePtP = obj.data.J*squeeze(obj.qP.variance)+repmat(eye(obj.data.M),1,1,obj.data.K);%repmat(eye(obj.data.M),1,1,obj.data.K);
            end
        end
        
        function compute_eFtPtPF(obj)
            
            diagFvec=reshape(obj.qF.mean',1,obj.data.M,1,obj.data.M);
            diagFvecT = permute(diagFvec,[2 1 4 3]);
            fullF=(obj.util.matrixProductPrSlab(diagFvecT,diagFvec));
            
            ePtPtensor=reshape(obj.ePtP,1,1,obj.data.M,obj.data.M,obj.data.K);
            ePtPdiag = reshape(sum(bsxfun(@times,eye(obj.data.M),obj.ePtP),2),1,1,obj.data.M,obj.data.K);
            
            value=squeeze(sum(sum(bsxfun(@times,fullF,ePtPtensor),3),4))+squeeze(sum(bsxfun(@times,obj.qF.variance,ePtPdiag),3));
            
            obj.eFtPtPF = value;
        end
        
        function compute_eCtC(obj)
            obj.eCtC = obj.qC.meanOuterProductSingle;
        end
        
        function compute_eDFtPtPFD(obj)
            
            value = bsxfun(@times,obj.eCtC,obj.eFtPtPF);
            obj.eDFtPtPFD = value;
        end
        
        function compute_eAiDFtPtPFDAi(obj)
            obj.eAiDFtPtPFDAi = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD,1);
        end
        
        function value = computeXInnerProduct(obj)
            
            value = squeeze(sum(sum(obj.data.X.^2,1),2))';
            
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





