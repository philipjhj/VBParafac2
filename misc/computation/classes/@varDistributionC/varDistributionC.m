classdef varDistributionC < handle
    properties %(Access = private)
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
        
        eD
        eA
        % Options
        opts    % Struct with options
    end
    
    properties (Dependent)
        % ELBO Terms
        ELBO
        ePxz
        eQz
    end
    
    properties (Access = private)
        
        % Computed values
        XInnerProductPrSlab
        qPvonmisesEntropy = 0
        
        % Shared terms between moments
        
        eFtPt
        eDeAtXk
        eDeAtXkeFtPtTrace
        eAtA
        eCsquared
        eDAtAD
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
        
        orderParams = {'qP','qF','qC','qA','qSigma','qAlpha'}
    end
    
    methods
        function obj = varDistributionC(modelobj,data)
            % Class controling the full variational distribution
            % Belongs to a 'varBayesModelParafac2' object
            
            obj.data = data;
            obj.opts = modelobj.opts;
            obj.util = modelobj.util;
        end
        
        function obj = initializeVariationalDististribution(obj)
            obj.initializeVariationalFactors;
            obj.initializeSufficientStatistics;
        end
        
        function initializeVariationalFactors(obj)
            obj.qA = multiNormalDist('qA',[obj.data.I obj.data.M],true,obj.util);
            
            obj.qC = multiNormalDist('qC',[obj.data.K obj.data.M],false,obj.util);
            obj.qF = multiNormalDist('qF',[obj.data.M obj.data.M],false,obj.util);
            obj.qP = multiNormalDist('qP',[obj.data.J obj.data.M obj.data.K],true,obj.util);
            obj.qSigma = GammaDist('qSigma',[1 obj.data.K]);
            obj.qAlpha = GammaDist('qAlpha',[1 obj.data.M]);
            
            obj.pSigma = GammaDist('pSigma',[1 1]);
            obj.pAlpha = GammaDist('pAlpha',[1 1]);
        end
        function initializeSufficientStatistics(obj)
            % Initialize Shared Terms
            obj.XInnerProductPrSlab = obj.compute_XInnerProductPrSlab;
            obj.updateSufficientStatistics;
            
            % Initialize ELBO terms
            allVariationalFactors = {'qA','qC','qF','qP','qSigma','qAlpha'};
            
            obj.updateStatistics(allVariationalFactors)
            obj.computeMeanLogValues(allVariationalFactors)
            obj.computeEntropyValues(allVariationalFactors)
        end
        
        % TODO: Refactor SufficientStatistics into seperate class?
        function updateSufficientStatistics(obj)
            obj.compute_eA;
            obj.compute_eAtA;
            obj.compute_eD;
            obj.compute_eFtPt;
            obj.compute_eDeAtXk;
            obj.compute_eDeAtXkeFtPtTrace;
            obj.compute_eCtC;
            obj.compute_eDAtAD;
            obj.compute_eCsquared;
            obj.compute_ePtP;
            obj.compute_eFtPtPF;
            obj.compute_eDFtPtPFD;
            obj.compute_eAiDFtPtPFDAi;
        end
        function updateSufficientStatisticsConditional(obj,VariationalFactorName)
            % Only updates shared terms depending on current updated
            % variational factor to minimize number of computations
            % TODO: refactor updates into
            if strcmp(VariationalFactorName,'qP')
                obj.compute_ePtP;
                
                if ~ismember('qF',obj.opts.activeParams)
                    obj.compute_eFtPtPF;
                    obj.compute_eFtPt;
                elseif all(~ismember({'qC','qF'},obj.opts.activeParams))
                    obj.compute_eDFtPtPFD;
                elseif all(~ismember({'qC','qF','qA'},obj.opts.activeParams))
                    obj.compute_eAiDFtPtPFDAi;
                end
                
            elseif strcmp(VariationalFactorName,'qF')
                obj.compute_eFtPtPF;
                obj.compute_eFtPt;
                
                if ~ismember('qC',obj.opts.activeParams)
                    obj.compute_eDFtPtPFD;
                elseif all(~ismember({'qC','qA'},obj.opts.activeParams))
                    obj.compute_eAiDFtPtPFDAi;
                end
                
            elseif strcmp(VariationalFactorName,'qC')
                obj.compute_eD;
                obj.compute_eCtC;
                obj.compute_eCsquared;
                obj.compute_eDFtPtPFD;
                
                if ~ismember('qA',obj.opts.activeParams)
                    obj.compute_eA;
                    obj.compute_eAiDFtPtPFDAi;
                    obj.compute_eDeAtXk;
                end
                
            elseif strcmp(VariationalFactorName,'qA')
                obj.compute_eA;
                obj.compute_eAiDFtPtPFDAi;
                obj.compute_eDeAtXk;
                obj.compute_eDeAtXkeFtPtTrace; % Term for qSigma update
            end
        end
        
        % #################################################################
        % # ELBO computations
        
        function value = get.ELBO(obj)
            value = obj.ePxz+obj.eQz;
        end
        
        function value = get.ePxz(obj)
            obj.computeMeanLogValues(obj.opts.activeParams)
            
            value = obj.qXMeanLog+obj.qAMeanLog+obj.qCMeanLog+...
                obj.qFMeanLog+obj.qPMeanLog+obj.qSigmaMeanLog+...
                obj.qAlphaMeanLog;
        end
        function value = get.eQz(obj)
            obj.computeEntropyValues(obj.opts.activeParams)
            
            value = obj.qAEntropy+obj.qCEntropy+obj.qFEntropy+...
                obj.qPEntropy+obj.qSigmaEntropy+obj.qAlphaEntropy;
        end
        
        function computeMeanLogValues(obj,variationalFactorNames)
            for i  = 1:numel(variationalFactorNames)
                methodStr = strcat('compute',variationalFactorNames{i},'MeanLog');
                obj.(methodStr);
            end
            
            obj.computeMeanLogSpecialCases(variationalFactorNames);
        end
        function computeMeanLogSpecialCases(obj,variationalFactorNames)
            
            % pXmeanLog Depends on everything but qAlpha
            if any(~ismember(variationalFactorNames,'qAlpha'))
                obj.computepXmeanLog
            end
            
            % C depends on Alpha, so update even if C is not active
            if ~ismember('qC',variationalFactorNames) && ismember('qAlpha',variationalFactorNames)
                obj.computeqCMeanLog;
            end
            
            % No expected value if hyperparameter are maximized
            if strcmp(obj.opts.estimationNoise,'max')
                obj.qSigmaMeanLog = 0;
            end
            if strcmp(obj.opts.estimationARD,'max')
                obj.qAlphaMeanLog = 0;
            end
        end
        
        function computeEntropyValues(obj,variationalFactorNames)
            for i  = 1:numel(variationalFactorNames)
                methodStr = strcat(variationalFactorNames{i},'Entropy');
                obj.(methodStr) = obj.(variationalFactorNames{i}).entropy;
            end
            
            obj.computeEntropySpecialCases;
        end
        function computeEntropySpecialCases(obj)
            if strcmp(obj.opts.estimationP,'vonmises')
                obj.qPEntropy = obj.qPvonmisesEntropy;
            end
            
            if strcmp(obj.opts.estimationNoise,'max')
                obj.qSigmaEntropy = 0;
            end
            if strcmp(obj.opts.estimationARD,'max')
                obj.qAlphaEntropy = 0;
            end
        end
        
        % #################################################################
        % expectations for ELBO
        % TODO: Refactor this code to more general functions
        
        function computepXmeanLog(obj)
            temp = obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(...
                obj.util.matrixProductPrSlab(obj.data.X,obj.qP.mean),...
                obj.qF.mean),obj.eD);
            
            tempSum = sum(sum(sum(obj.util.hadamardProductPrSlab(temp,...
                obj.util.transformToTensor(obj.qSigma.mean)),3).*obj.eA));
            
            obj.qXMeanLog = obj.data.J*obj.data.I/2*sum(obj.qSigma.MeanLog)-...
                1/2*sum(obj.qSigma.mean.*(...
                sum(obj.eAiDFtPtPFDAi)+obj.XInnerProductPrSlab))+tempSum;
        end
        function computeqAMeanLog(obj)
            obj.qAMeanLog = -1/2*obj.qA.meanInnerProductSumComponent;
        end
        function computeqCMeanLog(obj)
            if isempty(obj.qAlpha.entropy) && strcmp(obj.opts.estimationARD,'avg')
                obj.qAlpha.updateStatistics;
            end
            obj.qCMeanLog = obj.data.K/2*sum(obj.qAlpha.MeanLog)-1/2*(...
                trace(diag(obj.qAlpha.mean)*sum(obj.qC.variance,3))+...
                sum(sum(obj.qC.mean.^2*diag(obj.qAlpha.mean))));
            
        end
        function computeqFMeanLog(obj)
            obj.qFMeanLog = -1/2*obj.qF.meanInnerProductSumComponent;
        end
        function computeqPMeanLog(obj)
            obj.qPMeanLog = -1/2*obj.qP.meanInnerProductSumComponent;
        end
        function computeqSigmaMeanLog(obj)
            obj.qSigmaMeanLog = sum((obj.pSigma.alpha-1).*...
                obj.qSigma.MeanLog-obj.qSigma.mean.*1/obj.pSigma.beta);
        end
        function computeqAlphaMeanLog(obj)
            obj.qAlphaMeanLog = sum((obj.pAlpha.alpha-1)*...
                obj.qAlpha.MeanLog-obj.qAlpha.mean.*1/obj.pAlpha.beta);
        end
        
        % #################################################################
        % # Moment Updates
        
        function updateMoments(obj)
            obj.correctVariationalFactorUpdateOrder; %TODO: should be in a set function for activeParams
            
            for i = 1:numel(obj.opts.activeParams)
                if obj.opts.debugFlag > 1
                    ELBO_prev = obj.ELBO;
                end
                
                obj.updateVariationalFactor(obj.opts.activeParams{i})
                obj.updateSufficientStatisticsConditional(obj.opts.activeParams{i})
                
                if obj.opts.debugFlag > 1
                    obj.updateSufficientStatistics;
                    ELBO_new = obj.ELBO;
                    obj.checkELBOconvergenceOnUpdate(obj.opts.activeParams{i},ELBO_prev,ELBO_new)
                end
            end
        end
        
        function correctVariationalFactorUpdateOrder(obj)
            if obj.testVariationalFactorUpdateOrder
                obj.updateVariationalFactorUpdateOrder;
            end
        end
        function bool=testVariationalFactorUpdateOrder(obj)
            bool=~isequal(obj.orderParams,obj.opts.activeParams);
        end
        function updateVariationalFactorUpdateOrder(obj)
            obj.opts.activeParams = obj.orderParams(...
                ismember(obj.orderParams,obj.opts.activeParams));
        end
        
        function updateVariationalFactor(obj,variationalFactorName)
            if ~ismember(variationalFactorName,{'qAlpha','qSigma'}) || obj.data.iter>=25
                obj.(strcat('update',variationalFactorName));
                obj.updateStatistics({variationalFactorName})
            end
        end
        
        function updateStatistics(obj,variationalFactorNames)
            for i  = 1:numel(variationalFactorNames)
                obj.(variationalFactorNames{i}).updateStatistics;
            end
            
            obj.updateStatisticsSpecialCases(variationalFactorNames);
        end
        function updateStatisticsSpecialCases(obj,variationalFactorNames)
            if ismember('qAlpha',variationalFactorNames) && ...
                    strcmp(obj.opts.estimationARD,'max')
                
                obj.qAlpha.mean = obj.data.K./sum(obj.eCsquared,1);
                obj.qAlpha.MeanLog = log(obj.qAlpha.mean);
            end
            
            if ismember('qSigma',variationalFactorNames) && ...
                    strcmp(obj.opts.estimationNoise,'max')
                
                obj.qSigma.mean = 1./(1/(obj.data.J*obj.data.I)*(sum(obj.eAiDFtPtPFDAi,1)+...
                    obj.XInnerProductPrSlab-...
                    2*squeeze(sum(sum(...
                    obj.util.hadamardProductPrSlab(...
                    obj.util.matrixProductPrSlab(obj.eA,...
                    obj.util.matrixProductPrSlab(obj.eD,...
                    obj.util.matrixProductPrSlab(obj.qF.mean',...
                    permute(obj.qP.mean,[2 1 3]))))...
                    ,obj.data.X),1),2))'));
                obj.qSigma.MeanLog = log(obj.qSigma.mean);
            end
        end
        
        function checkELBOconvergenceOnUpdate(obj,variationalFactorName,ELBO_prev,ELBO_new)
            if (ELBO_new-ELBO_prev)/abs(ELBO_new) < -1e-12 && obj.data.iter > 1
                warning('off','backtrace')
                warning(sprintf('varDist:Update:%s',variationalFactorName),...
                    'Update problem; ELBO absolute/relative change for %s update is %f \t %d\n',variationalFactorName,...
                    ELBO_new-ELBO_prev,(ELBO_new-ELBO_prev)/abs(ELBO_new));
                warning('on','backtrace')
                obj.data.errorIters = [obj.data.errorIters obj.data.iter];
                obj.data.errorIters_parameter = [obj.data.errorIters_parameter {variationalFactorName}];
            end
            obj.data.ELBOall = [obj.data.ELBOall obj.ELBO];
        end
        
        % ## Normal distribution updates
        function updateqA(obj)
            obj.qA.variance = inv(sum(obj.util.hadamardProductPrSlab(...
                obj.eDFtPtPFD,obj.util.transformToTensor(obj.qSigma.mean)),3)...
                +eye(obj.data.M));
            
            sum_k = sum(obj.util.hadamardProductPrSlab(...
                obj.util.transformToTensor(obj.qSigma.mean),...
                obj.util.matrixProductPrSlab(...
                obj.util.matrixProductPrSlab(obj.eD,obj.eFtPt),...
                permute(obj.data.X,[2 1 3]))),3);
            obj.qA.mean = (obj.qA.variance*sum_k)';
        end
        function updateqC(obj)
            obj.qC.variance=obj.util.matrixInversePrSlab(...
                bsxfun(@plus,obj.util.hadamardProductPrSlab(...
                obj.util.transformToTensor(obj.qSigma.mean),...
                obj.util.hadamardProductPrSlab(obj.eAtA,...
                obj.eFtPtPF))...
                ,diag(obj.qAlpha.mean)));
            
            obj.qC.mean=squeeze(obj.util.matrixProductPrSlab(...
                obj.util.hadamardProductPrSlab(...
                obj.util.transformToTensor(obj.qSigma.mean),...
                sum(obj.util.hadamardProductPrSlab(eye(obj.data.M),...
                obj.util.matrixProductPrSlab(obj.eFtPt,...
                obj.util.matrixProductPrSlab(permute(obj.data.X,[2 1 3]),...
                obj.eA))),1)),...
                obj.qC.variance))';
        end
        function updateqF(obj)
            
            p = obj.data.M;
            n = obj.data.K;
            
            % the bsxfun as index for ePtP gets the diagonals
            ePtPdiag = reshape(obj.ePtP(bsxfun(@plus,[1:p+1:p*p]',[0:n-1]*p*p))',1,1,n,p);
            
            varInv = bsxfun(@plus,squeeze(sum(...
                obj.util.hadamardProductPrSlab(...
                obj.util.hadamardProductPrSlab(...
                obj.util.transformToTensor(obj.qSigma.mean),...
                obj.eDAtAD),...
                ePtPdiag),3)),eye(obj.data.M));
            
            obj.qF.variance = obj.util.matrixInversePrSlab(varInv);
            
            t2=squeeze(sum(obj.util.hadamardProductPrSlab(...
                obj.util.transformToTensor(obj.qSigma.mean),...
                obj.util.matrixProductPrSlab(permute(obj.qP.mean,[2 1 3]),...
                permute(obj.eDeAtXk,[2 1 3]))),3));
            
            ePtProws = reshape(permute(obj.ePtP,[2 1 3]),1,obj.data.M,obj.data.M,obj.data.K);
            
            for m = 1:obj.data.M
                allButM = repmat(1:obj.data.M,obj.data.M,1)~=m;
                
                a = obj.util.hadamardProductPrSlab(ePtProws,obj.qF.mean');
                b = obj.util.hadamardProductPrSlab(a,allButM);
                
                c = squeeze(sum(b,2));
                
                t1=sum(obj.util.hadamardProductPrSlab(...
                    obj.util.transformToTensor(obj.qSigma.mean),...
                    obj.util.matrixProductPrSlab(...
                    obj.eDAtAD,...
                    c(:,m,:))),3)';
                
                obj.qF.mean(m,:) = (t2(m,:)-t1)*obj.qF.variance(:,:,m);
            end
        end
        function updateqP(obj)
            if ~strcmp(obj.opts.estimationP,'vonmises')
                %
                obj.qP.variance = permute( obj.util.matrixInversePrSlab(bsxfun(@plus,obj.util.hadamardProductPrSlab(...
                    obj.util.transformToTensor(obj.qSigma.mean),...
                    (obj.qF.computeMeanInnerProductScaledSlabs(...
                    permute(obj.eDAtAD,[1 2 4 3])...
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
                    
                    gradconstant=(obj.qF.mean*obj.eD(:,:,k)*obj.eA'*obj.data.X(:,:,k))';
                    costconstant=obj.qF.mean*obj.eD(:,:,k)*obj.eA'*obj.data.X(:,:,k);%*obj.qP.mean(:,:,k);
                    
                    % cost = costFunc(obj.qP.mean(:,:,k));
                    
                    [x,~] = trustregions(problem,[],options);
                    % fprintf('\n Slab %d with cost diff: %.2f \n',k,cost-xcost)
                    obj.qP.mean(:,:,k) = x;
                    
                end
                
            elseif strcmp(obj.opts.estimationP,'vonmises')
                obj.qPvonmisesEntropy = 0;
                
                F=obj.util.hadamardProductPrSlab(obj.util.matrixProductPrSlab(obj.qF.mean,obj.eDeAtXk),...
                    obj.util.transformToTensor(obj.qSigma.mean));
                
                for k=1:obj.data.K
                    [UU,SS,VV]=svd(F(:,:,k),'econ');
                    [f,~,lF]=hyperg(obj.data.J,diag(SS),3);
                    E_Z=UU*diag(f)*VV';                 % Expectation of Z
                    H_Z= lF-sum(sum(F(:,:,k).*E_Z));           % Entropy
                    obj.qP.mean(:,:,k) = E_Z';
                    obj.qPvonmisesEntropy = obj.qPvonmisesEntropy+H_Z;
                end
                
            elseif strcmp(obj.opts.estimationP,'parafac2svd')
                qFeDeAtX=obj.util.matrixProductPrSlab(obj.qF.mean,...
                    obj.eDeAtXk);
                
                for k = 1:obj.data.K
                    [U,~,V] = svd(qFeDeAtX(:,:,k),'econ');
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
        
        % ## Gamma distribution updates
        function updateqSigma(obj)
            if strcmp(obj.opts.estimationNoise,'avg')
                obj.qSigma.alpha(:) = obj.pSigma.alpha+obj.data.I*obj.data.J/2;
                
                obj.qSigma.beta = 1./(1./obj.pSigma.beta+1/2*...
                    sum(obj.eAiDFtPtPFDAi,1)+1/2*obj.XInnerProductPrSlab-...
                    obj.eDeAtXkeFtPtTrace);
                
            elseif strcmp(obj.opts.estimationNoise,'max')
                obj.qSigma.mean = 1./(1/(obj.data.J*obj.data.I)*...
                    (sum(obj.eAiDFtPtPFDAi,1)+...
                    obj.XInnerProductPrSlab-...
                    2*obj.eDeAtXkeFtPtTrace));
                
                obj.qSigma.MeanLog = log(obj.qSigma.mean);
            end
        end
        function updateqAlpha(obj)
            if strcmp(obj.opts.estimationARD,'avg')
                obj.qAlpha.alpha(:) = obj.pAlpha.alpha+1/2*obj.data.K;
                obj.qAlpha.beta = 1./(1/obj.pAlpha.beta+1/2*sum(obj.eCsquared,1));
            elseif strcmp(obj.opts.estimationARD,'max')
                obj.qAlpha.mean = obj.data.K./sum(obj.eCsquared,1);
                obj.qAlpha.MeanLog = log(obj.qAlpha.mean);
            end
        end
        
        % #################################################################
        % # Shared Terms
        
        % ## First order
        function compute_eD(obj)
            obj.eD = obj.util.matrixDiagonalPrSlab(obj.qC.mean');
        end
        
        function compute_eA(obj,r0)
            % Expectation of A w.r.t. mode r0
            obj.eA = obj.qA.mean;
        end
        
        function compute_eFtPt(obj)
            obj.eFtPt = obj.util.matrixProductPrSlab(obj.qF.mean',...
                permute(obj.qP.mean,[2 1 3]));
        end
        
        function compute_eDeAtXk(obj)
            obj.eDeAtXk = obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(...
                obj.eD,obj.eA'),...
                obj.data.X);
        end
        
        function compute_eDeAtXkeFtPtTrace(obj)
            % Trace of XkteAeDkeFtePt computed by a Hadamard product
            obj.eDeAtXkeFtPtTrace=squeeze(sum(sum(obj.util.hadamardProductPrSlab(obj.eDeAtXk,obj.eFtPt),1),2))';
        end
        
        % ## Second or Higher Order
        function compute_eAtA(obj)
           obj.eAtA = obj.qA.meanOuterProduct;
        end
        
        function compute_eCsquared(obj)
            obj.eCsquared = obj.qC.distSquared;
        end
        
        function compute_eCtC(obj)
            obj.eCtC = obj.qC.meanOuterProductPrVec;
        end
        
        function compute_eDAtAD(obj)
           obj.eDAtAD = obj.util.hadamardProductPrSlab(obj.eCtC,obj.eAtA); 
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
            ePtPdiag = reshape(sum(...
                obj.util.hadamardProductPrSlab(eye(obj.data.M),obj.ePtP),2),...
                1,1,obj.data.M,obj.data.K);
            
            value=squeeze(sum(sum(...
                obj.util.hadamardProductPrSlab(fullF,ePtPtensor),3),4))+...
                squeeze(sum(obj.util.hadamardProductPrSlab(obj.qF.variance,ePtPdiag),3));
            
            obj.eFtPtPF = value;
        end
        
        function compute_eDFtPtPFD(obj)
            value = obj.util.hadamardProductPrSlab(obj.eCtC,obj.eFtPtPF);
            obj.eDFtPtPFD = value;
        end
        
        function compute_eAiDFtPtPFDAi(obj)
            obj.eAiDFtPtPFDAi = obj.qA.computeMeanInnerProductScaledSlabs(obj.eDFtPtPFD,1);
        end
        
        function value = compute_XInnerProductPrSlab(obj)
            value = squeeze(sum(sum(obj.data.X.^2,1),2))';
        end
        
        % #################################################################
        % # Statistics methods
        function nActive = nActiveComponents(obj,method)
            
            if nargin<2
                method = obj.opts.nActiveComponents;
            end
            
            if strcmp(method,'hard')
                if isa(obj.qC.mean,'gpuArray')
                    nActive = sum(sum(gather(obj.qC.mean))~=0);
                else
                    nActive = sum(sum(obj.qC.mean)~=0);
                end
            elseif strcmp(method,'threshold')
                nActive = find(cumsum(sort(1./obj.qAlpha.mean,'descend')/sum(1./obj.qAlpha.mean))>0.95,1);
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
