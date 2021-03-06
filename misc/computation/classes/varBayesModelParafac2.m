classdef varBayesModelParafac2 < parafac2BaseClass
    properties
        qDistTrain
        qDistTest
        
        fullData
        dataTest
        dataTrain
        
        cvRunsTest
        cvRunsTrain
        
        util
        opts
        
        CV_ELBOS_test
        CV_ELBOS_train
        % TODO: add statistics class
    end
    
     properties (Dependent)
        X
        A
        C
        D
        F
        P
        
        M
     end
    
    properties (Access = protected)
        currentPartition
        currentData
        currentqDist
        
        testParametersHeteroscedastic = {'qP','qC','qSigma'};
        testParametersHomoscedastic = {'qP','qC'};
        
        ListenerPartitionUpdates
    end
    
    events
        partitionUpdate
    end
    
    properties (Dependent)
        n_neg_ELBO_diff
        n_neg_ELBO_diff_idx
    end
    
    methods
        % Functions to extract the values for the variables of the base
        function value = get.X(obj)
            value = obj.fullData.Xunfolded;
        end
        function value = get.A(obj)
            value = obj.qDistTrain.qA.mean;
        end
        function value = get.C(obj)
            value = obj.qDistTrain.qC.mean;
        end
        function value = get.F(obj)
            value = obj.qDistTrain.qF.mean;
        end
        function value = get.P(obj)
            value = obj.qDistTrain.qP.mean;
        end
        function value = get.M(obj)
            value = obj.fullData.M;
        end
        
        function obj = varBayesModelParafac2(data,M, opts)
            if nargin < 3
                 obj.opts = optionsClass;
            else
                obj.opts = opts;
            end
            obj.util = utilitiesClass(obj);
            
            % TODO: Do not partition data here, gather full data first,
            % and partition in CV methods
            if isa(data,'double')
                % If data matrix input
                
            elseif isa(data,'dataClass')
                obj.fullData = data;
                obj.fullData.M = M;
            end
            
            obj.dataTrain = dataClass;
            obj.dataTest = dataClass;
            obj.dataTrain.M = M;
            obj.dataTest.M = M;
            obj.dataTrain.partitionName = 'Train';
            obj.dataTest.partitionName = 'Test';
            
            obj.qDistTrain = varDistributionC(obj.dataTrain, obj.opts, obj.util);
            obj.qDistTest = varDistributionC(obj.dataTest, obj.opts, obj.util);
            
            %obj.ListenerPartitionUpdates = addlistener(obj,'partitionUpdate',@obj.setCurrentPartitionNames);
        end
        
        function crossValidateM(obj,Minterval)
            M = numel(Minterval);
            
            T = 5;
            obj.cvRunsTest = cell(obj.fullData.K,T,M);
            obj.cvRunsTrain = cell(obj.fullData.K,T,M);
            obj.CV_ELBOS_train = zeros(obj.fullData.K,T,M);
            obj.CV_ELBOS_test = zeros(obj.fullData.K,T,M);
            for m = 1:M
                for k = 1:obj.fullData.K
                    
                    fprintf('[Start: %s]\t[M: %d]\t[Fold: %d/%d]\n',...
                        datestr(datetime('now')),Minterval(m),k,obj.fullData.K)
                    % init train T times
                    for t = 1:T
                        obj.dataTrain = dataClass;
                        obj.dataTrain.partitionName = 'Train';
                        obj.dataTrain.M = Minterval(m);
                        obj.partitionData(obj.fullData.X,k);
                        obj.qDistTrain = varDistributionC(obj,obj.dataTrain);
                        
                        %                         fprintf('[Train initialization: %d/%d]\n',t,T);
                        %                         obj.dataTrain.restartDataDiagnostics;
                        
                        obj.fitTrainingData;
                        
                        
                        obj.cvRunsTrain{k,t,m} = struct('qDist',obj.qDistTrain,'Data',obj.dataTrain);
                        obj.CV_ELBOS_train(k,t,m) = obj.qDistTrain.ELBO;
                        obj.dataTrain.Xunfolded = [];
                    end
                    
                    [~,ind]=max(obj.CV_ELBOS_train(k,:,m));
                    
                    obj.qDistTrain = obj.cvRunsTrain{k,ind,m}.qDist;
                    
                    % init test T times
                    for t = 1:T
                        obj.dataTest = dataClass;
                        obj.dataTest.partitionName = 'Test';
                        obj.dataTest.M = Minterval(m);
                        obj.partitionData(obj.fullData.X,k);
                        obj.qDistTest = varDistributionC(obj,obj.dataTest);
                        
                        %                         fprintf('[Test initialization: %d/%d]\n',t,T);
                        %                         obj.dataTest.restartDataDiagnostics;
                        obj.fitTestData;
                        
                        obj.cvRunsTest{k,t,m} = struct('qDist',obj.qDistTest,'Data',obj.dataTest);
                        obj.CV_ELBOS_test(k,t,m) = obj.qDistTest.ELBO;
                        obj.dataTest.Xunfolded = [];
                    end
                    %                         fprintf('\t[Fold: %f]\n',obj.CV_ELBOS(k,t,m))
                end
                
            end
        end
        
        function partitionData(obj,X,testSlices)
            dims = size(X);
            K=dims(end);
            nPrSlice = prod(dims(1:(end-1)));
            
            allSlices = 1:K;
            
            if nargin < 3
                testPct = 0;
                testSlices = randsample(K,floor(testPct*K));
            end
            trainSlices = allSlices(~ismember(allSlices,testSlices));
            
            X = reshape(X,nPrSlice,K);
            obj.dataTest.Xunfolded = reshape(X(:,testSlices),[dims(1:(end-1)) numel(testSlices)]);
            obj.dataTrain.Xunfolded = reshape(X(:,trainSlices),[dims(1:(end-1)) numel(trainSlices)]);
        end
        
        % Necessary? If yes, move to qDist
        function restartPartitions(obj)
            %             clear obj.qDistTrain;
            %             obj.qDistTrain = varDistributionC(obj);
            %                     obj.dataTrain.iter = 1;
            %
            %
            
        end
        
        function set.currentPartition(obj,name)
            obj.currentPartition = name;
            obj.setCurrentPartitionNames
            %notify(obj,'partitionUpdate');
        end
        function setCurrentPartitionNames(obj,~,~)
            obj.currentData = strcat('data',obj.currentPartition);
            obj.currentqDist = strcat('qDist',obj.currentPartition);
        end
        
        function fitTrainingData(obj)
            obj.currentPartition = obj.dataTrain.partitionName;
            
            rng(obj.opts.rngInput);
            obj.qDistTrain.initializeVariationalDististribution;
            
            obj.computeVarDistribution;
        end
        function fitTestData(obj)
            obj.currentPartition = obj.dataTest.partitionName;
            
            obj.qDistTest.initializeVariationalDististribution;
            
            trainParameters = obj.opts.activeParams;
            
            if ~isempty(strfind(obj.opts.estimationNoise,'Shared'))% Shared
                obj.opts.activeParams = obj.testParametersHomoscedastic(...
                    ismember(obj.testParametersHomoscedastic,obj.opts.activeParams));
            else
                obj.opts.activeParams = obj.testParametersHeteroscedastic(...
                    ismember(obj.testParametersHeteroscedastic,obj.opts.activeParams));
            end
            obj.qDistTest.qA = obj.qDistTrain.qA;
            obj.qDistTest.qF = obj.qDistTrain.qF;
            obj.qDistTest.qAlpha = obj.qDistTrain.qAlpha;
            
            obj.computeVarDistribution;
            obj.opts.activeParams = trainParameters;
        end
        
        function [stopReason,errorStatus] = computeVarDistribution(obj)
            stopReason = 0;
            errorStatus = 0;
            
            if obj.opts.verbose
                disp('Starting CAVI with:')
                disp(obj.(obj.currentqDist).opts.activeParams)
            end
            
            if obj.opts.verbose
                obj.displayHeader;
            end
            
            obj.initializeCAVI;
            obj.(obj.currentData).computeStartTimes;
            
            obj.(obj.currentData).ELBO = obj.(obj.currentqDist).ELBO;
            
            obj.(obj.currentData).ticCAVIwall = tic;
            obj.(obj.currentData).ticCAVIcpu = cputime;
            while obj.checkIfNotConvergenced
                obj.(obj.currentqDist).updateMoments;
                obj.(obj.currentData).ELBO = obj.(obj.currentqDist).ELBO;
                
                obj.storeDiagnostics;
                obj.debugCheckELBOconvergence;
                obj.displayResultsAll;
                
                obj.(obj.currentData).iter = obj.(obj.currentData).iter+1;
            end
            
            if obj.opts.verbose
                obj.displayHeader;
            end
            
            obj.printStopReason;
            obj.cleanUp;
        end
        
        function initializeCAVI(obj)
            if obj.(obj.currentData).iter == 1
                obj.(obj.currentData).n_components(obj.(obj.currentData).iter) = obj.(obj.currentqDist).nActiveComponents;
                obj.(obj.currentData).n_components_hard(obj.(obj.currentData).iter) = obj.(obj.currentqDist).nActiveComponents('hard');
                obj.(obj.currentData).evaltime(obj.(obj.currentData).iter) = 1e-10;
                obj.(obj.currentData).evaltime_cpu(obj.(obj.currentData).iter) = 1e-10;
            end
        end
        function bool = checkIfNotConvergenced(obj)
            bool = abs(obj.(obj.currentData).ELBO_diff)/abs(obj.(obj.currentData).ELBO) > obj.opts.tol && ...
                obj.opts.maxIter+1 > obj.(obj.currentData).iter && ...
                obj.opts.maxTime > obj.(obj.currentData).getLatestEvalTime || ...
                obj.(obj.currentData).iter<1.5*max(obj.opts.scaleLearningDelay,obj.opts.noiseLearningDelay);
        end
        
        function storeDiagnostics(obj)
            % Store progress information
            if isempty(obj.(obj.currentData).ELBO_chain) || numel(obj.(obj.currentData).ELBO_chain)<obj.(obj.currentData).iter+1
                obj.(obj.currentData).ELBO_chain = [obj.(obj.currentData).ELBO_chain zeros(1,100)];
                obj.(obj.currentData).fit_chain = [obj.(obj.currentData).fit_chain zeros(1,100)];
                obj.(obj.currentData).evaltime = [obj.(obj.currentData).evaltime zeros(1,100)];
                obj.(obj.currentData).evaltime_cpu = [obj.(obj.currentData).evaltime_cpu zeros(1,100)];
                obj.(obj.currentData).n_components = [obj.(obj.currentData).n_components zeros(1,100)];
                obj.(obj.currentData).n_components_hard = [obj.(obj.currentData).n_components_hard zeros(1,100)];
            end
            obj.(obj.currentData).ELBO_chain(obj.(obj.currentData).iter) = obj.(obj.currentData).ELBO;
            obj.(obj.currentData).fit_chain(obj.(obj.currentData).iter) = obj.Parafac2Fit(obj.(obj.currentqDist));
            obj.(obj.currentData).evaltime(obj.(obj.currentData).iter) = toc(obj.(obj.currentData).ticCAVIwall)+obj.(obj.currentData).startWallTime;
            obj.(obj.currentData).evaltime_cpu(obj.(obj.currentData).iter) = cputime-obj.(obj.currentData).ticCAVIcpu+obj.(obj.currentData).startCpuTime;
            obj.(obj.currentData).n_components(obj.(obj.currentData).iter) = obj.(obj.currentqDist).nActiveComponents;
            obj.(obj.currentData).n_components_hard(obj.(obj.currentData).iter) = obj.(obj.currentqDist).nActiveComponents('hard');
        end
        
        function debugCheckELBOconvergence(obj)
            if obj.opts.debugFlag >= 1 && obj.(obj.currentData).ELBO_diff/abs(obj.(obj.currentData).ELBO) < -1e-7 && obj.(obj.currentData).iter>0
                warning('off','backtrace')
                warning('At iter %d ELBO not converging; relativ diff. is %.10f, diff; %.4f \n',...
                    obj.(obj.currentData).iter,obj.(obj.currentData).ELBO_diff/abs(obj.(obj.currentData).ELBO),obj.(obj.currentData).ELBO_diff)
                warning('on','backtrace')
                
                fprintf('%5d',obj.(obj.currentData).iter);
                %                 obj.displayResultsAll(obj.(obj.currentData).ELBO_diff,obj.(obj.currentqDist));
            end
        end
        function printStopReason(obj)
            if obj.(obj.currentData).ELBO_diff/abs(obj.(obj.currentData).ELBO) < obj.opts.tol
                if obj.opts.verbose
                    fprintf('CAVI has converged with last change %f\n', abs(obj.(obj.currentData).ELBO_diff)/abs(obj.(obj.currentData).ELBO))
                end
                obj.(obj.currentData).stopReason = 'converged';
            elseif obj.opts.maxIter <= obj.(obj.currentData).iter
                if obj.opts.verbose
                    fprintf('CAVI has stopped at iteration %d (max iteration) with change %f\n',obj.(obj.currentData).iter-1,abs(obj.(obj.currentData).ELBO_diff)/abs(obj.(obj.currentData).ELBO))
                end
                obj.(obj.currentData).stopReason = 'maximum iteration';
            elseif obj.opts.maxTime <= obj.(obj.currentData).evaltime(obj.(obj.currentData).evaltime==max(obj.(obj.currentData).evaltime))
                if obj.opts.verbose
                    fprintf('CAVI has stopped after %f s. evaluation (max time) with change %f\n',obj.(obj.currentData).evaltime(obj.(obj.currentData).evaltime==max(obj.(obj.currentData).evaltime)),abs(obj.(obj.currentData).ELBO_diff)/abs(obj.(obj.currentData).ELBO))
                end
                obj.(obj.currentData).stopReason = 'maximum time';
            elseif obj.Parafac2Fit(obj.(obj.currentqDist))/100>(1-obj.opts.tol)
                if obj.opts.verbose
                    fprintf('CAVI has stopped with a fit of %f %% with change %f\n',obj.Parafac2Fit(obj.(obj.currentqDist)),abs(obj.(obj.currentData).ELBO_diff)/abs(obj.(obj.currentData).ELBO))
                end
                obj.(obj.currentData).stopReason = 4;
            end
        end
        function cleanUp(obj)
            % Trim preallocated memory
            obj.(obj.currentData).ELBO_chain = nonzeros(obj.(obj.currentData).ELBO_chain)';
            obj.(obj.currentData).fit_chain = nonzeros(obj.(obj.currentData).fit_chain)';
            obj.(obj.currentData).evaltime = nonzeros(obj.(obj.currentData).evaltime)';
            obj.(obj.currentData).evaltime_cpu = nonzeros(obj.(obj.currentData).evaltime_cpu)';
            obj.(obj.currentData).n_components = nonzeros(obj.(obj.currentData).n_components)';
            obj.(obj.currentData).n_components_hard = nonzeros(obj.(obj.currentData).n_components_hard)';
            
            % Gather results to CPU
            if strcmpi(obj.opts.matrixProductPrSlab,'gpu')
                obj.(obj.currentqDist).qA.mean = gather(obj.(obj.currentqDist).qA.mean);
                obj.(obj.currentqDist).qA.variance = gather(obj.(obj.currentqDist).qA.variance);
                obj.(obj.currentqDist).qC.mean = gather(obj.(obj.currentqDist).qC.mean);
                obj.(obj.currentqDist).qC.variance = gather(obj.(obj.currentqDist).qC.variance);
                obj.(obj.currentqDist).qF.mean = gather(obj.(obj.currentqDist).qF.mean);
                obj.(obj.currentqDist).qF.variance = gather(obj.(obj.currentqDist).qF.variance);
                obj.(obj.currentqDist).qP.mean = gather(obj.(obj.currentqDist).qP.mean);
                obj.(obj.currentqDist).qP.variance = gather(obj.(obj.currentqDist).qP.variance);
                obj.(obj.currentqDist).qSigma.alpha = gather(obj.(obj.currentqDist).qSigma.alpha);
                obj.(obj.currentqDist).qSigma.beta= gather(obj.(obj.currentqDist).qSigma.beta);
                obj.(obj.currentqDist).qAlpha.alpha = gather(obj.(obj.currentqDist).qAlpha.alpha);
                obj.(obj.currentqDist).qAlpha.beta= gather(obj.(obj.currentqDist).qAlpha.beta);
            end
        end

        function initialize_from_existing_model(obj, model)

            if model.M > obj.M
                disp('number of components is too large in previous model')
            end

            
            obj.partitionData(obj.fullData.X)
            obj.currentPartition = obj.dataTrain.partitionName;
            obj.qDistTrain.initializeVariationalDististribution;

            M_old = size(model.qDistTrain.qA.mean,2);

            obj.qDistTrain.qA.mean(:,1:M_old)=model.qDistTrain.qA.mean;
            obj.qDistTrain.qA.variance(1:M_old,1:M_old)=model.qDistTrain.qA.variance;
            obj.qDistTrain.qF.mean(1:M_old,1:M_old)=model.qDistTrain.qF.mean;
            obj.qDistTrain.qF.variance(1:M_old,1:M_old,1:M_old)=model.qDistTrain.qF.variance;

            obj.qDistTrain.qAlpha.mean(1:M_old)=model.qDistTrain.qAlpha.mean;
            obj.qDistTrain.qAlpha.variance(1:M_old)=model.qDistTrain.qAlpha.variance;
            obj.qDistTrain.qAlpha.beta(1:M_old)=model.qDistTrain.qAlpha.beta;
            obj.qDistTrain.qAlpha.alpha(1:M_old)=model.qDistTrain.qAlpha.alpha;
            
            if strcmp(obj.opts.estimationNoise,'avgShared')
                obj.qDistTrain.qSigma.mean=model.qDistTrain.qSigma.mean;
                obj.qDistTrain.qSigma.variance=model.qDistTrain.qSigma.variance;
                obj.qDistTrain.qSigma.beta=model.qDistTrain.qSigma.beta;
                obj.qDistTrain.qSigma.alpha=model.qDistTrain.qSigma.alpha;
            end

            obj.qDistTrain.initializeSufficientStatistics;
        end
        
        function [ELBO, best_model] = evaluateELBO(obj, X_test, n_init)

            old_verbose=obj.opts.verbose;
            old_showIter=obj.opts.showIter;
            old_maxIter=obj.opts.maxIter;
            old_rngInput=obj.opts.rngInput; %randi(10);
            old_activeParams = obj.opts.activeParams;

            test_data = dataClass;
            test_data.Xunfolded = X_test;

            
            test_opts=obj.opts;
            test_opts.verbose = 0;
            test_opts.maxIter = 1e4;
            test_opts.initMethod = 'mle';

            if strcmp(obj.opts.estimationNoise,'avg')
                test_opts.activeParams = {'qC','qP','qSigma'};
            elseif strcmp(obj.opts.estimationNoise,'avgShared')
                test_opts.activeParams = {'qC','qP'};
            end
            
            % NOT USED
            %test_opts.showIter = 100;
            %test_opts.debugFlag = 0;
            
            if nargin < 3
               n_init=10;
            end
            ELBO=-realmax;
            for i = 1:n_init
                test_model=varBayesModelParafac2(test_data, obj.fullData.M, test_opts);

                trainTic=tic;
                rng(i);
                test_model.initialize_from_existing_model(obj)

                test_model.computeVarDistribution;
                trainToc=toc(trainTic);

                if ELBO < test_model.qDistTrain.ELBO
                    ELBO=test_model.qDistTrain.ELBO;
                    best_model=test_model;
                    fprintf('Current best model: %d components, ELBO=%f\n', ...
                        test_model.qDistTrain.nActiveComponents('hard'),ELBO)
                end
                fprintf('%d / %d completed, took %.2f seconds\n',i,n_init,trainToc)
            end

            obj.opts.verbose=old_verbose;
            obj.opts.showIter=old_showIter;
            obj.opts.maxIter=old_maxIter;
            obj.opts.rngInput=old_rngInput; %randi(10);
            obj.opts.activeParams=old_activeParams;
        end
        
        % TODO: refactor all code below!
        
        function displayHeader(obj)
            %'ELBO','ePxz','eQz'
            names = {'# comps','# strict','ELBO Diff','ELBO','Fit'};
            %'eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n');
            fprintf('\n\n');
        end
        function displayResults(obj,qDist)
            
            %             disp(repmat('*',1,20))
            %             fprintf('ELBO: %f \n',obj.qDistTrain.ELBO);
            %             fprintf('ePxz: %f \n',obj.qDistTrain.ePxz);
            %             fprintf('eQz: %f \n',obj.qDistTrain.eQz);
            fprintf('%15.2e',...
                qDist.ELBO)%,obj.qDistTrain.ePxz,obj.qDistTrain.eQz,...
            %                   obj.qDistTrain.XqMeanLog,obj.qDistTrain.AqMeanLog,obj.qDistTrain.CqMeanLog,...
            %                   obj.qDistTrain.FqMeanLog,obj.qDistTrain.PqMeanLog,obj.qDistTrain.SigmaqMeanLog,...
            %                   obj.qDistTrain.AlphaqMeanLog,obj.qDistTrain.AqEntropy,obj.qDistTrain.CqEntropy,...
            %                   obj.qDistTrain.FqEntropy,obj.qDistTrain.PqEntropy,obj.qDistTrain.SigmaqEntropy,...
            %                   obj.qDistTrain.AlphaqEntropy)
            fprintf('\n')
            
        end
        function displayResultsAll(obj)
            if obj.opts.verbose && obj.(obj.currentData).iter ~= 0 && mod(obj.(obj.currentData).iter,obj.opts.showIter) == 0
                qDist = obj.(obj.currentqDist);
                fprintf('%5d',obj.(obj.currentData).iter);
                ELBO = qDist.data.ELBO_chain(qDist.data.iter);
                fprintf('%10d',qDist.data.n_components(qDist.data.iter),qDist.data.n_components_hard(qDist.data.iter));
                fprintf('%10.2e',...
                    obj.(obj.currentData).ELBO_diff,...
                    ELBO,qDist.data.fit_chain(qDist.data.iter));
                fprintf('\n')
            end
        end
      
        function xRecon = computeReconstructedData(obj)
            xRecon = zeros(obj.dataTrain.I,obj.dataTrain.J,obj.dataTrain.K);
            for k = 1:obj.dataTrain.K
                xRecon(:,:,k) = obj.qDistTrain.qA.mean*diag(obj.qDistTrain.qC.mean(k,:))*obj.qDistTrain.qF.mean'*obj.qDistTrain.qP.mean(:,:,k)';
            end
        end
        
        % #### Plot functions
        function plotSolutionSynthK(obj,k,MLEflag)
            
            
            %             MLEflag = 1;
            [~,index]= sort(obj.qDistTrain.qAlpha.mean,'ascend');
            
            plotParafac2SolutionK(k,bsxfun(@minus,obj.fullData.X,obj.fullData.Etrue),obj.qDistTrain.qA.mean,obj.qDistTrain.qC.mean,...
                obj.qDistTrain.qF.mean,obj.qDistTrain.qP.mean,obj.fullData.Atrue,obj.fullData.Ctrue,...
                obj.fullData.Ftrue,obj.fullData.Ptrue,MLEflag,obj.qDistTrain.nActiveComponents,index);
            
            
        end
        function plotSolutionRealMatrixK(obj,k)
            xRecon = obj.qDistTrain.qA.mean*diag(obj.qDistTrain.qC.mean(k,:))*obj.qDistTrain.qF.mean'*obj.qDistTrain.qP.mean(:,:,k)';
            
            subplot(1,3,1)
            imagesc(obj.dataTrain.X(:,:,k))
            colorbar
            subplot(1,3,2)
            imagesc(xRecon)
            colorbar
            subplot(1,3,3)
            imagesc((obj.dataTrain.X(:,:,k)-xRecon)./abs(obj.dataTrain.X(:,:,k)))
            colorbar
            
            avgError=sum(sum((obj.dataTrain.X(:,:,k)-xRecon)./abs(obj.dataTrain.X(:,:,k))))/numel(xRecon);
            fprintf('Avg. Error on relative estimate: %f\n',avgError)
        end
        function plotSolutionReal3D(obj,k,m,mask)
            set(0,'DefaultFigureWindowStyle','docked')
            %             sortedAIDX = obj.sortComponents(obj.qDistTrain.qA.mean);
            
            [~,sortedAIDX]=sort(abs(obj.qDistTrain.qC.mean(k,:)),'descend');
            
            disp(sortedAIDX(m))
            U = obj.qDistTrain.qA.mean(:,sortedAIDX(m));
            V = mean(obj.qDistTrain.qP.mean(:,sortedAIDX(m),:),3);%*obj.qDistTrain.qF.mean(sortedAIDX(m),:);
            
            plotComponent(U,V', mask, [ 53    63    46])
        end
        function plotELBO(obj,plotinterval)
            % Trim preallocated memory
            obj.ELBO_chain = nonzeros(obj.ELBO_chain)';
            obj.fit_chain = nonzeros(obj.fit_chain)';
            obj.evaltime = nonzeros(obj.evaltime)';
            obj.n_components = nonzeros(obj.n_components)';
            obj.n_components_hard = nonzeros(obj.n_components_hard)';
            
            
            subplot(3,1,1)
            semilogy(nonzeros(obj.ELBO_chain))
            
            title('ELBO')
            subplot(3,1,2)
            semilogy(plotinterval(1):(plotinterval(2)-1),diff(nonzeros(obj.ELBO_chain(plotinterval(1):plotinterval(2)))))
            xlim([0 obj.data.iter])
            title('ELBO difference between iterations')
            subplot(3,1,3)
            plot(obj.n_components_hard)
            if ~isempty(obj.data.Mtrue)
                hold on
                plot([0 numel(obj.n_components_hard)],[obj.data.Mtrue obj.data.Mtrue],'r--','LineWidth',2)
                hold off
                axis tight
            end
            ylim([0 max(obj.n_components)+2])
            ylabel('N active components')
            title('Number of active components (red=true)')
        end
        function plotHinton(obj)
            set(0,'DefaultFigureWindowStyle','normal')
            H = hinton(obj.qDistTrain.qC.mean);% ones(obj.data.K,obj.data.M));
            %             SET(H, 'INVERTHARDCOPY', 'OFF') (if printed as hardcopy, need to opposit
            %             colors)
            
        end
        
        function [fit,fit_true] = Parafac2Fit(obj,qDist,Xtrue)
            if nargin < 3
               Xtrue = [];
            end
            residual=bsxfun(@minus,qDist.data.X,obj.util.matrixProductPrSlab(...
                qDist.qA.mean,obj.util.matrixProductPrSlab(qDist.eD,...
                obj.util.matrixProductPrSlab(qDist.qF.mean',permute(...
                qDist.qP.mean,[2 1 3])))));
            
            sum_res = norm(residual(:),'fro')^2;
            sum_x = norm(qDist.data.X(:),'fro')^2;
            
            fit=gather((1-sum_res/sum_x)*100); 
            
            if ~isempty(Xtrue)
            residual_true=bsxfun(@minus,Xtrue,obj.util.matrixProductPrSlab(...
                qDist.qA.mean,obj.util.matrixProductPrSlab(qDist.eD,...
                obj.util.matrixProductPrSlab(qDist.qF.mean',permute(...
                qDist.qP.mean,[2 1 3])))));
                sum_res_true = norm(residual_true(:),'fro')^2;
                sum_x_true = norm(Xtrue(:),'fro')^2;
                fit_true=gather((1-sum_res_true/sum_x_true)*100);
                
                sum_res = norm(residual(:),'fro')^2;
                sum_x = norm(Xtrue(:),'fro')^2;
            end                       
            
              
            
        end
        
        function X_components=plotComponents(obj,xaxis)
            
            X=obj.dataTrain.X;
            K=obj.dataTrain.K;
            M=obj.dataTrain.M;
            
            A=obj.qDistTrain.qA.mean;
            C=obj.qDistTrain.qC.mean;
            F=obj.qDistTrain.qF.mean;
            P=obj.qDistTrain.qP.mean;
            
            PF = mtimesx(P,F);
            
            X_components = zeros([size(X),M]);
            for k = 1:K
                for m = 1:M
                    if any(m==[])
                        cont=-1;
                    else
                        cont=1;
                    end
                    X_components(:,:,k,m)=cont*A(:,m)*C(k,m)*PF(:,m)';
                end
            end
            
            for m = 1:M
                subplot(1,3,m)
                for k = 1:K
                    plot(xaxis,X_components(:,:,k,m));
                    %                 hold on
                end
                grid on
                axis tight
                
                if m == 2
                    xlabel('Emission wavelength')
                end
                set(gca,'fontsize',42)
            end
            
        end
        function value = SNR(obj,qDist)
            value = norm(qDist.data.Xtrue(:))^2/norm(qDist.data.Etrue(:))^2;
            disp((sum(1./qDist.data.Alphatrue))/(1/qDist.data.Sigmatrue(1)))
            disp(value)
        end
        
        
        
        function value = get.n_neg_ELBO_diff(obj)
            value = sum(diff(nonzeros(obj.ELBO_chain))<0);
        end
        function value = get.n_neg_ELBO_diff_idx(obj)
            value = find(diff(nonzeros(obj.ELBO_chain))<0);
        end
        
    end
    
    methods (Static)
        
        function sortedIdx = sortComponents(parameter)
            [~,sortedIdx] = sort(var(parameter));
        end
        
        function generatedData = generateDataFromModel(dataSettings)
            
            % Init data
            I = dataSettings.dimensions(1);
            J = dataSettings.dimensions(2);
            K = dataSettings.dimensions(3);
            M = dataSettings.dimensions(4);
            
            generatedData = dataClass;
            
            generatedData.Xunfolded = zeros([I J K]);
            generatedData.Xtrue = zeros([I J K]);
            generatedData.Mtrue= M;
            
            generatedData.Atrue = mvnrnd(zeros(generatedData.I,generatedData.Mtrue),eye(generatedData.Mtrue));
            
            if strcmp(dataSettings.initMethod,'kiers')
                F = ones(M,M)*dataSettings.congruence;
                
                for m = 1:M
                    F(m,m) = 1;
                end
                disp('Perform Cholesky factorization on F')
                generatedData.Ftrue = chol(F);
                

                if M>1
                    C = zeros(K,M);
                    C(:,1) = rand(K,1);
                    m=1;
                    n_failures=0;
                    
                    while m < M  
                        col_new=rand(K,1);
                        scores=zeros(1,m);
                        for i_comp = 1:m
                            scores(i_comp)=congruenceScore(C(:,i_comp),col_new);
                        end
                        if all(scores<0.8)
                            m=m+1;
                            C(:,m)=col_new;
                        else
                            n_failures=n_failures+1;
                        end
                        fprintf('Found %d / %d components; %d failed attempts\n',m,M,n_failures)
                    end
                else
                    C = rand(K,M);
                end
                generatedData.Ctrue = C;
                generatedData.Ptrue = zeros(generatedData.J,generatedData.Mtrue,generatedData.K);
                generatedData.Etrue = zeros(I,J,K);
                
                for k = 1:generatedData.K
                    fprintf('Processing slab %d / %d\n',k , generatedData.K)
                    generatedData.Ptrue(:,:,k) = orth(mvnrnd(zeros(generatedData.J,generatedData.Mtrue),eye(generatedData.Mtrue)));
                    generatedData.Xtrue(:,:,k) = generatedData.Atrue*diag(generatedData.Ctrue(k,:))*...
                        generatedData.Ftrue'*generatedData.Ptrue(:,:,k)';
                end
                
                ssq = generatedData.computeNoiseLevel(generatedData, dataSettings.SNR);
                
                if strcmp(dataSettings.noiseType,'heteroscedastic')
                    pp=rand(1, generatedData.K);
                elseif strcmp(dataSettings.noiseType,'homoscedastic')
                    pp=ones(1, generatedData.K);
                end
                
                %bar(pp/sum(pp))
                pp=pp/sum(pp)*generatedData.K;
                
                generatedData.Sigmatrue = ssq;
                
                for k = 1:generatedData.K
                    generatedData.Etrue(:,:,k) = mvnrnd(zeros(generatedData.I,generatedData.J)...
                        ,eye(generatedData.J)*ssq*pp(k));
                end
                
                generatedData.Xunfolded = generatedData.Xtrue+generatedData.Etrue;
                
                snr=(norm(generatedData.Xtrue(:),'fro')^2/norm(generatedData.Etrue(:),'fro')^2);
                 disp('SNR')
                 disp(snr)
                 disp('SNR (dB)')
                 disp(10*log10(snr))
%                 disp('Noise / Signal strength')
%                 disp((norm(generatedData.Etrue(:),'fro')^2/norm(generatedData.Xtrue(:),'fro')^2)*100)
%                 disp('Signal %')
%                 disp((1-norm(generatedData.Etrue(:),'fro')^2/norm(generatedData.Xunfolded(:),'fro')^2)*100)
                
            elseif strcmp(dataSettings.initMethod,'generative')
                
                if ~isfield(dataSettings,'precision')
                    SigmaPrecision = 1e12;
                    AlphaPrecision = 1e-3;
                else
                    SigmaPrecision = dataSettings.precision(1);
                    AlphaPrecision = dataSettings.precision(2);
                end
                
                generatedData.Sigmatrue = repmat(SigmaPrecision,1,generatedData.K);
                generatedData.Alphatrue = repmat(AlphaPrecision,1,generatedData.Mtrue);
                
                generatedData.Atrue = mvnrnd(zeros(generatedData.I,generatedData.Mtrue),eye(generatedData.Mtrue));
                generatedData.Ftrue = mvnrnd(zeros(generatedData.Mtrue,generatedData.Mtrue),eye(generatedData.Mtrue));
                generatedData.Ctrue = mvnrnd(zeros(generatedData.K,generatedData.Mtrue),diag(1./generatedData.Alphatrue));
                generatedData.Etrue = zeros(generatedData.I,generatedData.J,generatedData.K);
                generatedData.Ptrue = zeros(generatedData.J,generatedData.Mtrue,generatedData.K);
                
                for k = 1:generatedData.K
                    
                    generatedData.Ptrue(:,:,k) = orth(mvnrnd(zeros(generatedData.J,generatedData.Mtrue),eye(generatedData.Mtrue)));
                    
                    generatedData.Etrue(:,:,k) = mvnrnd(zeros(generatedData.I,generatedData.J)...
                        ,eye(generatedData.J)*1./generatedData.Sigmatrue(k));
                    
                    generatedData.Xtrue(:,:,k) = generatedData.Atrue*diag(generatedData.Ctrue(k,:))*...
                        generatedData.Ftrue'*generatedData.Ptrue(:,:,k)';
                    
                    generatedData.X(:,:,k) = generatedData.Xtrue(:,:,k)+generatedData.Etrue(:,:,k);
                end
            end
        end
    end
end
