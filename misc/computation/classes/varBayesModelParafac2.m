classdef varBayesModelParafac2 < handle
    
    properties
        ELBO_chain
        fit_chain
        n_components
        n_components_hard
        evaltime
        evaltime_cpu
        
        test_ELBO_chain
        test_fit_chain
        test_n_components
        test_n_components_hard
        test_evaltime
        test_evaltime_cpu
        
        % Variational Distributions
        qDistTrain
        qDistTest
        
        testData
        trainData
        
        util
        
        opts
        % TODO: add statistics class
    end
    
    properties (Dependent)
        n_neg_ELBO_diff
        n_neg_ELBO_diff_idx
    end
    
    methods
        
        function value = get.n_neg_ELBO_diff(obj)
            value = sum(diff(nonzeros(obj.ELBO_chain))<0);
        end
        
        function value = get.n_neg_ELBO_diff_idx(obj)
            value = find(diff(nonzeros(obj.ELBO_chain))<0);
        end
        
        function partitionData(obj,X)
            disp('Partitioning data slices into 80% training set, 20% test set');
            
            testPct = 0;
            
            obj.trainData = dataClass;
            obj.testData = dataClass;
            
            dims = size(X);
            K=dims(end);
            nPrSlice = prod(dims(1:(end-1)));
            
            allSlices = 1:K;
            testSlices = randsample(K,floor(testPct*K));
            trainSlices = allSlices(~ismember(allSlices,testSlices));
            
            X = reshape(X,nPrSlice,K);
            obj.testData.Xunfolded = reshape(X(:,testSlices),[dims(1:(end-1)) numel(testSlices)]);
            obj.trainData.Xunfolded = reshape(X(:,trainSlices),[dims(1:(end-1)) numel(trainSlices)]);
        end
        
        function obj = varBayesModelParafac2(X,M)
            obj.opts = optionsClass;
            obj.util = utilitiesClass(obj);
            
            obj.opts.maxTime = realmax;
            
            
            % Load given data set
            obj.partitionData(X);
            
            obj.trainData.M = M;
            obj.testData.M = M;
            
            % Create variational distribution object
            obj.qDistTrain = varDistributionC(obj,obj.trainData);
            obj.qDistTest = varDistributionC(obj,obj.testData);
        end
        
        
        function restartqDistTrain(obj)
            clear obj.qDistTrain;
            obj.qDistTrain = varDistributionC(obj);
            obj.trainData.iter = 1;
            
            obj.ELBO_chain = [];
            obj.fit_chain = [];
            obj.n_components = [];
            obj.n_components_hard = [];
            obj.evaltime = [];
            obj.evaltime_cpu = [];
            
        end
        
        function fitTrainingData(obj,maxiter)
            [obj.qDistTrain,obj.trainData]=obj.computeVarDistribution(maxiter,obj.qDistTrain,obj.trainData);
        end
        
        function fitTestData(obj,maxiter)
            previous_opt = obj.opts.activeParams;
            obj.opts.activeParams = {'qP','qC','qSigma'};
            obj.qDistTest.qA = obj.qDistTrain.qA;
            obj.qDistTest.qF = obj.qDistTrain.qF;
            [obj.qDistTest,obj.testData]=obj.computeVarDistribution(maxiter,obj.qDistTest,obj.testData);
            obj.opts.activeParams = previous_opt;
        end
        
        function [qDist,data,stopReason,errorStatus] = computeVarDistribution(obj,maxiter,qDist,data)
            
            if obj.opts.verbose
                disp('Starting CAVI with:')
                disp(qDist.opts.activeParams)
                %             disp('\n')
            end
            
            % Set options and meta data
            if nargin < 2
                obj.opts.maxiter = intmax;
            elseif isempty(obj.opts.maxiter) || obj.opts.maxiter < maxiter
                obj.opts.maxiter = maxiter;
            end
            
            stopReason = 0;
            errorStatus = 0;
            
            if isempty(data.iter)
                data.iter = 1;
            end
            
            
            % Initialize on first iteration
            if data.iter == 1
                rng(obj.opts.rngInput);
                qDist = qDist.initializeVariationalDististribution;
                
                obj.n_components(1) = qDist.nActiveComponents;
                obj.n_components_hard(1) = qDist.nActiveComponents('hard');
                obj.evaltime(1) = 1e-10;
                obj.evaltime_cpu(1) = 1e-10;
            end
            
            startTime = obj.evaltime(obj.evaltime==max(obj.evaltime));
            startTime_cpu = obj.evaltime_cpu(obj.evaltime_cpu==max(obj.evaltime_cpu));
            
            % Compute Initial ELBO
            ELBO = qDist.ELBO;
            ELBO_prev = 0;
            diff = ELBO-ELBO_prev;
            
            if obj.opts.verbose
                obj.displayHeader;
            end
            
            % Update Variational Factors until ELBO has converged
            ticCAVI=tic;
            t_CAVI_cpu = cputime;
            while abs(diff)/abs(ELBO) > obj.opts.tol && obj.opts.maxiter+1 > data.iter...
                    && obj.opts.maxTime > obj.evaltime(obj.evaltime==max(obj.evaltime))
                
                % Update active (see options) variational factors
                qDist.updateMoments;
                
                
                % Compute ELBO
                ELBO_prev = ELBO;
                ELBO = qDist.ELBO;
                
                diff = ELBO-ELBO_prev;
                
                % Store progress information
                if isempty(obj.ELBO_chain) || numel(obj.ELBO_chain)<data.iter+1
                    obj.ELBO_chain = [obj.ELBO_chain zeros(1,100)];
                    obj.fit_chain = [obj.fit_chain zeros(1,100)];
                    obj.evaltime = [obj.evaltime zeros(1,100)];
                    obj.evaltime_cpu = [obj.evaltime_cpu zeros(1,100)];
                    obj.n_components = [obj.n_components zeros(1,100)];
                    obj.n_components_hard = [obj.n_components_hard zeros(1,100)];
                end
                obj.ELBO_chain(data.iter) = ELBO;
                obj.fit_chain(data.iter) = obj.Parafac2Fit(qDist);
                obj.evaltime(data.iter) = toc(ticCAVI)+startTime;
                obj.evaltime_cpu(data.iter) = cputime-t_CAVI_cpu+startTime_cpu;
                obj.n_components(data.iter) = qDist.nActiveComponents;
                obj.n_components_hard(data.iter) = qDist.nActiveComponents('hard');
                
                % Check convergence
                if obj.opts.debugFlag >= 1 && diff/abs(ELBO) < -1e-7 && data.iter>0
                    warning('off','backtrace')
                    warning('At iter %d ELBO not converging; relativ diff. is %.10f, diff; %.4f \n',...
                        data.iter,diff/abs(ELBO),diff)
                    warning('on','backtrace')
                    %                     keyboard
                    errorStatus=-1;
                    fprintf('%5d',data.iter);
                    obj.displayResultsAll(diff,qDist);
                    %                     break
                end
                
                
                
                % Display Progress
                if obj.opts.verbose && data.iter ~= 0 && mod(data.iter,obj.opts.showIter) == 0
                    
                    % Output progress
                    % ...
                    fprintf('%5d',data.iter);
                    obj.displayResultsAll(diff,qDist);
                end
                
                if obj.Parafac2Fit(qDist)/100>=1%(1-obj.opts.tol)
                    break;
                end
                
                data.iter = data.iter+1;
            end
            
            if obj.opts.verbose
                obj.displayHeader;
                fprintf('\n\n');
                
            end
            
            % Stop Criteria Message
            if obj.opts.verbose
                if diff/abs(ELBO) < obj.opts.tol
                    fprintf('CAVI has converged with last change %f\n', abs(diff)/abs(ELBO))
                    stopReason = 1;
                elseif obj.opts.maxiter <= data.iter
                    fprintf('CAVI has stopped at iteration %d (max iteration) with change %f\n',data.iter-1,abs(diff)/abs(ELBO))
                    stopReason = 2;
                elseif obj.opts.maxTime <= obj.evaltime(obj.evaltime==max(obj.evaltime))
                    fprintf('CAVI has stopped after %f s. evaluation (max time) with change %f\n',obj.evaltime(obj.evaltime==max(obj.evaltime)),abs(diff)/abs(ELBO))
                    stopReason = 3;
                elseif obj.Parafac2Fit(qDist)/100>(1-obj.opts.tol)
                    fprintf('CAVI has stopped with a fit of %f %% with change %f\n',obj.Parafac2Fit(qDist),abs(diff)/abs(ELBO))
                    stopReason = 4;
                end
                
            end
            
            % Trim preallocated memory
            obj.ELBO_chain = nonzeros(obj.ELBO_chain)';
            obj.fit_chain = nonzeros(obj.fit_chain)';
            obj.evaltime = nonzeros(obj.evaltime)';
            obj.evaltime_cpu = nonzeros(obj.evaltime_cpu)';
            obj.n_components = nonzeros(obj.n_components)';
            obj.n_components_hard = nonzeros(obj.n_components_hard)';
            
            % Gather results to CPU
            if strcmpi(obj.opts.matrixProductPrSlab,'gpu')
                qDist.qA.mean = gather(qDist.qA.mean);
                qDist.qA.variance = gather(qDist.qA.variance);
                qDist.qC.mean = gather(qDist.qC.mean);
                qDist.qC.variance = gather(qDist.qC.variance);
                qDist.qF.mean = gather(qDist.qF.mean);
                qDist.qF.variance = gather(qDist.qF.variance);
                qDist.qP.mean = gather(qDist.qP.mean);
                qDist.qP.variance = gather(qDist.qP.variance);

                qDist.qSigma.alpha = gather(qDist.qSigma.alpha);
                qDist.qSigma.beta= gather(qDist.qSigma.beta);
                qDist.qAlpha.alpha = gather(qDist.qAlpha.alpha);
                qDist.qAlpha.beta= gather(qDist.qAlpha.beta);
            end
        end
        
        
        function displayHeader(obj)
            %'ELBO','ePxz','eQz'
            names = {'# comps','# strict','ELBO Diff','ELBO','Fit'};
            %'eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n');
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
        
        function displayResultsAll(obj,diff,qDist)
            
            %             disp(repmat('*',1,20))
            %             fprintf('ELBO: %f \n',obj.qDistTrain.ELBO);
            %             fprintf('ePxz: %f \n',obj.qDistTrain.ePxz);
            %             fprintf('eQz: %f \n',obj.qDistTrain.eQz);
            %                   obj.qDistTrain.ELBO,obj.qDistTrain.ePxz,obj.qDistTrain.eQz,...
            ELBO = obj.ELBO_chain(qDist.data.iter);
            fprintf('%10d',obj.n_components(qDist.data.iter),obj.n_components_hard(qDist.data.iter));
            fprintf('%10.2e',...
                diff,...
                ELBO,obj.fit_chain(qDist.data.iter));
            %                 obj.qDistTrain.qXMeanLog,obj.qDistTrain.qAMeanLog,obj.qDistTrain.qCMeanLog,...
            %                 obj.qDistTrain.qFMeanLog,obj.qDistTrain.qPMeanLog,obj.qDistTrain.qSigmaMeanLog,...
            %                 obj.qDistTrain.qAlphaMeanLog,obj.qDistTrain.qAEntropy,obj.qDistTrain.qCEntropy,...
            %                 obj.qDistTrain.qFEntropy,obj.qDistTrain.qPEntropy,obj.qDistTrain.qSigmaEntropy,...
            %                 obj.qDistTrain.qAlphaEntropy)
            fprintf('\n')
            
        end
        
        
        
        function obj = compute_reconstruction(obj)
            
            obj.data.Xrecon_m = zeros(obj.data.I,obj.data.J,obj.data.K,obj.data.M);
            
            A = obj.qDistTrain.qA.mean;
            D = obj.qDistTrain.eD;
            F = obj.qDistTrain.qF.mean;
            P = obj.qDistTrain.qP.mean;
            
            for m = 1:obj.data.M
                obj.data.Xrecon_m(:,:,:,m) = obj.util.matrixProductPrSlab(...
                    obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(...
                    A(:,m),D(m,m,:)),F(:,m)'),permute(P,[2 1 3]));
            end
            
            obj.data.Xrecon = sum(obj.data.Xrecon_m,4);
            
            if ~isempty(obj.data.Xtrue) && isempty(obj.data.Xtrue_m)
                
                obj.data.Xtrue_m = zeros(obj.data.I,obj.data.J,obj.data.K,obj.data.M);
                
                A = obj.data.Atrue;
                D = bsxfun(@mtimes,reshape(obj.data.Ctrue',1,...
                    obj.data.Mtrue,obj.data.K),...
                    repmat(eye(obj.data.Mtrue),1,1,obj.data.K));
                F = obj.data.Ftrue;
                P = obj.data.Ptrue;
                
                for m = 1:obj.data.Mtrue
                    obj.data.Xtrue_m(:,:,:,m) = obj.util.matrixProductPrSlab(...
                        obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(...
                        A(:,m),D(m,m,:)),F(:,m)'),permute(P,[2 1 3]));
                end
                
                
            end
            
        end
        
        
        
        % #### Plot functions
        function plotSolutionSynthK(obj,k,MLEflag)
            
            
            %             MLEflag = 1;
            
            [~,index]= sort(obj.qDistTrain.qAlpha.mean,'ascend');
            
            plotParafac2SolutionK(k,bsxfun(@minus,obj.data.X,obj.data.Etrue),obj.qDistTrain.qA.mean,obj.qDistTrain.qC.mean,...
                obj.qDistTrain.qF.mean,obj.qDistTrain.qP.mean,obj.data.Atrue,obj.data.Ctrue,...
                obj.data.Ftrue,obj.data.Ptrue,MLEflag,obj.qDistTrain.nActiveComponents,index);
            
            
        end
        
        
        function plotSolutionRealMatrixK(obj,k)
            xRecon = obj.qDistTrain.qA.mean*diag(obj.qDistTrain.qC.mean(k,:))*obj.qDistTrain.qF.mean'*obj.qDistTrain.qP.mean(:,:,k)';
            
            subplot(1,3,1)
            imagesc(obj.data.X(:,:,k))
            colorbar
            subplot(1,3,2)
            imagesc(xRecon)
            colorbar
            subplot(1,3,3)
            imagesc((obj.data.X(:,:,k)-xRecon)./abs(obj.data.X(:,:,k)))
            colorbar
            
            avgError=sum(sum((obj.data.X(:,:,k)-xRecon)./abs(obj.data.X(:,:,k))))/numel(xRecon);
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
        %
        
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
        
        function [fit,fit_true] = Parafac2Fit(obj,qDist)
            
            residual=bsxfun(@minus,qDist.data.X,obj.util.matrixProductPrSlab(...
                qDist.qA.mean,obj.util.matrixProductPrSlab(qDist.eD,...
                obj.util.matrixProductPrSlab(qDist.qF.mean',permute(...
                qDist.qP.mean,[2 1 3])))));
            
            sum_res = norm(residual(:),'fro')^2;
            sum_x = norm(qDist.data.X(:),'fro')^2;
            
            fit=gather((1-sum_res/sum_x)*100);
            fit_true=gather((1-sum_res/norm(qDist.data.Xtrue(:))^2)*100);
            
        end
        
        
        function value = SNR(obj,qDist)
            value = norm(qDist.data.Xtrue(:))^2/norm(qDist.data.Etrue(:))^2;
            disp((sum(1./qDist.data.Alphatrue))/(1/qDist.data.Sigmatrue(1)))
            disp(value)
        end
        
    end
    
    methods (Static)
        
        function sortedIdx = sortComponents(parameter)
            [~,sortedIdx] = sort(var(parameter));
        end
        
        %function obj = loadobj(obj)
        %             if obj.data.I>1e4
        %m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
        %obj.data.X = m.Y;
        %                 obj.qDistTrain.compute_eAiDFtPtPFDAi;
        %             end
        %end
        
        
        function generatedData = generateDataFromModel(options)
            
            % Init data
            I = options.dimensions(1);
            J = options.dimensions(2);
            K = options.dimensions(3);
            M = options.dimensions(4);
            
            generatedData = dataClass;
            
            generatedData.Xunfolded = zeros([I J K]);
            generatedData.Xtrue = zeros([I J K]);
            generatedData.Mtrue= M;
            
            generatedData.Atrue = mvnrnd(zeros(generatedData.I,generatedData.Mtrue),eye(generatedData.Mtrue));
            
            if strcmp(options.initMethod,'kiers')
                F = ones(M,M)*options.congruence;
                
                for m = 1:M
                    F(m,m) = 1;
                end
                
                generatedData.Ftrue = chol(F);
                
                score = zeros(M);
                score(1,2) = 1;
                i=0;
                while any(any(nonzeros(score)>0.8))
                    C = rand(K,M);
                    for m1 = 1:M
                        for m2 = (m1+1):M
                            score(m1,m2) = congruenceScore(C(:,m1),C(:,m2));
                        end
                    end
                    i = i+1;
                end
                
                generatedData.Ctrue = 30*C;
                
                generatedData.Ptrue = zeros(generatedData.J,generatedData.Mtrue,generatedData.K);
                
                generatedData.Etrue = zeros(I,J,K);
                
                
                
                
                
                
                
                for k = 1:generatedData.K
                    
                    generatedData.Ptrue(:,:,k) = orth(mvnrnd(zeros(generatedData.J,generatedData.Mtrue),eye(generatedData.Mtrue)));
                    
                    
                    
                    generatedData.Xtrue(:,:,k) = generatedData.Atrue*diag(generatedData.Ctrue(k,:))*...
                        generatedData.Ftrue'*generatedData.Ptrue(:,:,k)';
                    
                    
                end
                
                %                 SNR = linspace(-4,0,generatedData.K);
                SNR = [repmat(-20,1,5) zeros(1,generatedData.K-5)];
                
                
                ssq = generatedData.computeNoiseLevel(generatedData,SNR);
                
                generatedData.Sigmatrue = ssq;
                
                for k = 1:generatedData.K
                    
                    %                     generatedData.Etrue(:,:,k) = mvnrnd(zeros(generatedData.I,generatedData.J)...
                    %                         ,eye(generatedData.J)*ssq);%(k));
                    %
                    generatedData.Etrue(:,:,k) = mvnrnd(zeros(generatedData.I,generatedData.J)...
                        ,eye(generatedData.J)*ssq(k));
                    
                end
                
                
                generatedData.Xunfolded = generatedData.Xtrue+generatedData.Etrue;
                
            elseif strcmp(options.initMethod,'generative')
                
                if ~isfield(options,'precision')
                    SigmaPrecision = 1e12;
                    AlphaPrecision = 1e-3;
                else
                    SigmaPrecision = options.precision(1);
                    AlphaPrecision = options.precision(2);
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


