classdef varBayesModelParafac2 < handle
    
    
    properties
        ELBO_chain
        n_components
        evaltime
        % Variational Distribution
        qDist
        
        data
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
        
        function obj = saveobj(obj)
            if obj.data.I > 1e4
                obj.data.X = [];
                obj.qDist.eAiDFtPtPFDAi = [];
                %             else
                %                 obj = obj;
            end
        end
        
        function obj = varBayesModelParafac2(X,M)
            % Summary of constructor
            
            %             mtimesx('SPEED');
            
            obj.data = dataClass;
            obj.opts = optionsClass;
            obj.util = utilitiesClass(obj);
            
            obj.opts.maxTime = realmax;
            
            % Some dims to test
            I = 100;
            J = I;
            K = 10;
            
            Mtrue = 10;
            Mesti = Mtrue;
            
            
            
            if  nargin < 1
                % Use default precision values
                obj.data = obj.generateDataFromModel([I J K Mtrue]);
                obj.data.M = Mesti;
            elseif numel(X) == 2
                % Set precision with input in X
                obj.data = obj.generateDataFromModel([I J K Mtrue],X);
                obj.data.M = Mesti;
            elseif isa(X,'dataClass')
                % Store/reference already generated data
                obj.data = X;
                obj.data.M = M;
            else
                % Load given data set
                obj.data.X = X;
                obj.data.M = M;
                [obj.data.I, obj.data.J, obj.data.K] = size(obj.data.X);
            end
            
            
            % Create variational distribution object
            obj.qDist = varDistributionC(obj);
        end
        
        
        function restartqDist(obj)
            clear obj.qDist;
            obj.qDist = varDistributionC(obj);
            obj.data.iter = 1;
            
            obj.ELBO_chain = [];
            obj.n_components = [];
            obj.evaltime = [];
            
        end
        
        function [stopReason,errorStatus] = computeVarDistribution(obj,maxiter)
            
            if obj.opts.verbose
                disp('Starting CAVI with:')
                disp(obj.qDist.opts.activeParams)
                %             disp('\n')
            end
            
            % Implementation of CAVI to compute the variational distribution
            % of the probabilistic Parafac2 model
            
            % Set options and meta data
            if nargin < 2
                obj.opts.maxiter = intmax;
            elseif isempty(obj.opts.maxiter) || obj.opts.maxiter < maxiter
                obj.opts.maxiter = maxiter;
            end
            
            stopReason = 0;
            errorStatus = 0;
            
            if isempty(obj.data.iter)
                obj.data.iter = 1;
            end
            
            
            % Initialize on first iteration
            if obj.data.iter == 1;
                rng(obj.opts.rngInput);
                obj.qDist = obj.qDist.initDist;
                
                obj.n_components(1) = obj.qDist.nActiveComponents;
                obj.evaltime(1) = 1e-10;
            end
            
            startTime = obj.evaltime(obj.evaltime==max(obj.evaltime));
            
            % Compute Initial ELBO
            ELBO = obj.qDist.ELBO;
            ELBO_prev = 0;
            diff = ELBO-ELBO_prev;
            
            if obj.opts.verbose
                obj.displayHeader;
            end
            
            % Update Variational Factors until ELBO has converged
            ticCAVI=tic;
            while abs(diff)/abs(ELBO) > obj.opts.tol && obj.opts.maxiter+1 > obj.data.iter...
                    && obj.opts.maxTime > obj.evaltime(obj.evaltime==max(obj.evaltime))
                
                % Update active (see options) variational factors
                obj.qDist.updateMoments;
                
                
                % Compute ELBO
                ELBO_prev = ELBO;
                ELBO = obj.qDist.ELBO;
                
                diff = ELBO-ELBO_prev;
                
                % Store progress information
                if isempty(obj.ELBO_chain) || numel(obj.ELBO_chain)<obj.data.iter+1
                    obj.ELBO_chain = [obj.ELBO_chain zeros(1,100)];
                    obj.evaltime = [obj.evaltime zeros(1,100)];
                    obj.n_components = [obj.n_components zeros(1,100)];
                end
                obj.ELBO_chain(obj.data.iter) = ELBO;
                obj.evaltime(obj.data.iter) = toc(ticCAVI)+startTime;
                obj.n_components(obj.data.iter) = obj.qDist.nActiveComponents;
                
                % Check convergence
                if obj.opts.debugFlag >= 1 && diff/abs(ELBO) < -1e-7 && obj.data.iter>0
                    warning('off','backtrace')
                    warning('At iter %d ELBO not converging; relativ diff. is %.10f, diff; %.4f \n',...
                        obj.data.iter,diff/abs(ELBO),diff)
                    warning('on','backtrace')
                    %                     keyboard
                    errorStatus=-1;
                    fprintf('%5d',obj.data.iter);
                    obj.displayResultsAll(diff);
                    %                     break
                end
                
                
                
                % Display Progress
                if obj.opts.verbose && obj.data.iter ~= 0 && mod(obj.data.iter,obj.opts.showIter) == 0
                    
                    % Output progress
                    % ...
                    fprintf('%5d',obj.data.iter);
                    obj.displayResultsAll(diff);
                end
                
                if obj.Parafac2Fit/100>=1%(1-obj.opts.tol)
                    break;
                end
                
                obj.data.iter = obj.data.iter+1;
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
                elseif obj.opts.maxiter <= obj.data.iter
                    fprintf('CAVI has stopped at iteration %d (max iteration) with change %f\n',obj.data.iter-1,abs(diff)/abs(ELBO))
                    stopReason = 2;
                elseif obj.opts.maxTime <= obj.evaltime(obj.evaltime==max(obj.evaltime))
                    fprintf('CAVI has stopped after %f s. evaluation (max time) with change %f\n',obj.evaltime(obj.evaltime==max(obj.evaltime)),abs(diff)/abs(ELBO))
                    stopReason = 3;
                elseif obj.Parafac2Fit/100>(1-obj.opts.tol)
                    fprintf('CAVI has stopped with a fit of %f %% with change %f\n',obj.Parafac2Fit,abs(diff)/abs(ELBO))
                    stopReason = 4;    
                end
            end
            
            % Trim preallocated memory
            obj.ELBO_chain = nonzeros(obj.ELBO_chain)';
            obj.evaltime = nonzeros(obj.evaltime)';
            obj.n_components = nonzeros(obj.n_components)';
        end
        
        
        function displayHeader(obj)
            %'ELBO','ePxz','eQz'
            names = {'# comps','ELBO Diff','ELBO','Fit'};
                %'eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n');
        end
        
        
        function displayResults(obj)
            
            %             disp(repmat('*',1,20))
            %             fprintf('ELBO: %f \n',obj.qDist.ELBO);
            %             fprintf('ePxz: %f \n',obj.qDist.ePxz);
            %             fprintf('eQz: %f \n',obj.qDist.eQz);
            fprintf('%15.2e',...
                obj.qDist.ELBO)%,obj.qDist.ePxz,obj.qDist.eQz,...
            %                   obj.qDist.XqMeanLog,obj.qDist.AqMeanLog,obj.qDist.CqMeanLog,...
            %                   obj.qDist.FqMeanLog,obj.qDist.PqMeanLog,obj.qDist.SigmaqMeanLog,...
            %                   obj.qDist.AlphaqMeanLog,obj.qDist.AqEntropy,obj.qDist.CqEntropy,...
            %                   obj.qDist.FqEntropy,obj.qDist.PqEntropy,obj.qDist.SigmaqEntropy,...
            %                   obj.qDist.AlphaqEntropy)
            fprintf('\n')
            
        end
        
        function displayResultsAll(obj,diff)
            
            %             disp(repmat('*',1,20))
            %             fprintf('ELBO: %f \n',obj.qDist.ELBO);
            %             fprintf('ePxz: %f \n',obj.qDist.ePxz);
            %             fprintf('eQz: %f \n',obj.qDist.eQz);
            %                   obj.qDist.ELBO,obj.qDist.ePxz,obj.qDist.eQz,...
            ELBO = obj.ELBO_chain(obj.data.iter);
            fprintf('%10d',obj.n_components(obj.data.iter));
            fprintf('%10.2e',...
                diff,...
                ELBO,obj.Parafac2Fit);
%                 obj.qDist.qXMeanLog,obj.qDist.qAMeanLog,obj.qDist.qCMeanLog,...
%                 obj.qDist.qFMeanLog,obj.qDist.qPMeanLog,obj.qDist.qSigmaMeanLog,...
%                 obj.qDist.qAlphaMeanLog,obj.qDist.qAEntropy,obj.qDist.qCEntropy,...
%                 obj.qDist.qFEntropy,obj.qDist.qPEntropy,obj.qDist.qSigmaEntropy,...
%                 obj.qDist.qAlphaEntropy)
            fprintf('\n')
            
        end
        
        
        
        % #### Plot functions
        function plotSolutionSynthK(obj,k,MLEflag)
            
            
            %             MLEflag = 1;
            
            [~,index]= sort(obj.qDist.qAlpha.mean,'ascend');
            
            plotParafac2SolutionK(k,bsxfun(@minus,obj.data.X,obj.data.Etrue),obj.qDist.qA.mean,obj.qDist.qC.mean,...
                obj.qDist.qF.mean,obj.qDist.qP.mean,obj.data.Atrue,obj.data.Ctrue,...
                obj.data.Ftrue,obj.data.Ptrue,MLEflag,obj.qDist.nActiveComponents,index);
            
            
        end
        
        
        function plotSolutionRealMatrixK(obj,k)
            xRecon = obj.qDist.qA.mean*diag(obj.qDist.qC.mean(k,:))*obj.qDist.qF.mean'*obj.qDist.qP.mean(:,:,k)';
            
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
            %             sortedAIDX = obj.sortComponents(obj.qDist.qA.mean);
            
            [~,sortedAIDX]=sort(abs(obj.qDist.qC.mean(k,:)),'descend');
            
            disp(sortedAIDX(m))
            U = obj.qDist.qA.mean(:,sortedAIDX(m));
            V = mean(obj.qDist.qP.mean(:,sortedAIDX(m),:),3);%*obj.qDist.qF.mean(sortedAIDX(m),:);
            
            plotComponent(U,V', mask, [ 53    63    46])
        end
        %
        
        function plotELBO(obj,plotinterval)
            % Trim preallocated memory
            obj.ELBO_chain = nonzeros(obj.ELBO_chain)';
            obj.evaltime = nonzeros(obj.evaltime)';
            obj.n_components = nonzeros(obj.n_components)';
            
            
            subplot(3,1,1)
            semilogy(nonzeros(obj.ELBO_chain))
            title('ELBO')
            subplot(3,1,2)
            semilogy(plotinterval(1):(plotinterval(2)-1),diff(nonzeros(obj.ELBO_chain(plotinterval(1):plotinterval(2)))))
            xlim([0 obj.data.iter])
            title('ELBO difference between iterations')
            subplot(3,1,3)
            plot(obj.n_components)
            if ~isempty(obj.data.Mtrue)
                hold on
                plot([0 numel(obj.n_components)],[obj.data.Mtrue obj.data.Mtrue],'r','LineWidth',2)
                hold off
                axis tight
            end
            ylim([0 max(obj.n_components)+2])
            ylabel('N active components')
            title('Number of active components (red=true)')
        end
        
        function plotHinton(obj)
            set(0,'DefaultFigureWindowStyle','normal')
            H = hinton(obj.qDist.qC.mean);% ones(obj.data.K,obj.data.M));
            %             SET(H, 'INVERTHARDCOPY', 'OFF') (if printed as hardcopy, need to opposit
            %             colors)
            
        end
        
        function fit = Parafac2Fit(obj)
            
            residual=bsxfun(@minus,obj.data.X,obj.util.matrixProductPrSlab(...
                obj.qDist.qA.mean,obj.util.matrixProductPrSlab(obj.qDist.eD,...
                obj.util.matrixProductPrSlab(obj.qDist.qF.mean',permute(...
                obj.qDist.qP.mean,[2 1 3])))));
            
            sum_res = 0;
            sum_x = 0;
            for k = 1:obj.data.K
                sum_res = sum_res+norm(residual(:,:,k))^2;
                sum_x = sum_x + norm(obj.data.X(:,:,k))^2;
            end
            
            fit=(1-sum_res/sum_x)*100;
            
        end
        
        
        function SNR(obj)
            disp(norm(obj.data.Xtrue(:))^2/norm(obj.data.Etrue(:))^2)
            disp((sum(1./obj.data.Alphatrue))/(1/obj.data.Sigmatrue(1)))
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
        %                 obj.qDist.compute_eAiDFtPtPFDAi;
        %             end
        %end
        
        
        function generatedData = generateDataFromModel(dimensions,precision)
            %             rng('default')
            %             rng('shuffle')
            
            if nargin < 2
%                 SigmaPrecision = 1e4;
%                 AlphaPrecision = 1e-8;
            else
                SigmaPrecision = precision(1);
                AlphaPrecision = precision(2);
            end
            
            I = dimensions(1);
            J = dimensions(2);
            K = dimensions(3);
            M = dimensions(4);
            
            
            generatedData = dataClass;
            
            generatedData.X = zeros([I J K]);
            generatedData.Xtrue = zeros([I J K]);
            generatedData.Mtrue= M;
            
            %             generatedData.SigmaAtrue = 1e-1;
            %             generatedData.SigmaBtrue = 1;
            %             generatedData.AlphaAtrue = 1e2;
            %             generatedData.AlphaBtrue = 1;
            
            %             generatedData.Sigmatrue = gamrnd(generatedData.SigmaAtrue,generatedData.SigmaBtrue,1,generatedData.K);
            %             generatedData.Alphatrue = gamrnd(generatedData.AlphaAtrue,generatedData.AlphaBtrue,1,generatedData.Mtrue);
            
            generatedData.Sigmatrue = repmat(SigmaPrecision,1,generatedData.K);
            generatedData.Alphatrue = repmat(AlphaPrecision,1,generatedData.Mtrue);
            
            generatedData.Atrue = mvnrnd(zeros(generatedData.I,generatedData.Mtrue),eye(generatedData.Mtrue));
            generatedData.Ftrue = mvnrnd(zeros(generatedData.Mtrue,generatedData.Mtrue),eye(generatedData.Mtrue));
            
            generatedData.Ctrue = mvnrnd(zeros(generatedData.K,generatedData.Mtrue),diag(1./generatedData.Alphatrue));
            
            generatedData.Etrue = zeros(generatedData.I,generatedData.J,generatedData.K);
            generatedData.X = zeros(generatedData.I,generatedData.J,generatedData.K);
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
