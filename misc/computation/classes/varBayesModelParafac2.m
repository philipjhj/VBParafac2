classdef varBayesModelParafac2 < handle
    
    
    properties
        ELBO_chain
        evaltime
        % Variational Distribution
        qDist
        
        % Settings
        verbose = 1; % 1, display, 0 hide everything
        maxiter
    end
    properties
        data = dataClass;
    end
    
    methods
        
        function obj = saveobj(obj)
            obj.data.X = [];
            obj.qDist.eAiDFtPtPFDAi = [];
        end
        
        function obj = varBayesModelParafac2(X,M)
            % Summary of constructor
            
%             mtimesx('SPEED');
            
            if nargin < 1
                % Some dims to test
                I = 5;
                J = I;
                K = 5;
                
                Mtrue = 2;
                Mesti = 2*Mtrue;
                
                obj.data = obj.generateDataFromModel([I J K Mtrue]);
                obj.data.M = Mesti;
                
            elseif isa(X,'dataClass')
                obj.data = X;
                obj.data.M = M;
            else
                obj.data.X = X;
                obj.data.M = M;
                [obj.data.I, obj.data.J, obj.data.K] = size(obj.data.X); 
            end
            
            
            % Create variational distribution object
            obj.qDist = varDistributionC(obj);
            
        end
        
        
        function restartqDist(obj)
            
            obj.qDist = varDistributionC(obj);
            
            
        end
        
        function status = computeVarDistribution(obj,maxiter)
            % Implementation of CAVI to compute the variational distribution
            % of the probabilistic Parafac2 model
            
            if nargin < 2
                obj.maxiter = intmax;
            else
                obj.maxiter = maxiter;
            end
            
            status = 0;
            
            % Compute ELBO
            ELBO = obj.qDist.ELBO;
            ELBO_prev = 0;
            
            if obj.verbose
            %'ELBO','ePxz','eQz'
            names = {'ELBO Diff','eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n') ;
            end
            % Update Variational Factors until ELBO has converged
            if isempty(obj.data.iter)
                obj.data.iter = 0;
%             else
%                 obj.data.iter = obj.data.iter+1;
            end
            
            diff = ELBO-ELBO_prev;
            tic;
            while abs(diff)/abs(ELBO) > 1e-7 && obj.maxiter > obj.data.iter
                
                % Update all variational factors
                obj.qDist.updateMoments;
                
                
                % Compute ELBO
                ELBO_prev = ELBO;
                ELBO = obj.qDist.ELBO;
                
                if isempty(obj.ELBO_chain) || numel(obj.ELBO_chain)<obj.data.iter+1
                    obj.ELBO_chain = [obj.ELBO_chain zeros(1,100)];
                    obj.evaltime = [obj.evaltime zeros(1,100)];
                end
                obj.ELBO_chain(obj.data.iter+1) = ELBO;
                obj.evaltime(obj.data.iter+1) = toc;

                diff = ELBO-ELBO_prev;
                if obj.verbose && obj.data.iter ~= 0
                if diff < -1e-6 
                   fprintf('ELBO not converging; difference is %.4f \n',diff)
                   status=-1;
                end
                
                
                % Output progress
                % ...
                fprintf('%5d',obj.data.iter);
                obj.displayResultsAll(diff);
                end
                obj.data.iter = obj.data.iter+1;
            end
            %fprintf('%5d',obj.iter);
            %obj.displayResults;
                
            if obj.verbose
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n');
            end
            obj.ELBO_chain = nonzeros(obj.ELBO_chain)';
            obj.evaltime = nonzeros(obj.evaltime)';
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
            fprintf('%10.2e',...
                diff,obj.qDist.ELBO)%...
%             obj.qDist.qXMeanLog,obj.qDist.qAMeanLog,obj.qDist.qCMeanLog,...
%             obj.qDist.qFMeanLog,obj.qDist.qPMeanLog,obj.qDist.qSigmaMeanLog,...
%             obj.qDist.qAlphaMeanLog,obj.qDist.qAEntropy,obj.qDist.qCEntropy,...
%             obj.qDist.qFEntropy,obj.qDist.qPEntropy,obj.qDist.qSigmaEntropy,...
%             obj.qDist.qAlphaEntropy)
        fprintf('\n')
        
        end
        
        
        
        % #### Plot functions
        function plotSolutionSynthK(obj,k,MLEflag)
            
            
%             MLEflag = 1;
           
            plotParafac2SolutionK(k,obj.data.X,obj.qDist.qA.mean,obj.qDist.qC.mean,...
                obj.qDist.qF.mean,obj.qDist.qP.mean,obj.data.Atrue,obj.data.Ctrue,...
                obj.data.Ftrue,obj.data.Ptrue,MLEflag);
            
            
        end
        
        
        function plotSolutionMatrixRealK(obj,k)
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
            U = obj.qDist.qA.mean*diag(obj.qDist.qC.mean(k,:));
            V = obj.qDist.qP.mean(:,:,k)*obj.qDist.qF.mean;

            plotComponent(U(:,m),V(:,m)', mask, [ 53    63    46])
        end
        
    end
    
    methods (Static)
        
        function obj = loadobj(obj)
            m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
            obj.data.X = m.Y;
            obj.qDist.compute_eAiDFtPtPFDAi;
        end
        
        
        function generatedData = generateDataFromModel(dimensions)
            rng('default')
            rng(2)
            
            I = dimensions(1);
            J = dimensions(2);
            K = dimensions(3);
            M = dimensions(4);
            
            
            generatedData = dataClass;
            
            generatedData.X = zeros([I J K]);
            generatedData.Mtrue= M;
            
            generatedData.SigmaAtrue = 1;
            generatedData.SigmaBtrue = 1e4;
            generatedData.AlphaAtrue = 1;
            generatedData.AlphaBtrue = 1e-4;
            
%             generatedData.Sigmatrue = repmat(1e1,1,generatedData.K);%gamrnd(generatedData.SigmaAtrue,generatedData.SigmaBtrue,1,generatedData.K);
%             generatedData.Alphatrue = repmat(1e-4,1,generatedData.Mtrue);%gamrnd(generatedData.AlphaAtrue,generatedData.AlphaBtrue,1,generatedData.Mtrue);
            
            generatedData.Sigmatrue = gamrnd(generatedData.SigmaAtrue,generatedData.SigmaBtrue,1,generatedData.K);
            generatedData.Alphatrue = gamrnd(generatedData.AlphaAtrue,generatedData.AlphaBtrue,1,generatedData.Mtrue);
            
            
            generatedData.Atrue = 10*mvnrnd(zeros(generatedData.I,generatedData.Mtrue),eye(generatedData.Mtrue));
            generatedData.Ftrue = 10*mvnrnd(zeros(generatedData.Mtrue,generatedData.Mtrue),eye(generatedData.Mtrue));
            
            generatedData.Ctrue = 10*mvnrnd(zeros(generatedData.K,generatedData.Mtrue),sqrt(1./generatedData.Alphatrue)*eye(generatedData.Mtrue));
            
            generatedData.Etrue = zeros(generatedData.I,generatedData.J,generatedData.K);
            generatedData.X = zeros(generatedData.I,generatedData.J,generatedData.K);
            generatedData.Ptrue = zeros(generatedData.J,generatedData.Mtrue,generatedData.K);
            
            for k = 1:generatedData.K
                
                generatedData.Ptrue(:,:,k) = orth(mvnrnd(zeros(generatedData.J,generatedData.Mtrue),eye(generatedData.Mtrue)));
                
                generatedData.Etrue(:,:,k) = mvnrnd(zeros(generatedData.I,generatedData.J)...
                    ,eye(generatedData.J)*sqrt(1./generatedData.Sigmatrue(k)));
                generatedData.X(:,:,k) = generatedData.Atrue*diag(generatedData.Ctrue(k,:))*...
                    generatedData.Ftrue'*generatedData.Ptrue(:,:,k)'+generatedData.Etrue(:,:,k);
            end
            
        end
        
        
    end
    
    
    
end