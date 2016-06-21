classdef varBayesModelParafac2 < handle
    
    
    properties
        ELBO_chain
        % Variational Distribution
        qDist
        
        % Settings
        verbose = 1; % 1, display, 0 hide everything
        maxiter = 500;%intmax;
    end
    properties
        data = dataClass;
    end
    
    methods
        
        function obj = saveobj(obj)
            obj.data.X = [];
            disp(obj.data.X)
        end
        
        function obj = loadobj(obj)
            m=matfile('/media/data/DataAndResults/Thesis/motor_normalized_all_subs.mat');
            obj.data.X = m.Y;
        end
        
        
        function obj = varBayesModelParafac2(X,M)
            % Summary of constructor
            
%             mtimesx('SPEED');
            rng(2)
            if nargin < 1
                % Some dims to test
                m = 2;
                dim = 5;
                k = 5;
                obj.data.M = m;
                obj.data.Mtrue= m;
                
                obj.data.X = zeros([dim dim k]);
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
        
        function status = computeVarDistribution(obj)
            % Implementation of CAVI to compute the variational distribution
            % of the probabilistic Parafac2 model
            
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
            obj.data.iter = 0;
            
            diff = ELBO-ELBO_prev;
            while abs(diff)/abs(ELBO) > 1e-6 && obj.maxiter > obj.data.iter
                
                % Update all variational factors
                obj.qDist.updateMoments;
                
                
                % Compute ELBO
                ELBO_prev = ELBO;
                ELBO = obj.qDist.ELBO;
                
                if isempty(obj.ELBO_chain) || numel(obj.ELBO_chain)<obj.data.iter+1
                    obj.ELBO_chain = [obj.ELBO_chain zeros(1,100)];
                end
                obj.ELBO_chain(obj.data.iter+1) = ELBO;
                

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
        function plotSolution(obj,k,MLEflag)
            
            
%             MLEflag = 1;
           
            plotParafac2SolutionK(k,obj.data.X,obj.qDist.qA.mean,obj.qDist.qC.mean,...
                obj.qDist.qF.mean,obj.qDist.qP.mean,obj.data.Atrue,obj.data.Ctrue,...
                obj.data.Ftrue,obj.data.Ptrue,MLEflag);
            
            
        end
        
        
        
        
    end
    
    
    
end