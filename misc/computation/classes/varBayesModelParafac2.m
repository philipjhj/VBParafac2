classdef varBayesModelParafac2 < handle
    
    
    properties
        % Variational Distribution
        qDist
        
        % Settings
        verbose = 1; % 1, display, 0 hide everything
    end
    properties (Constant)
        data = dataClass;
    end
    
    methods
        function obj = varBayesModelParafac2(X,M)
            % Summary of constructor
            
            rng(2)
            if nargin < 1
                % Some dims to test
                m = 5;
                dim = 5;
                k = 5;
                obj.data.M = 1*m;
                obj.data.Mtrue=1*m;
                
                obj.data.X = zeros([dim dim k]);
            else
                obj.data.X = X;
                obj.data.M = M;
                [obj.data.I, obj.data.J, obj.data.K] = size(obj.data.X); 
            end
            
            
            % Create variational distribution object
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
            
            names = {'ELBO','ePxz','eQz','eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%15s',names{i})
            end
            fprintf('\n') ;
            end
            % Update Variational Factors until ELBO has converged
            iter = 0;
            
            diff = ELBO-ELBO_prev;
            while abs(diff)/abs(ELBO) > 1e-6
                
                iter = iter+1;
                % Update all variational factors
                obj.qDist.updateMoments;
                
                
                % Compute ELBO
                ELBO_prev = ELBO;
                ELBO = obj.qDist.ELBO;
                
                
                diff = ELBO-ELBO_prev;
                if obj.verbose
                if diff < -1e-6
                   fprintf('ELBO not converging; difference is %.4f \n',diff)
                   status=-1;
                end
                
                
                
                % Output progress
                % ...
                fprintf('%5d',iter);
                obj.displayResultsAll;
                end
            end
            %fprintf('%5d',iter);
            %obj.displayResults;
                
            if obj.verbose
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%10s',names{i})
            end
            fprintf('\n');
            end
        end
        
        function displayResults(obj)
            
%             disp(repmat('*',1,20))
%             fprintf('ELBO: %f \n',obj.qDist.ELBO);
%             fprintf('ePxz: %f \n',obj.qDist.ePxz);
%             fprintf('eQz: %f \n',obj.qDist.eQz);
              fprintf('%25.2f',...
                  obj.qDist.ELBO)%,obj.qDist.ePxz,obj.qDist.eQz,...
%                   obj.qDist.XqMeanLog,obj.qDist.AqMeanLog,obj.qDist.CqMeanLog,...
%                   obj.qDist.FqMeanLog,obj.qDist.PqMeanLog,obj.qDist.SigmaqMeanLog,...
%                   obj.qDist.AlphaqMeanLog,obj.qDist.AqEntropy,obj.qDist.CqEntropy,...
%                   obj.qDist.FqEntropy,obj.qDist.PqEntropy,obj.qDist.SigmaqEntropy,...
%                   obj.qDist.AlphaqEntropy)
               fprintf('\n')

        end
        
        function displayResultsAll(obj)
            
%             disp(repmat('*',1,20))
%             fprintf('ELBO: %f \n',obj.qDist.ELBO);
%             fprintf('ePxz: %f \n',obj.qDist.ePxz);
%             fprintf('eQz: %f \n',obj.qDist.eQz);
              fprintf('%15.2f',...
                  obj.qDist.ELBO,obj.qDist.ePxz,obj.qDist.eQz,...
                  obj.qDist.XqMeanLog,obj.qDist.AqMeanLog,obj.qDist.CqMeanLog,...
                  obj.qDist.FqMeanLog,obj.qDist.PqMeanLog,obj.qDist.SigmaqMeanLog,...
                  obj.qDist.AlphaqMeanLog,obj.qDist.AqEntropy,obj.qDist.CqEntropy,...
                  obj.qDist.FqEntropy,obj.qDist.PqEntropy,obj.qDist.SigmaqEntropy,...
                  obj.qDist.AlphaqEntropy)
               fprintf('\n')

        end
    end
end