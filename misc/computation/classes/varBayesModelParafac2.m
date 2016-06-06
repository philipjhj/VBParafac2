classdef varBayesModelParafac2 < handle
    
    
    properties
        % Variational Distribution
        qDist
    end
    properties (Constant)
        data = dataClass;
    end
    
    methods
        function obj = varBayesModelParafac2(X,M)
            % Summary of constructor
            
            
            if nargin < 1
                % Some dims to test
                obj.data.M = 4;
                obj.data.Mtrue=2;
                obj.data.X = zeros([5 5 2]);
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
            names = {'ELBO','ePxz','eQz','eX','eA','eC','eF','eP','eSigma','eAlpha','hA','hC','hF','hP','hSigma','hAlpha'};
            
            fprintf('%5s','Iter');
            for i = 1:numel(names)
                fprintf('%15s',names{i})
            end
            fprintf('\n') ;
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
                if diff < -1e-6
                   fprintf('ELBO not converging; difference is %.4f \n',diff)
                   status=-1;
                end
                
                
                
                % Output progress
                % ...
                fprintf('%5d',iter);
                obj.displayResultsAll;
                
            end
            %fprintf('%5d',iter);
            %obj.displayResults;
                
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