classdef varBayesModelParafac2 < handle
    properties
        
        
        varDistribution = varDistributionC % Consist of all probDist needed for this
        ELBO
        
        % Dimensions
        I % 1st Dim
        J % 2nd Dim
        K % 3rd Dim
        M % Latent Dim
        
        
        
        
    end
    
    methods
        function computeVarDistribution(obj)
            % Implementation of CAVI to compute the variational distribution
            % of the probabilistic Parafac2 model
            
            % Initialize variational factors
            % ...
            
            % Compute ELBO
            obj.ELBO = computeELBO(obj);
            ELBO_prev = 0;
            
            % Update Variational Factors until ELBO has converged
            while obj.ELBO-ELBO_prev > 1e-6
                
                % Update all variational factors except qP
                % ...
                
                % Approximate qP's parameters
                computeVarFacP
                
                % Compute ELBO
                ELBO_prev = obj.ELBO;
                obj.ELBO = computeELBO(obj);
                
                % Output progress
                % ...
                
            end
        end
        
        function computeELBO(obj)
            % Compute the expected of log p(x,z) w.r.t q(z)
            %E_pxz = obj.varDistribution.
                %compute_ElogpX(qDistribution)+compute_ElogpA(qDistribution)+...
                %compute_ElogpC(qDistribution)+compute_ElogpF(qDistribution)+...
                %compute_ElogpP(qDistribution)+compute_ElogpSigma(qDistribution)+...
                %compute_ElogpAlpha(qDistribution);
            
            % Compute the expected of log q(z) w.r.t q(z)
            %E_qz = ...;
                
            obj.ELBO = obj.varDistribution.ePxz-obj.varDistribution.eQz;
        end
        
        function displayResults(obj)
            
        end
    end
    
end