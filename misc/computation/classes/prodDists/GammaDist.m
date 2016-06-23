classdef GammaDist < probabilityDist
    % Summary of help
    % description
    
    properties
        alpha
        beta
        MeanLog
    end
    
    
    methods
        function obj = GammaDist(varname,arrayDim)
            % Summary of constructor
            
            type = 'Gamma';
            
            if nargin < 2
                arrayDim = 1;
            end
            obj@probabilityDist(varname,type,arrayDim);
            
            obj.alpha = 2*ones(obj.arrayDim);
            obj.beta = 2*ones(obj.arrayDim);
            
        end
        
        function updateStatistics(obj)
           obj.computeMean;
           obj.computeVariance;
           obj.computeEntropy;
           obj.computeMeanLog;
        end
    end
    methods (Access = protected)
        function computeMean(obj)
            obj.mean = obj.alpha.*obj.beta;
        end
        
        function computeVariance(obj)
            obj.variance = obj.alpha.*obj.beta.^2;
        end
        
        function computeEntropy(obj)
              obj.entropy = sum(obj.alpha+log(obj.beta)+gammaln(obj.alpha)...
                  +(1-obj.alpha).*psi(obj.alpha));
        end
        
        function computeMeanLog(obj)
              obj.MeanLog = psi(obj.alpha)+log(obj.beta);
        end
        
        
    end
end