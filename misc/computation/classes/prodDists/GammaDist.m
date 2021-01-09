classdef GammaDist < probabilityDist
    properties
        alpha
        beta
        MeanLog
    end
    
    
    methods
        function obj = GammaDist(varname,arrayDim)
            type = 'Gamma';
            
            if nargin < 2
                arrayDim = 1;
            end
            obj@probabilityDist(varname,type,arrayDim);
            
            if strcmpi(varname(2:end),'Alpha')
                obj.alpha = 1e-3*ones(obj.arrayDim);
                obj.beta = 1e3*ones(obj.arrayDim);
            elseif strcmpi(varname(2:end),'Sigma')
                obj.alpha = 1e-3*ones(obj.arrayDim);
                obj.beta = 1e3*ones(obj.arrayDim);
            end
            
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
