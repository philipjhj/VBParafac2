classdef gammaDist < probabilityDist
    % Summary of help
    % description
    
    properties
       alpha
       beta
    end
    
    methods
        function obj = gammaDist(varname,arrayDim)
            % Summary of constructor
            
            type = 'gamma';
            
            if nargin < 2
                arrayDim = 1;
            end
            obj@probabilityDist(varname,type,arrayDim);
            
        end
        
        function obj = computeMean(obj)
            % Test if alpha > 1
            obj.mean = obj.beta/obj.beta;
        end
        
        function obj = computeVariance(obj)
            
            obj.variance = 
            
        end
        
        function value = computeEntropy(obj)
            % Computes non constant terms of the entropy
            value = obj.aa+log(b*gamma(a))-(1+a)*psi(a);
        end
    end
end