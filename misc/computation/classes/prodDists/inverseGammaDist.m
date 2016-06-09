classdef inverseGammaDist < probabilityDist
    % Summary of help
    % description
    
    properties
        alpha
        beta
    end
    
    properties (Dependent)
        meanGamma
    end
    
    methods
        function obj = inverseGammaDist(varname,arrayDim)
            % Summary of constructor
            
            type = 'inverseGamma';
            
            if nargin < 2
                arrayDim = 1;
            end
            obj@probabilityDist(varname,type,arrayDim);
            
            obj.alpha = 2*ones(obj.arrayDim);
            obj.beta = 2*ones(obj.arrayDim);
            
        end
        
        function value = get.meanGamma(obj)
            % Test if alpha, beta > 0
            value = obj.alpha.*obj.beta;
%             value = computeMean(obj);
            %value = 1./obj.mean;
        end
        
        
%         function value = get.mean(obj)
%             % Compute entropy on demand
%             value = computeMean(obj);
%         end
    end
    methods (Access = protected)
        function value = computeMean(obj)
            % Test if alpha > 1
            value = obj.beta./(obj.alpha-1);
        end
        
        
        function obj = computeVariance(obj)
            
            %obj.variance =
            
        end
        
        function value = computeEntropy(obj)
            % Computes non constant terms of the entropy
%             value = sum(obj.alpha+log(obj.beta)+gammaln(obj.alpha)...
%                 -(1+obj.alpha).*psi(obj.alpha));
              value = sum(obj.alpha+log(obj.beta)+gammaln(obj.alpha)...
                  +(1-obj.alpha).*psi(obj.alpha));
        end
    end
end