classdef uniNormalDist < probabilityDist
    % Summary of help
    % description
    properties (Dependent)
        meanDistSquared % X^2
        covMatrix
    end
    
    methods
        function obj = uniNormalDist(varname,arrayDim)
            % Summary of constructor
            
            type = 'UniNormal';
            
            if nargin < 2
                arrayDim = 1;
            end
            obj@probabilityDist(varname,type,arrayDim);
        end
        
        
        function value = get.meanDistSquared(obj)
            value = obj.variance + obj.mean.^2;
        end
        
        function covMatrix = get.covMatrix(obj)
            covMatrix=obj.computeCovMatrix;
        end
        
    end
    methods (Access=protected)
        
        function value = computeMean(obj)
            value = obj.mean;
        end
        
        function covMatrix=computeCovMatrix(obj)
            K = obj.arrayDim(1);
            %covMatrix=cell(1,K);
            covMatrix = zeros(obj.arrayDim(2),obj.arrayDim(2),K);
            for k = 1:K
                covMatrix(:,:,k) = obj.mean(k,:)'*obj.mean(k,:) +...
                    diag(obj.variance(k,:));
            end
        end
        
        
        function computeVariance(obj)
            
        end
        
        function value = computeEntropy(obj)
            % Computes non constant terms of the entropy
            value = sum(sum(1/2*log(2*pi*exp(1)*obj.variance)));
        end
    end
end