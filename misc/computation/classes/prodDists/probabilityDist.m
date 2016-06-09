classdef probabilityDist
    % Summary
    % Rows are independent distributions
    
    properties
        varname % Name of variable
        type % Type of distribution
        
        mean % First moment
        variance % Second moment
    end
    properties (Access=protected)
        arrayDim % Dimensions of the distribution
        I
        J
        K
    end
    
    properties (Dependent)
        entropy % entropy
    end
    
    
    methods
        function obj = probabilityDist(varname,type,arrayDim)
            % Summary of constructor
            
            
            obj.varname=varname;
            obj.type=type;
            obj.arrayDim = arrayDim;
            obj.I = arrayDim(1);
            obj.J = arrayDim(2);
            if numel(arrayDim) == 3
                obj.K = arrayDim(3);
            else
                obj.K =1;
            end
            
            if strcmp(obj.type,'multiNormal')
                obj.mean = rand(obj.arrayDim);
                
                %varIDX = obj.arrayDim([1 3:end]);
                %if numel(varIDX) > 1
                if ismatrix(obj.mean)
                    obj.variance = zeros(obj.J,obj.J,obj.I);
                    for i = 1:obj.I
                        obj.variance(:,:,i) = eye(obj.J); % cov(obj.mean);
                    end
                elseif numel(obj.arrayDim) == 3
                    obj.variance = zeros(obj.J,obj.J,obj.I,obj.K);
                    for i = 1:obj.I;
                        for k = 1:obj.K;
                            obj.variance(:,:,i,k) =  eye(obj.J); % cov(obj.mean(:,:,k));
                        end
                    end
                else
                    error('Too many dimensions');
                end
                %[obj.variance{:}]=deal(2*eye(obj.arrayDim(2)));
            else
                obj.mean = rand(obj.arrayDim);
                obj.variance = ones(obj.arrayDim);
            end
        end
        
        function value = get.entropy(obj)
            % Compute entropy on demand
            value = computeEntropy(obj);
        end
        
        
        
    end
    
    methods (Abstract, Access=protected)
        computeMean(obj) % To update Mean
        computeVariance(obj) % To update variance
        computeEntropy(obj) % To compute none constant terms of entropy
        
        %computeExpectedOuterProduct(obj)
        %computeExpectedInnerProduct(obj)
    end
end