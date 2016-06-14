classdef probabilityDist < handle
    % Summary
    % Rows are independent distributions
    
    properties
        varname % Name of variable
        type    % Type of distribution
        
        mean
        variance
        entropy
    end
    properties (Access=protected)
        arrayDim % Dimensions of the distribution
        I
        J
        K
        VarDim
    end
    
    
    methods
        function obj = probabilityDist(varname,type,arrayDim,VarEqualBool)
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
                
                if VarEqualBool
                    obj.VarDim = 1;
                else
                    obj.VarDim = obj.I;
                end
                
                if ismatrix(obj.mean)
                    obj.variance = zeros(obj.J,obj.J,obj.VarDim);
                    for i = 1:obj.VarDim
                        obj.variance(:,:,i) = eye(obj.J); % cov(obj.mean);
                    end
                elseif numel(obj.arrayDim) == 3
                    obj.variance = zeros(obj.J,obj.J,obj.VarDim,obj.K);
                    for i = 1:obj.VarDim;
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
    end
    
    methods (Abstract, Access=protected)
        computeMean(obj) % To update Mean
        computeVariance(obj) % To update variance
        computeEntropy(obj) % To compute none constant terms of entropy
    end
    methods (Abstract)
        updateStatistics(obj)
    end
end