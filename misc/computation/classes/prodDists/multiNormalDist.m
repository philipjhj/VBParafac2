classdef multiNormalDist < probabilityDist
    % Summary of help
    % description
    
    properties
        VarEqual
        util
    end
    
    properties (Dependent)
        meanOuterProduct
        meanOuterProductSingle
        meanInnerProductMatrix % E[X^T.X] (Constant value, if correct)
        %         meanInnerProductTransposed % E[X.X^T] (Constant value, if correct)
        %meanInnerProductSumComponent
        distSquared
    end
    
    methods
        function obj = multiNormalDist(varname,arrayDim,VarEqualBool,util)
            % Summary of constructor
            
            type = 'multiNormal';
            
            if nargin < 3
                VarEqualBool = false;
            end
            
            if nargin < 2
                arrayDim = 1;
            end
            
            obj@probabilityDist(varname,type,arrayDim,VarEqualBool);
            
            obj.util = util;
            
            obj.VarEqual = VarEqualBool;
            
        end
        
        function value = get.distSquared(obj)
            value = obj.computeElementWiseSquared;
        end
        
        function value = meanInnerProductSumComponent(obj,A)
            
            if ismatrix(obj.mean)
                value = obj.computeMeanInnerProduct(eye(obj.J));
               
                %                 %                   value = trace(obj.mean*obj.mean')+trace(sum(obj.variance,3));
                %                 if obj.VarEqual
                %                     value = sum(sum(obj.mean.^2))+obj.I*trace(obj.variance);
                %                 else
                %                     value = sum(sum(obj.mean.^2))+trace(sum(obj.variance,3));
                %                 end
            else
                value = sum(sum(sum(obj.mean.^2)))+obj.I*...
                    sum(sum(sum(sum(bsxfun(@times,eye(obj.J),obj.variance)))));
            end
        end
        
        function value = get.meanOuterProduct(obj)
            % Computes the mean of the outer product (x'*x) for all vectors
            % x in distribution
            
            value = obj.computeMeanOuterProduct;
        end
        
        function value = get.meanOuterProductSingle(obj)
            % Computes the mean of the outer product (x'*x) for all vectors
            % x in distribution
            
            value = obj.computeMeanOuterProductSingle;
        end
        
        
        function value = get.meanInnerProductMatrix(obj)
            % Computes the inner product (x*x') for all vectors x in
            % distribuiton
            value = obj.computeMeanInnerProductScaled;
        end
        
        function updateStatistics(obj)
            obj.computeEntropy;
        end
    end
    
    % Functions with input from outside class
    methods (Access = public)
        function matrixIPS = computeMeanInnerProductScaledSlabs(obj,A,diagFlag)
            % Computes the inner product <x'*A*x> scaled by A for all
            % vectors x or only diag of <x'*A*x> if diagFlag is set
            if nargin < 2
                A = eye(obj.J);
            end
            
            if nargin < 3
                diagFlag = 0;
            end
            
            if diagFlag
                a = squeeze(sum(bsxfun(@times,obj.util.matrixProductPrSlab(obj.mean,A),obj.mean),2));
                b = squeeze(sum(sum(bsxfun(@times,eye(obj.J),obj.util.matrixProductPrSlab(A,obj.variance)),1),2))';
            else
                a = squeeze(obj.util.matrixProductPrSlab(obj.util.matrixProductPrSlab(obj.mean,A),obj.mean'));
                
                x = squeeze(sum(sum(bsxfun(@times,eye(obj.J),obj.util.matrixProductPrSlab(A,obj.variance)),1),2));
                b = obj.util.matrixDiagonalPrSlab(x);
                
            end
            
            matrixIPS = bsxfun(@plus,a,b);
            
        end
        
    end
    
    % Functions not requiring input from outside the class
    methods (Access = protected)
        function value = computeElementWiseSquared(obj)
            % E[X.^2]
            variance=permute(squeeze(sum(bsxfun(@times,eye(obj.J),obj.variance),2)),[2 1 3]);
            
            value=bsxfun(@plus,variance,obj.mean.^2);
        end
        
        function value = computeMeanInnerProduct(obj,A)
            % E[xAx'] = Tr(A*Sigma)+mu*A*mu' =
            
            if obj.VarEqual
                value = obj.I*sum(sum(A'.*obj.variance))+sum(sum((obj.mean'*obj.mean).*A));
            else
                value = sum(sum(A'.*sum(obj.variance,3)))+sum(sum((obj.mean'*obj.mean).*A));
            end
        end
        
        
        function value = computeMeanOuterProduct(obj)
            % E[x'x] = Sigma + mu'*mu
            
            if ismatrix(obj.mean)
                % sum of Sigma+mu*mu' for all mean vectors and Sigmas
                
                if obj.VarEqual
                    value = obj.I*obj.variance + obj.mean'*obj.mean;
                else
                    value = sum(obj.variance,3) + obj.mean'*obj.mean;
                end
            end
            
        end
        
        function value = computeMeanOuterProductSingle(obj)
            if ismatrix(obj.mean)
                % Sigma + mu*mu' for each mean vector and sigma pair
                
                mu = reshape(obj.mean',1,obj.arrayDim(2),obj.arrayDim(1));
                mumu = obj.util.matrixProductPrSlab(permute(mu,[2 1 3]),...
                    mu);
                value = bsxfun(@plus,mumu,obj.variance);
                
            end
            
        end
        
        function computeMean(obj)
            obj.mean;
        end
        
        function computeVariance(obj)
            obj.variance;
        end
        
        function computeEntropy(obj)
            % Computes non constant part of the entropy
            
            if ndims(obj.mean) == 3
                K = obj.K;
            else
                K = 1;
            end
            
            gpuflag = 0;
            
            if isa(obj.variance,'gpuArray')
                %[~, U, P] = pagefun(@lu,obj.variance);
                %du = pagefun(@diag,U);
                %c = pagefun(det,P) * pagefun(@prod,pagefun(@sign,du));
                %v = pagefun(@log,c) + pagefun(@sum,pagefun(@log,pagefun(@abs,du)));
                
                %logdetValue = gather(v);
                %logdetValue = sum(sum(logdetValue,4),3);
                
                %value = obj.VarDim*(obj.J/2*(1+log(2*pi)))...
                %           +1/2*logdetValue;
                obj.variance = gather(obj.variance);
                gpuflag = 1;
            end
            
            %         else
            value = 0;
            for i=1:obj.VarDim
                for k=1:K
                    value = value+obj.J/2*(1+log(2*pi))...
                        +1/2*logdet((obj.variance(:,:,i,k)));
                end
            end
            if gpuflag
                obj.variance = gpuArray(obj.variance);
            end
            
            %        end
            
            
            if obj.VarEqual
                obj.entropy = obj.I*value;
            else
                obj.entropy = value;
            end
        end
    end
end
