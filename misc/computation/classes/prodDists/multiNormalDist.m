classdef multiNormalDist < probabilityDist
    % Summary of help
    % description
    
    properties
        VarEqual
    
    end
    
    properties (Dependent)
        meanOuterProduct
        meanOuterProductSingle
        meanInnerProductMatrix % E[X^T.X] (Constant value, if correct)
%         meanInnerProductTransposed % E[X.X^T] (Constant value, if correct)
        meanInnerProductSumComponent
        distSquared
    end
    
    methods
        function obj = multiNormalDist(varname,arrayDim,VarEqualBool)
            % Summary of constructor
            
            type = 'multiNormal';
            
            if nargin < 3
                VarEqualBool = false;
            end
            
            if nargin < 2
                arrayDim = 1;
            end
            
            obj@probabilityDist(varname,type,arrayDim,VarEqualBool);
            
            obj.VarEqual = VarEqualBool;
            
        end
        
        function value = get.distSquared(obj)
           value = obj.computeElementWiseSquared; 
        end
        
        function value = get.meanInnerProductSumComponent(obj)
            %                 value = sum(sum(sum(obj.distSquared)));
            
            if ismatrix(obj.mean)
                %                   value = trace(obj.mean*obj.mean')+trace(sum(obj.variance,3));
                if obj.VarEqual
                    value = sum(sum(obj.mean.^2))+obj.I*trace(obj.variance);
                else
                    value = sum(sum(obj.mean.^2))+trace(sum(obj.variance,3));
                end
            else
                value = 0;
                for k = 1:obj.arrayDim(3)
                    value = value+sum(sum(obj.mean(:,:,k).^2))+obj.I*trace(obj.variance(:,:,:,k));
                    %                          value = value+trace(obj.mean(:,:,k)*obj.mean(:,:,k)')+trace(sum(obj.variance(:,:,:,k),3));
                end
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
            
            if ismatrix(obj.mean)
                value = obj.computeMeanInnerProductScaled;
            else
                K=obj.arrayDim(3);
                value = zeros(obj.arrayDim(1));
                for k = 1:K
                    value(:,:,k) = obj.computeMeanInnerProductScaled(k);
                end
            end
        end
        
%         function value = get.meanInnerProductTransposed(obj)
%             
%             % Covariance is square and symmetric
%             
%             if ismatrix(obj.mean)
%                 meanT = obj.mean';
%             else
% %                 obj.mean = permute(obj.mean,[2 1 3]);
%                 K=obj.arrayDim(3);
%                 value = zeros([size(obj.mean,2) size(obj.mean,2) size(obj.mean,3)]);
%                 for k = 1:K
%                     value(:,:,k) = obj.mean(:,:,k)'*obj.mean(:,:,k);
%                     value(:,:,k) = value(:,:,k) +...
%                     diag(diag(sum(obj.variance(:,:,:,k),3)));
%                 end
%             end
%         end
        
        function updateStatistics(obj)
           obj.computeEntropy;
        end
    end
    
    % Functions with input from outside class
    methods (Access = public)
        function matrixIPS = computeMeanInnerProductScaledSlabs(obj,A,diagFlag)
            % Computes the inner product <x'*A*x> scaled by A for all 
            % vectors x or only diag of <x'*A*x> if diagFlag is set
            
            if nargin < 3
               diagFlag = 0; 
            end
            
            I = obj.I;
            
            if numel(obj.arrayDim) == 2
                K_dist = 1;
            else
                K_dist = obj.arrayDim(3);
            end
            
            if diagFlag
                matrixIPSsize = [I 1];
            else
                matrixIPSsize = [I I];
            end
                
            checkDimA = ndims(A);
            if checkDimA == 3 % Init for per slab in A
                K_A = size(A,3);
                matrixIPS = zeros([matrixIPSsize K_A]);
                %[matrixIPS{:}] = deal(zeros(I,J));%,K_dist));
            else % Init for one slab A
                K_A = 1;
                matrixIPS = zeros([matrixIPSsize K_dist]);
            end
            
            for k = 1:max(K_dist,K_A)
                % Make the k'th slab in distribution match with k'th slab
                % in A
                
                % Compute the inner product for all vectors in the
                % distribution matrix
                % either for same A for all slabs in the distribution matrix
                % or k'th A for the k'th slab of the distribution
                if checkDimA
                    matrixIPS(:,:,k) = obj.computeMeanInnerProductScaled(1,A(:,:,k),diagFlag);
                else
                    matrixIPS(:,:,k) = obj.computeMeanInnerProductScaled(k,A,diagFlag);
                end
            end
            if diagFlag
                matrixIPS = squeeze(matrixIPS);
            end
        end
        
    end
    
    % Functions not requiring input from outside the class
    methods (Access = protected)
        function value = computeElementWiseSquared(obj)
            
            value = zeros(obj.I,obj.J,obj.K);
            
            for k = 1:obj.K;
                variance=zeros(obj.I,obj.J);
                for i = 1:obj.I
                    variance(i,:) = diag(obj.variance(:,:,i,k))';
                end
                value(:,:,k) = variance+obj.mean(:,:,k).^2;
            end
        end
        
        function value = computeMeanOuterProduct(obj)
           
            if ismatrix(obj.mean)
                % 2D
                
                if obj.VarEqual
                    value = obj.I*obj.variance + obj.mean'*obj.mean;
                else
                    value = sum(obj.variance,3) + obj.mean'*obj.mean;
                end
            else
                % 3D
                K = obj.arrayDim(3);
                value = zeros(obj.arrayDim(2),obj.arrayDim(2),K);
                for k = 1:K
                    if obj.varEqual
                        value(:,:,k) = obj.I*obj.variance(:,:,:,k) + ...
                            obj.mean(:,:,k)'*obj.mean(:,:,k);    
                    else
                        value(:,:,k) = sum(obj.variance(:,:,:,k),3) + ...
                            obj.mean(:,:,k)'*obj.mean(:,:,k);
                    end
                end
            end 
            
        end
        
        function value = computeMeanOuterProductSingle(obj)
            
            value = zeros(obj.arrayDim(2),obj.arrayDim(2),obj.arrayDim(1));
            
            if ismatrix(obj.mean)
                % 2D
                for i = 1:obj.arrayDim(1);
                    value(:,:,i) = obj.variance(:,:,i) + ...
                        obj.mean(i,:)'*obj.mean(i,:);
                end
            else
                % 3D
                K = obj.arrayDim(3);
                value = zeros(obj.arrayDim(2),obj.arrayDim(2),K);
                for k = 1:K
                    value(:,:,k) = sum(obj.variance(:,:,:,k),3) + ...
                        obj.mean(:,:,k)'*obj.mean(:,:,k);
                end
            end
            
        end
        
        function meanIPS = computeMeanInnerProductScaled(obj,k_dist,A,diagFlag)
            % Computes <x'*A*y> for x == y
            % Computes <x'>*A*<y> for x ~= y
            % If A is not given, the identity matrix is applied
            
            if nargin < 2
                A = eye(size(obj.mean,2));
                k_dist=1;
            elseif nargin < 3
                A = eye(size(obj.mean,2));
            end
            
            if ~diagFlag
            meanIPS=obj.mean(:,:,k_dist)*A*obj.mean(:,:,k_dist)';
            
            % Add variance term to diagonal
            if obj.VarEqual
                meanIPS = meanIPS+...
                    diag(ones(1,obj.I)*trace(A*obj.variance(:,:,1,k_dist)));
            else
                for i=1:size(meanIPS,1)
                    covariance = obj.variance(:,:,i,k_dist);
                    meanIPS(i,i) = meanIPS(i,i)+trace(A*covariance);
                end
            end
            else
%                 meanIPS = zeros(obj.I,1);
%                 for i=1:obj.I
%                 meanIPS(i) = obj.mean(i,:,k_dist)*A*obj.mean(i,:,k_dist)';
%                 end
%                 
                meanIPS = sum(obj.mean(:,:,k_dist)*A.*obj.mean(:,:,k_dist),2);
                
                if obj.VarDim == 1
                    meanIPS = bsxfun(@plus,meanIPS,trace(A*obj.variance(:,:,k_dist)));
                    
                else
                    for i=1:size(meanIPS,1)
                        covariance = obj.variance(:,:,i,k_dist);
                        meanIPS(i) = meanIPS(i)+trace(A*covariance);
                    end
                    
                end
                
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
            
            value = 0;
            for i=1:obj.VarDim
                for k=1:K
                    value = value+obj.arrayDim(2)/2*(1+log(2*pi))...
                        +1/2*logdet((obj.variance(:,:,i,k)));
%                     value = value+log(sqrt(det(2*pi*obj.variance(:,:,i,k))))...
%                         +obj.J/2;
                end
            end
            
            if obj.VarEqual
                obj.entropy = obj.I*value;
            else
                obj.entropy = value;
            end
        end
    end
end