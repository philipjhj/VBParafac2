classdef normalParafac2 < handle
    
    
    properties
        X
        A
        C
        D
        F
        P
        
        M % # components
        
        Xtrue_m
        Xrecon_m
        Xrecon
        fit
    end
    
    
    methods
        function obj = normalParafac2(X,A,C,F,P)
            
            obj.X = X; %permute(X,[2 1 3]);
            
            if nargin > 1
                
                obj.A = A;
                obj.C = C;
                obj.F = F;
                obj.P = P;
                
                
                if iscell(obj.P)
                    obj.P = cat(3,obj.P{:});
                end
            end
        end
        
        function obj = fitParafac2(obj,M)
            
            obj.M = M;
            [A,F,C,P,modelFit]=parafac2(obj.X,obj.M,[0 0],[0 0 0 0 1]);
            
            
            obj.A = A;
            obj.F = F;
            obj.C = C;
            obj.P = cat(3,P{:});
            obj.fit = modelFit;
            
            
        end
        
        function [fit,fit_true] = Parafac2Fit(obj,Xtrue)
            
            obj.D = zeros(size(obj.C,2),size(obj.C,2),size(obj.C,1));
            
            for k = 1:size(obj.C,1)
                obj.D(:,:,k) = diag(obj.C(k,:));
            end
            
            residual=bsxfun(@minus,obj.X,mtimesx(...
                obj.A,mtimesx(obj.D,...
                mtimesx(obj.F',permute(...
                obj.P,[2 1 3])))));
            
            sum_res = 0;
            sum_x = 0;
            for k = 1:size(obj.C,1)
                sum_res = sum_res+norm(residual(:,:,k),'fro')^2;
                sum_x = sum_x + norm(obj.X(:,:,k),'fro')^2;
            end
            
            fit=(1-sum_res/sum_x)*100;
            fit_true=(1-sum_res/norm(Xtrue(:))^2)*100;
            
        end
        
        function [fit] = SSE_Fit(obj,Xtrue)
            
            obj.D = zeros(size(obj.C,2),size(obj.C,2),size(obj.C,1));
            
            for k = 1:size(obj.C,1)
                obj.D(:,:,k) = diag(obj.C(k,:));
            end
            
            residual=bsxfun(@minus,Xtrue,mtimesx(...
                obj.A,mtimesx(obj.D,...
                mtimesx(obj.F',permute(...
                obj.P,[2 1 3])))));
            
            
            fit=norm(residual(:),'fro')^2/norm(Xtrue(:),'fro')^2;
            
            
        end
       
        
        function obj = compute_reconstruction(obj,data)
            
            obj.Xrecon_m = zeros(size(obj.X,1),size(obj.X,2),size(obj.X,3),size(obj.C,2));
            
            
            obj.D = zeros(size(obj.C,2),size(obj.C,2),size(obj.C,1));
            
            for k = 1:size(obj.C,1)
                obj.D(:,:,k) = diag(obj.C(k,:));
            end
            
            
            A = obj.A;
            D = obj.D;
            F = obj.F;
            P = obj.P;
            
            for m = 1:size(A,2)
                obj.Xrecon_m(:,:,:,m) = mtimesx(...
                    mtimesx(mtimesx(...
                    A(:,m),D(m,m,:)),F(:,m)'),permute(P,[2 1 3]));
            end
            
            obj.Xrecon = sum(obj.Xrecon_m,4);
            
            if nargin>1
            if ~isempty(data.Xtrue) && isempty(obj.Xtrue_m) 
                
                obj.Xtrue_m = zeros(size(obj.X,1),size(obj.X,2),size(obj.X,3),size(obj.C,1));
                
                A = data.Atrue;
                D = bsxfun(@mtimes,reshape(data.Ctrue',1,...
                data.Mtrue,data.K),...
                repmat(eye(data.Mtrue),1,1,data.K));
                F = data.Ftrue;
                P = data.Ptrue;
                
                for m = 1:data.Mtrue
                    obj.Xtrue_m(:,:,:,m) = mtimesx(...
                        mtimesx(mtimesx(...
                        A(:,m),D(m,m,:)),F(:,m)'),permute(P,[2 1 3]));
                end
                
                
            end
            
            end
        end
        
        
    end
        
        
%         function nActive = nActiveComponents(obj,method)
%             
%             if nargin<2
%                method = obj.opts.nActiveComponents; 
%             end
%             
%             if strcmp(method,'hard')
%                 if isa(obj.qC.mean,'gpuArray')
%                     nActive = sum(sum(gather(obj.qC.mean))~=0);
%                 else
%                     nActive = sum(sum(obj.qC.mean)~=0);
%                 end
%             elseif strcmp(method,'threshold')
%                 nActive = find(cumsum(sort(1./obj.qAlpha.mean,'descend')/sum(1./obj.qAlpha.mean))>0.95,1);
%             end
%         end
        
    
end
