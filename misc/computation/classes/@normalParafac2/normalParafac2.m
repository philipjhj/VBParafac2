classdef normalParafac2 < parafac2BaseClass
    
    
    properties
        X
        A
        C
        F
        P
        
        M % # components
        
        D
        
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
            [A,F,C,P,modelFit]=parafac2(obj.X,obj.M,[0 0],[1e-12 10000 0 0 1]);
            obj.A = A;
            obj.F = F;
            obj.C = C;
            obj.P = cat(3,P{:});
            obj.fit = modelFit;
        end
        
        function value = CCDParafac2(obj)
           
            factors{1} = obj.A;
            factors{2} = obj.F;
            factors{3} = obj.C;
            
            XPk = mtimesx(obj.X,obj.P);
            
            value = corcond(XPk,factors,[],0);
        end
        
        function [fit,fit_true] = Parafac2Fit(obj,Xtrue)
            if nargin < 2
                Xtrue = [];
            end
            
            
            obj.D = zeros(size(obj.C,2),size(obj.C,2),size(obj.C,1));
            
            for k = 1:size(obj.C,1)
                obj.D(:,:,k) = diag(obj.C(k,:));
            end
            
            residual=bsxfun(@minus,obj.X,mtimesx(...
                obj.A,mtimesx(obj.D,...
                mtimesx(obj.F',permute(...
                obj.P,[2 1 3])))));
            
            
            if ~isempty(Xtrue)
                residual_true=bsxfun(@minus,Xtrue,mtimesx(...
                    obj.A,mtimesx(obj.D,...
                    mtimesx(obj.F',permute(...
                    obj.P,[2 1 3])))));
                
                
                sum_res_true = norm(residual_true(:),'fro')^2;
                sum_x_true = norm(Xtrue(:),'fro')^2;
                fit_true=(1-sum_res_true/sum_x_true)*100;
            end
            
            
            sum_res = 0;
            sum_x = 0;
            for k = 1:size(obj.C,1)
                sum_res = sum_res+norm(residual(:,:,k),'fro')^2;
                sum_x = sum_x + norm(obj.X(:,:,k),'fro')^2;
            end
            
            
            fit=(1-sum_res/sum_x)*100;
            
            
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
       
        
        function X_components=plotComponents(obj,xaxis)
            
            X=obj.X;
            K=size(X,3);
            M=obj.M;
            
            A=obj.A;
            C=obj.C;
            F=obj.F;
            P=obj.P;
            
            PF = mtimesx(P,F);
            
            X_components = zeros([size(X),M]);
            for k = 1:K
                for m = 1:M
                X_components(:,:,k,m)=A(:,m)*C(k,m)*PF(:,m)';
                end
            end
            Ms=[2 3 1];
            for m = 1:M
                subplot(1,3,m)
                for k = 1:K
                plot(xaxis,X_components(:,:,k,Ms(m)));
%                 hold on
                end
                grid on
                axis tight
                if m == 2
                xlabel('Emission wavelength')
                end
                set(gca,'fontsize',42)
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
end

