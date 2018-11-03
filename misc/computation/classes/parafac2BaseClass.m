classdef parafac2BaseClass < handle
    
    properties (Abstract = true)
        X
        A
        C
        F
        P
        
        M % # components
    end
    
    methods
        function plot_profiles(obj, include_scale, normalize, attempt_flip)
            if nargin < 4
               attempt_flip=false; 
            end
            if nargin < 3
                normalize=false;
            end
            if nargin < 2
                include_scale=false;
            end
            
            if normalize
                [~,C_local,~,~,profiles] = obj.normalize();
            else
                profiles=mtimesx(obj.F',permute(obj.P,[2 1 3]));
                C_local=obj.C;
            end
            
            if include_scale
                profiles = bsxfun(@times,permute(C_local,[2 3 1]),profiles);
            end
            
            if attempt_flip
                for m=1:obj.M
                    profiles_m=squeeze(profiles(m,:,:));
                    [~,idx_flip1]=max(abs(profiles_m));
                    for Nidx = 1:numel(idx_flip1)
                        profiles(m,:,Nidx)=profiles_m(:,Nidx)*sign(profiles_m(idx_flip1(Nidx),Nidx));
                    end
                end
            end
            
            colors=linspecer(obj.M,'qualitative');
            
            for m = 1:obj.M
                hold on
                plot(squeeze(profiles(m,:,:)),'color',colors(m,:))
            end
            hold off
        end
    
        function [A, C, F, P, FPk] = normalize(obj)
            [A,C,F,P,FPk] = normalizeParafac2(obj.A, obj.C,obj.F,obj.P);
        end
        
    end
    
end