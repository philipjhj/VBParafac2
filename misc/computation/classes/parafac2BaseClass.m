classdef parafac2BaseClass < handle
    
    properties (Abstract = true)
        X
        A
        C
        F
        P
        
        M % # components
        
    end
    
    properties (Access=protected)
       profiles_values = nan
    end
    
    methods
        
        function profiles = compute_profiles(obj, include_scale, normalize, attempt_flip)
            if isnan(obj.profiles_values)
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
                obj.profiles_values=profiles;
                
            else
                profiles=obj.profiles_values;
            end
        end
        
        function plts=plot_profiles(obj, include_scale, normalize, attempt_flip, colors, subplots)
            
            if nargin < 6
                subplots=false;
            end
            
            profiles=obj.compute_profiles(include_scale, normalize, attempt_flip);
            
            Ms=[];
            for i = unique(colors,'rows')'
                Ms=[Ms find(ismember(colors,i','rows'))']; 
            end
%             disp(Ms)
%             disp(obj.M)


            
            n_M=length(Ms);
            if subplots
                for m = fliplr(Ms)
                        subplot(n_M,1,m)
                        plts(:,m)=plot(squeeze(profiles(m,:,:)),'color',colors(m,:),'LineWidth',1.7);
                end                
            else
                for m = fliplr(Ms)
    %                 disp(m)
                    plts(:,m)=plot(squeeze(profiles(m,:,:)),'color',colors(m,:),'LineWidth',1.7);


                    hold on
                end
                hold off
            end
        end
    
        function [A, C, F, P, FPk] = normalize(obj)
            [A,C,F,P,FPk] = normalizeParafac2(obj.A, obj.C,obj.F,obj.P);
        end
        
        function X_reconstructed = compute_reconstruction(obj)
            [A_norm,C_norm,~,~,FPk_norm] = obj.normalize();
            CFPk = bsxfun(@times,permute(C_norm, [2 3 1]),FPk_norm);
            
            X_reconstructed=zeros(size(A_norm,1)*size(CFPk,2)*size(CFPk,3),obj.M);
            for t=1:obj.M
                Xtemp=mtimesx(A_norm(:,t),CFPk(t,:,:));
                X_reconstructed(:,t) = Xtemp(:);
            end
        end
        
    end
    
end