classdef normalParafac2 < handle
    
    
    properties
        X
        A
        C
        F
        P
        
        M % # components
        
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
        
        
        
    end
    
end

