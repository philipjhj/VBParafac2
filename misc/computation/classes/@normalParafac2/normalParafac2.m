classdef normalParafac2 < handle
    
    
    properties
        X
        A
        C
        F
        P
    end
    
    
    methods
        function obj = normalParafac2(X,A,C,F,P)
            obj.X = permute(X,[2 1 3]);
            obj.A = A;
            obj.C = C;
            obj.F = F;
            obj.P = P;
            
            
            if iscell(obj.P)
                obj.P = cat(3,obj.P{:});
            end
        end
        
        
    end
    
end

