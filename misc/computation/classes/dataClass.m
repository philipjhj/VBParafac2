classdef dataClass < handle
    properties
        % Data
        X
        M
        
        % True latent variables (if generated)
        Atrue
        Ctrue
        Ftrue
        Ptrue
        Etrue
        Sigmatrue
        SigmaAtrue
        SigmaBtrue
        Alphatrue
        AlphaAtrue
        AlphaBtrue
        Mtrue % True number of components
    end
    
    properties
        % Dims
        I
        J
        K
    end
    methods
        function set.X(obj,value)
            obj.X = value;
            obj.I = size(obj.X,1);
            obj.J = size(obj.X,2);
            obj.K = size(obj.X,3);
        end
        function set.M(obj,value)
            % Implment error checking here
            obj.M = value;
        end
    end
end