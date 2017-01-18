classdef dataClass < handle
    properties
        % Data
        X
        M
        
        % True latent variables (if generated)
        Xtrue = [] %without noise, X above contains noise
        Xtrue_m = [] % X split into components
        Atrue = []
        Ctrue = []
        Ftrue = []
        Ptrue = []
        Etrue = []
        Sigmatrue = []
        SigmaAtrue = []
        SigmaBtrue = []
        Alphatrue = []
        AlphaAtrue = []
        AlphaBtrue = []
        Mtrue = [] % True number of components
        
        iter
        ELBOall
        
        errorIters = []
        errorIters_parameter = {};
        
        % Reconstruction
        Xrecon
        Xrecon_m
        
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
    
    
    methods (Static)
        function ssq = computeNoiseLevel(obj,SNR)
        
           
            if numel(SNR) == 1
                ssq = norm(obj.Xtrue(:),'fro')^2/(obj.I*obj.J*obj.K*(10^(SNR/10)));
                
            elseif numel(SNR) == obj.K
                
                ssq = zeros(1,obj.K);
                for k = 1:obj.K
                   ssq(k) =  norm(obj.Xtrue(:,:,k),'fro')^2/(obj.I*obj.J*(10^(SNR(k)/10)));
                end
                
            else
                disp('SNR need exactly K elements')
            end
            
        end
    end
end