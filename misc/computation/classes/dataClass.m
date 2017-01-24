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
        
        % Listener
        ListenerHandle
    end
    
    properties
        % Dims
        I
        J
        K
        n_dims
        F
    end
    events
       dataUpdated
    end
    
    methods
        function obj = dataClass()
            obj.ListenerHandle = addlistener(obj,'dataUpdated',@obj.setDimensions);
        end
        
        function setDimensions(obj,~,~)
            obj.n_dims = ndims(obj.X);
            obj.F = obj.n_dims-2;
            
            obj.I = zeros(1,obj.F);
            for f = 1:obj.F
                obj.I(f) = size(obj.X,f);
            end
            
            obj.J = size(obj.X,obj.n_dims-1);
            obj.K = size(obj.X,obj.n_dims);
        end
        
        function set.X(obj,value)
            obj.X = value;
            notify(obj,'dataUpdated');
        end
        
        function set.M(obj,value)
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