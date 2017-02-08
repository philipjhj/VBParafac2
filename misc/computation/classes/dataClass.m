classdef dataClass < handle
    properties
        % Data
        Xunfolded
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
        R
    end
    events
       dataUpdated
    end
    
    methods
        function obj = dataClass()
            obj.ListenerHandle = addlistener(obj,'dataUpdated',@obj.setDimensions);
        end
        
        function setDimensions(obj,~,~)
            obj.n_dims = ndims(obj.Xunfolded);
            obj.R = obj.n_dims-2;
            
            obj.I = zeros(1,obj.R);
            for r = 1:obj.R
                obj.I(r) = size(obj.Xunfolded,r);
            end
            
            obj.J = size(obj.Xunfolded,obj.n_dims-1);
            obj.K = size(obj.Xunfolded,obj.n_dims);
            
            obj.X = reshape(obj.Xunfolded,[prod(obj.I) obj.J obj.K]);
        end
        
        function set.Xunfolded(obj,value)
            obj.Xunfolded = value;
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