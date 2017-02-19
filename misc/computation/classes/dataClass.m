classdef dataClass < handle
    properties
        Xunfolded
        X % X folded for the first N-2 dimensions (IxJxK)
        M % # components to estimate
        
        Xtrue = [] % without noise, X above contains noise
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
        
        % Reconstruction
        Xrecon
        Xrecon_m
        
        % Diagnostics
        ELBO_chain
        fit_chain
        n_components
        n_components_hard
        evaltime
        evaltime_cpu
        iter = 1
        ELBOall
        errorIters = []
        errorIters_parameter = {};
        stopReason
        
        % Temp Diagnostics
        ELBO = realmin;
        ELBO_prev = 0;
        ELBO_diff
        
        ticCAVIwall
        ticCAVIcpu
        startWallTime
        startCpuTime
        
        
        partitionName
        ListenerDataUpdates
        ListenerELBOUpdates
    end
    
    properties
        I
        J
        K
        n_dims
        R
    end
    events
        dataUpdated
        ELBOUpdated
    end
    
    methods
        function obj = dataClass()
            obj.ListenerDataUpdates = addlistener(obj,'dataUpdated',@obj.setDimensions);
            obj.ListenerELBOUpdates = addlistener(obj,'ELBOUpdated',@obj.updateELBOvalues);
        end
        
        function set.ELBO(obj,ELBO_new)
            ELBOupdate = newValueELBOevent(ELBO_new);
            notify(obj,'ELBOUpdated',ELBOupdate)
            obj.ELBO = ELBO_new;
        end
        function updateELBOvalues(obj,~, ELBOupdate)
            obj.ELBO_prev = obj.ELBO;
            obj.ELBO_diff = ELBOupdate.newValue-obj.ELBO_prev;
        end
        
        function set.Xunfolded(obj,value)
            obj.Xunfolded = value;
            notify(obj,'dataUpdated');
        end
        function setDimensions(obj,~,~)
            obj.n_dims = ndims(obj.Xunfolded);
            if obj.n_dims > 2
                obj.R = obj.n_dims-2;
                obj.I = zeros(1,obj.R);
                for r = 1:obj.R
                    obj.I(r) = size(obj.Xunfolded,r);
                end
                obj.J = size(obj.Xunfolded,obj.n_dims-1);
                obj.K = size(obj.Xunfolded,obj.n_dims);
            else
                obj.I = size(obj.Xunfolded,1);
                obj.J = size(obj.Xunfolded,2);
                obj.K = size(obj.Xunfolded,3);
            end
            
            
            obj.X = reshape(obj.Xunfolded,prod(obj.I),obj.J,obj.K);
        end
        
        function set.M(obj,value)
            obj.M = value;
        end
        
        function computeStartTimes(obj)
            obj.startWallTime = obj.evaltime(obj.evaltime==max(obj.evaltime));
            obj.startCpuTime = obj.evaltime_cpu(obj.evaltime_cpu==max(obj.evaltime_cpu));
        end
        function evalTime = getLatestEvalTime(obj)
            evalTime = obj.evaltime(obj.evaltime==max(obj.evaltime));
        end
        
        function restartDataDiagnostics(obj)
            obj.ELBO_chain = [];
            obj.fit_chain = [];
            obj.n_components = [];
            obj.n_components_hard = [];
            obj.evaltime = [];
            obj.evaltime_cpu = [];
            obj.iter = 1;
            obj.ELBOall = [];
            obj.errorIters = [];
            obj.errorIters_parameter = {};
            %
            %             % Temp Diagnostics
            obj.ELBO_prev = 0;
            obj.ELBO = realmin;
            
            %             obj.ELBO_diff
            %
            %             obj.ticCAVIwall
            %             obj.ticCAVIcpu
            %             obj.startWallTime
            %             obj.startCpuTime
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