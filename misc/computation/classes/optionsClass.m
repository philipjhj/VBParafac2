classdef optionsClass < handle
    properties
        % Metadata
        interactive = 0; % if 1; user will be prompted for option values
        % TODO: implement prompts
        debugFlag = 1; % debug level: 1. CAVI; 2. Updates;
        verbose = 1; % 1, display, 0 hide everything
        tol = 1e-9;
        showIter = 500;
        maxIter = intmax;
        maxTime = realmax;
        
        
        initMethod='random'
        
        % qDist options (with default options)
        estimationARD ='max'; % max or avg
        estimationNoise ='avg'; % max or avg
        estimationP ='parafac2svd'; % Method used to approximate E(qP)
        activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'};
        nActiveComponents = 'threshold'; % hard / threshold
        
        % Utilities options
        matrixProductPrSlab = 'mtimesx'; % mtimesx, mmx, gpu
        
        % RNG
        rngInput = 'shuffle'
    end
    
    
    methods
        function set.rngInput(obj,value)
            if isa(value,'double')
                obj.rngInput = value;
            elseif isa(value,'char')
                if strcmp(value,'shuffle')
                    obj.rngInput = value;
                else
                    obj.rngInput = str2double(value);
                end
            end
        end
    end
end
