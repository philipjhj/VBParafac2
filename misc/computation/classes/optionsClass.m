classdef optionsClass < handle
    properties
        % Metadata
        debugFlag = 1; % debug level: 1. CAVI; 2. Updates;
        verbose = 1; % 1, display, 0 hide everything
        tol = 1e-7;
        showIter = 500;
        maxiter
        maxTime
        
        % qDist options (with default options)
        estimationP='parafac2svd'; % Method used to approximate E(qP)
        activeParams = {'qP','qF','qC','qA','qSigma','qAlpha'};
        
        % Utilities options
        matrixProductPrSlab = 'mtimesx'; % mtimesx, mmx, gpu
        
    end
end