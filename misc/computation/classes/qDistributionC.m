classdef qDistributionC
    properties
        qA
        qC
        qF
        qP
        qSigma
        qAlpha
        
        % Dimensions
        I % 1st Dim
        J % 2nd Dim
        K % 3rd Dim
        M % Latent Dim
        
    end
    
    methods
        function obj = qDistributionC(dims)
            obj.I = dims(1);
            obj.J = dims(2);
            obj.K = dims(3);
            obj.M = dims(4);
            
            obj.qA = nan(1,2,obj.I);
            obj.qC = nan(1,2,obj.K,obj.M);
            obj.qF = nan(1,2,obj.K,obj.M);
            obj.qP = nan(1,2,obj.K,obj.J);
            obj.qSigma = nan(1,2,obj.K);
            obj.qAlpha = nan(1,2,obj.M); 
            
            
        end
        
        function obj = compute_Pmoments(obj)
           % Superclass function to allow specialization in subclasses 
        end
        
    end
    
end