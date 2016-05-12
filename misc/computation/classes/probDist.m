classdef probDist < handle
    
    properties
        mean
        meanLog
        variance
        type
    end
    
    properties (Dependent)
        entropy
        expectedInnerProduct
        expectedOuterProduct
    end
    
    
    methods
        function probDist(obj)
           
            
        end
        
       computeEntropy(obj)
       computeMean
       computeMeanLog
       computeVariance
       computeExpectedOuterProduct
       computeExpectedInnerProduct
       
   
    end 
end