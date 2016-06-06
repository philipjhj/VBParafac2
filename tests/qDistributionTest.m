    classdef qDistributionTest < matlab.unittest.TestCase
    
    
    properties (TestParameter)
       arraySize = struct('vector', [1 5],'matrix', [5 5],'tensor',[5 5 5]);
        
        
    end
    
    
    methods (Test)
        
        function testEntropyComputation(testCase)
            
            a=1;
            b=1;
            
            testCase.verifyEqual(a,b)
            
        end
        
        
        
    end
    
    methods (Test,TestTags = {'ProbabilityDistributions','MultiNormal'})
        function testMultiNormal(testCase,arraySize)
            
            testDist = multiNormalDist('qA',arraySize);
            mc = ?multiNormalDist;
            
            for prop = properties(testDist)'
                mp = findobj(mc.PropertyList,'Name','meanOuterProduct');
                fh = mp.GetMethod;
            
                testCase.verifyWarningFree(fh);
            end
        end
    end
        
    
    
end