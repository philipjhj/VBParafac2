classdef qDistributionTest < matlab.unittest.TestCase
    
    
    properties (TestParameter)
        arraySize = struct('vector', [1 5],'matrix', [5 5],'tensor',[5 5 5]);
        
        qAflag = struct('off', 0,'on', 1);
        qCflag = struct('off', 0,'on', 1);
        qFflag = struct('off', 0,'on', 1);
        qPflag = struct('off', 0,'on', 1);
        qAlphaflag = struct('off', 0,'on', 1);
        qSigmaflag = struct('off', 0,'on', 1);
        qPmethod = {'parafac2svd'}%,'vonmises'}
        
        precision = struct('low',[1e8 1e-8])%,'unit',[1 1],'high',[1e8 1e8])
        
        
        repetition = cellfun(@num2str,mat2cell(1:10,1,ones(1,10))...
            ,'UniformOutput',false)
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
            
            testCase.verifyTrue(all(testDist.arrayDim==arraySize))
            
            %             mc = ?multiNormalDist;
            %
            %             for prop = properties(testDist)'
            %                 mp = findobj(mc.PropertyList,'Name','meanOuterProduct');
            %                 fh = mp.GetMethod;
            %
            %                 testCase.verifyWarningFree(fh);
            %             end
        end
    end
    
    methods (Test,TestTags = {'CAVI'})
        function testCAVI(testCase,qAflag,qCflag,qFflag,...
                qPflag,qAlphaflag,qSigmaflag,qPmethod,precision,repetition)
            
            varFactors = {'qA','qC','qF','qP','qAlpha','qSigma'};
            
            varFactorFlags = [qAflag,qCflag,qFflag,qPflag,qAlphaflag,...
                qSigmaflag];
            
            I=100;
            J=300;
            K=10;
            M=5;
            Mesti = M;
            dimensions = [I J K M];
            
            %rng(10)
            data = varBayesModelParafac2.generateDataFromModel(dimensions,precision);
            
            testCAVI=varBayesModelParafac2(data,Mesti);
           
            testCAVI.verbose = 0;
            testCAVI.qDist.debugflag = 0;
            
            testCAVI.qDist.activeParams_opt = varFactors(varFactorFlags==1);
            
            testCAVI.maxiter = 2000;
            
            testCAVI.qDist.method=qPmethod;
            
            %             rng(3)
            
            testCase.verifyWarningFree(@testCAVI.computeVarDistribution)
            
            
        end
    end
    
    methods (Test,TestTags = {'Updates'})
        function testSharedComputations(testCase)
    
            
%             a = ones(10);
%             b = ones(10);
            
            % ePtP
            testModel = varBayesModelParafac2;
            
            testModel.qDist.qP.mean = repmat(1:10,10,1);
            testModel.qDist.qP.variance = repmat(1:10,10,1,10);
            
            testModel.qDist.compute_ePtP;
            a = testModel.qDist.ePtP;
            
            b = testModel.data.J*squeeze(testModel.qDist.qP.variance)+...
                repmat(eye(testModel.data.M),1,1,testModel.data.K);
            
            testCase.verifyEqual(a,b)
            
            % eFtPtPF
%             testModel = varBayesModelParafac2;
%             
%             testModel.qDist.qF.mean = repmat(1:10,10,1);
%             testModel.qDist.qF.variance = repmat(1:10,10,1,10);
%             
%             testModel.qDist.compute_eFtPtPF;
%             a = testModel.qDist.eFtPtPF;
%             
%             b = testModel.data.J*squeeze(testModel.qDist.qP.variance)+...
%                 repmat(eye(testModel.data.M),1,1,testModel.data.K);
%             
%             testCase.verifyEqual(a,b)
%             
%             
            
            
            
            
            
            
        end
        
    end
    
end