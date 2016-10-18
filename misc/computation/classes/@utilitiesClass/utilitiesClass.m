classdef utilitiesClass < handle
    
   properties
       % opt_FunctionName // options
       opts
   end
   
   
   methods
       
       %Constructor
       function obj = utilitiesClass(modelobj)
           obj.opts = modelobj.opts;
       end
       
       % Error handling
       
       function status = checkOptionSet(obj)
           % Returns;
           % 0 if option is set
           % -1 if option is not set
           
           status = 0;
           
           callerFunc = dbstack;
           callerFunc = callerFunc(2).name;
           
           optValue = obj.opt.(callerFunc);
           
           if isempty(optValue)
              disp(['Option not set for function; ',(callerFunc)]);
              status = -1;
           end
       end
       
       % Utility functions
       C = matrixProductPrSlab(obj,A,B)
       invM = matrixInversePrSlab(obj,M)
       diagM = matrixDiagonalPrSlab(obj,M)
       
       
       % Display functions
       
       function displayOptions(obj)
          
           % TODO: write method for printing options
           
       end
   end
    
end