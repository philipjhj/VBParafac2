function C = matrixProductPrSlab(obj,A,B)
% Computes the matrix product C = A*B with the defined method


% if obj.opts.debugFlag && obj.checkOptionSet
%     disp('Please set options first, terminating computations')
%     return
% end    

callerFunc = dbstack;
callerFunc = callerFunc(1).name;
           
           


method = obj.opts.matrixProductPrSlab;

switch method
    case 'gpu'
	if ~isa(A,'gpuArray')
        disp('Input A converted to gpuArray with input:')
        disp(inputname(2))
        disp(size(A))
		disp(callerFunc)
        A = gpuArray(A);
	end
	if ~isa(B,'gpuArray')
        disp('Input B converted to gpuArray with input:')
        disp(inputname(3))
        disp(size(B))
        disp(callerFunc)
		B = gpuArray(B);
	end
        C = pagefun(@mtimes,A,B);
        % C = gather(C);
        
    case 'naive'
        K = max(size(A,3),size(B,3));
        C = zeros(size(A,1),size(B,2),K);
        for k = 1:K
            if size(A,3) > 1 && size(B,3) > 1
                C(:,:,k) = A(:,:,k)*B(:,:,k);
            elseif size(A,3) > 1
                C(:,:,k) = A(:,:,k)*B;
            else
                C(:,:,k) = A*B(:,:,k);
            end
            
        end
        
    case 'mtimesx'
        C = mtimesx(A,B);
        
    case 'mmx'
        C = mmx('mult',A,B);
    otherwise
        C = mtimesx(A,B);
end
end
