function C = matrixProductPrSlab(obj,A,B)
% Computes the matrix product C = A*B with the defined method


% if obj.opts.debugFlag && obj.checkOptionSet
%     disp('Please set options first, terminating computations')
%     return
% end    

method = obj.opts.matrixProductPrSlab;

switch method
    case 'gpu'
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