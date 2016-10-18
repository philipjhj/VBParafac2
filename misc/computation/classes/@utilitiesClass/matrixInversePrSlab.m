function invM = matrixInversePrSlab(obj,M)




method = obj.opts.matrixProductPrSlab;

switch method
    case 'gpu'
        invM = pagefun(@inv,M);
        
    case 'naive'
        invM = zeros(size(M));
        
        for k = 1:size(M,3)
           invM(:,:,k) = inv(M(:,:,k)); 
        end
    otherwise
        invM = zeros(size(M));
        
        for k = 1:size(M,3)
           invM(:,:,k) = inv(M(:,:,k)); 
        end
end






end