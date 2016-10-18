function diagM = matrixDiagonalPrSlab(obj,M)

if ndims(M) == 2
   %2D (Every column to be diag of a slab) 
   diagM = arrayfun(@(idx) diag(M(:,idx)),1:size(M,2),'UniformOutput',0);
   diagM = cat(3,diagM{:});
elseif ndim(M) == 3
    %3D (Retrieve diagonal of each slab and put into diag of slab pr.
    %column)
%     warning('DiagonalPrSlab not implemented for 3D')
    diagM = arrayfun(@(idx) diag(M(:,:,idx)),1:size(M,3),'UniformOutput',0);
    diagM = cat(2,diagM{:});
end

end