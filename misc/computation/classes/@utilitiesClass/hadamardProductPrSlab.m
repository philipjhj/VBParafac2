function C = hadamardProductPrSlab(obj,A,B)
if iscell(A)
   C = bsxfun(@times,A{:});
else
    C = bsxfun(@times,A,B);
end
end