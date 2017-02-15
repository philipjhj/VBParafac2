function B = transformToTensor(obj,A)
dimensionality = sum(1<size(A));
if dimensionality == 1
    B = reshape(A,1,1,numel(A));
else
    B = A;
end
end