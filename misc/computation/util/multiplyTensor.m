function tensorNew = multiplyTensor(tensor,vec)

tensorNew = zeros(size(tensor));
K = size(tensor,3);
for k = 1:K
    tensorNew(:,:,k) = tensor(:,:,k)*vec(k);
end
    
end