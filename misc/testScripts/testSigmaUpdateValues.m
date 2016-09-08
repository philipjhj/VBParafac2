squeeze(sum(sum(...
                bsxfun(@times,...
                mtimesx(obj.qA.mean,mtimesx(obj.eD,mtimesx(obj.qF.mean',permute(obj.qP.mean,[2 1 3]))))...
                ,obj.data.X),1),2))');
            
%%


size(obj.qA.mean)
size(obj.eD)
size(obj.qF.mean')
size(obj.qP.mean)
%%
K = size(obj.eD,3);
I = size(obj.qA.mean,1);

mysum = zeros(1,K);
for k = 1:K
    for i = 1:I
        mysum(k) = mysum(k) + obj.qA.mean(i,:)*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(i,:,k)';
    end
end
mysum
%%


obj.XInnerProduct

%%
mysum = zeros(1,K);
for k = 1:K
    for i = 1:I
        mysum(k) = mysum(k) + obj.data.X(i,:,k)*obj.data.X(i,:,k)';
    end
end
mysum

%%

sum(obj.eAiDFtPtPFDAi,1)
%             size(obj.eAiDFtPtPFDAi)
            
            

%%

K = size(obj.eD,3);
I = size(obj.qA.mean,1);

mysum = zeros(1,K);
for k = 1:K
    for i = 1:I
        mysum(k) = mysum(k) + obj.qSigma.mean(k)*obj.qA.mean(i,:)*obj.eD(:,:,k)*obj.qF.mean'*obj.qP.mean(:,:,k)'*obj.data.X(i,:,k)';
    end
end
sum(mysum)

