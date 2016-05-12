clear;
load('/media/data/Dropbox/Studie/Kandidat/ThesisPhilipJHJ/code/other/02581 exercise/Matlab/claus.mat')
%%
m=3;
%%
[Pm,Fm,Cm,Am,fitm] = Parafac2(X,m);
%%

[A,H,C,P,fit] = parafac2(permute(X,[2 1 3]),m);
%%
i=1;
while ishandle(i)
    figure(i);
    clf;
    i=i+1;
end

%%

plotCP(mean(tmult(cat(3,P{:}),H,2),3),A,C,EmAx,ExAx,y,2)
%%
plotCP(mean(tmult(Pm,Fm,2),3),Am,Cm,EmAx,ExAx,y,3)

%%
M = size(H,2);
K = size(X,3);
Fkm = zeros(size(X,1),M,K);
Fk = Fkm;
dFk = 0;
for k = 1:K
    Fkm(:,:,k) = Pm(:,:,k)*Fm;
    Fk(:,:,k) = P{k}*H;
    dFk = dFk + sum(sum(Fk(:,:,k)-Fkm(:,:,k)));
end