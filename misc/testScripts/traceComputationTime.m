
n = 100:100:2500;
K=5;

t1 = zeros(numel(n),K);
t2 = zeros(numel(n),K);

for i = 1:numel(n)
    for k = 1:K
        A = rand(n(i));
        X = rand(n(i));
        
        f1 = @() trace(A*X);
        f2 = @() sum(sum(A'.*X));
        
        t1(i,k)=timeit(f1);
        t2(i,k)=timeit(f2);
        disp(n(i))
    end
end

%%
subplot(1,2,1)
plot(mean(t1,2))
subplot(1,2,2)
plot(mean(t2,2))
