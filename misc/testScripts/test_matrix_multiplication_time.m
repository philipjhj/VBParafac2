function test_matrix_multiplication_time

for i = 4
dim = 1000;
slabs = 20;
disp(strcat('Dim; ',num2str(dim),' slabs; ',num2str(slabs)))
A = rand(dim,dim,slabs);
B = rand(dim,dim,slabs);
C = rand(dim,dim,slabs);

timeit(@()LOOP_multi(A,B,C))
disp('Loop done')

timeit(@()mtimesx_multi(A,B,C))
disp('mtimesx done')

timeit(@()mmx_multi(A,B,C))
disp('mmx done')
% 
% A = rand(dim,dim,slabs,'gpuArray');
% B = rand(dim,dim,slabs,'gpuArray');
% C = rand(dim,dim,slabs,'gpuArray');
% timeGPU = gputimeit(@()GPU_multi(A,B,C))
% disp(timeGPU)
% disp('GPU done')
end
disp('All done')

end
function D = GPU_multi(A,B,C)
	D = pagefun(@mtimes,pagefun(@mtimes,A,B),C);
end

function D = LOOP_multi(A,B,C)
	D = zeros(size(A,1),size(C,2),size(A,3));
	for k = 1:size(A,3)
		D(:,:,k) = A(:,:,k)*B(:,:,k)*C(:,:,k);
	end
end

function D = mtimesx_multi(A,B,C)
	D = mtimesx(mtimesx(A,B),C);
end

function D = mmx_multi(A,B,C)
	D = mmx('mult',mmx('mult',A,B),C);
end