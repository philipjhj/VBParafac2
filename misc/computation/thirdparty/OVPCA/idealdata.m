function D=idealdata(var,disp);
% function returns data with ideal normal noise

N	=100;
p	=10;

randn('seed',0);
E	=randn(p,N);
mE	=mean(E')';

E	=E-mE*ones(1,N);	% remove mean value
% ch	=chol(E*E');
% E	=inv(ch)'*E*N;

if nargin<2
	disp=0;
end;


% eval
A	=[1 0 0;
	1 0 0;
	1 1 0;
	1 0 1;
	0 1 1;
	0 0 1;
	0 1 1;
	0 0 1;
	0 1 0;
	0 0 0];

X	=[sin(0:2*pi/(N-1):2*pi);
	1:-1/(N-1):0;
	log(1:N)];
X	= X-mean(X,2)*ones(1,N); % to keep zero mean!
D	= A*X + var*E;

if disp
	for i=1:3
		subplot(5,2,2*i-1);
		bar(A(:,i));
		title(sprintf('simulated values: $\\bm{a}_{%d}$',i));
%		xlabel('dimension');
		set(gca,'XLim',[0 10]);

		subplot(5,2,2*i);
		plot(X(i,:)');
		title(sprintf('simulated values: $\\bm{x}_{%d}$',i));
%		xlabel('realization');
	end
	subplot(3,1,3);
	plot(D');
	title('Example of data realization');
end

% save idealdata D A X E var
