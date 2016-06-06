function [mA,mL,mX,om,hist]=opca_fast(D,approx,vr,ir);
% function [mA,mL,mX,ms,hist]=opca_fast(D,);
%
% function returns variational values of decomposed components
% using approximations of 0F1 from hyperg.m
%
% approx 1       Mardia 1977
%        2       Bessel/Bessel
%        3       Woods 2002 (calibrated) [default]
%        4       Woods 2002 (non-calibrated)
%        5       Bessel/Bessel (with decreasing v)
%
% vr   vector of ranks for which a Bayesian posterior is evaluated
%      (e.g. vr=[3:10])
% ir   initial guess of rank,


% start with svd
[p,N]	=size(D);
if p>N
	error('Perhaps: opca(D'')!');
end
if nargin<2
	approx	= 3;
end
if size(approx,2)<2, approx=approx*[1 1]; end
% if no vector of ranks is specified, then consider it to be full
if nargin<4
	ir	= 0;
end
if nargin<3
	vr	= 1:p-1;
end

if nargout>4,	HIST	= 1; else HIST	= 0; end;

% MEAN VALUE
if 0
	m	= mean(D,2);
	D	= D-m*ones(1,N);
	hist.m	= m;
end
%

DD	= D*D';

% scaling
D	= D/sqrt(trace(DD));
DD	= D*D';
lamu= trace(DD)

% init
tic;				%% TIME started
[X0,L0,A0]	= svd(D',0);
hist.t.svd	= toc;
tic;

l0	= diag(L0);

for r = vr
	% Step 1, ML estimate
	lD	= l0(1:r);
	Ka	= ones(r,1);
	Kx	= ones(r,1);
	ml	= lD;

	if ir>0
		l0i	= l0(ir:end);
	else
		l0i	= l0(r+1:end);
	end
	om0	= N*(p-r)/(l0i'*l0i);
	lu	= sqrt(lamu./(1:r)');
	om	= om0;

	if HIST
		hist.ml	= l0(1:r);
	end
	hist.om	= om0;

	% Step 2, iterate
	i	= 2;
	diffom	= 1;
	while (diffom~=0 & i<500)	% converge really fast, ~= is affordable

		Ka	= hyperg(p/2, om*Kx.*lD.*ml,approx(1));
		Kx	= hyperg(N/2, om*Ka.*lD.*ml,approx(2));
		m	= Ka.*Kx.*lD;
		[ml,mvar]	= momtrun(m,sqrt(1/om),  0,lu);
		om	= (p*N)/(l0'*l0 - 2*sum(Kx.*ml.*Ka.*lD) + sum(mvar));

		if HIST
			hist.ml(:,i)	= ml;
		end
		hist.om(i)	= om;
		diffom	= om-hist.om(i-1);
		i	= i+1;
	end
	r

	% save values for further use
	B.Ka{r}	= Ka;
	B.Kx{r}	= Kx;
	B.ml{r} = ml;
	B.mvar{r}	= mvar;
	B.om{r}	= om;
	B.it{r}	= i;
	B.lu{r}	= lu;
	% Bayes rank sel
	if length(vr)>1
		% additional info
		% values of: Ka, Kx, ml, om, are assumed to be fixes
		[poma,VKa,lFa]	= hyperg(p/2, om*Kx.*lD.*ml,3);
		[pomx,VKx,lFx]	= hyperg(N/2, om*Ka.*lD.*ml,3);
		[poml,pomvar,nk]	= momtrun(m,sqrt(1/om),  0,lu);
		b	= N*p/2/om;

		lnB(r)	= (-r/2*log(pi)+gammaln(r)+gammaln(r/2+1)+r*log(2)...
				  +lFa +lFx -2*om*sum(lD.*Ka.*Kx.*ml)...
				  +sum(log(nk)) + 1/2*(m'*m - 2*ml'*m + sum(mvar))...
				  -(p*N/2)*log(b));
		lnBdetail	= [-r/2*log(lamu)+1/2*gammaln(r),...
				  +lFa+lFx-2*om*sum(lD.*Ka.*Kx.*ml),...
				  +sum(log(nk)), + 1/2*(m'*m - 2*ml'*m + sum(mvar)),...
				  -(p*N/2)*log(b)];
		if r==1, lnBdet	= lnBdetail; else lnBdet(r,:)=lnBdetail; end;
		if HIST
			B.VKa{r}	= VKa;
			B.VKx{r}	= VKx;
		end
	end

end % of ranks




if length(vr)>1
	% rank selection
	exB	= exp(lnB - max(lnB));
	fr	= exB/sum(exB);

	% results
	[pr,r]	= max(fr);
else
	[poma,VKa,lFa]	= hyperg(p/2, om*Kx.*lD.*ml,3);
	[pomx,VKx,lFx]	= hyperg(N/2, om*Ka.*lD.*ml,3);
end


hist.t.iter	= toc;
mA	= A0(:,1:r)*diag(B.Ka{r});
mX	= diag(B.Kx{r})*X0(:,1:r)';
mL	= diag(B.ml{r});

% Second moments
h	= 0.00001;
VKa2	= (hyperg(p/2,om*Kx.*lD.*ml+h,approx(1)) -hyperg(p/2,om*Kx.*lD.*ml,approx(1))) / (h);

VKx2	= (hyperg(N/2,om*Ka.*lD.*ml+h,approx(1)) -hyperg(N/2,om*Ka.*lD.*ml,approx(1))) / (h);

if HIST

	hist.A0	= A0;
	hist.X0	= X0;
	hist.Ka	= Ka;
	hist.Kx	= Kx;
	hist.mvar	= mvar;
	hist.om(1)	= om0;	% done here! hist used as stoping rule before!

	hist.VKa	= VKa;
	hist.VKx	= VKx;
	hist.VKa2	= VKa2;
	hist.VKx2	= VKx2;
	hist.B	= B;

	if length(vr)>1
		hist.fr	= fr;
		hist.lnB	= lnB;
		hist.lnBdet	= lnBdet;
		hist.r	= r;
	end
end
